"""Meeting 2: audyt adnotacji, FASTA, MMseqs2 i podział całych grup."""
from pathlib import Path
from collections import defaultdict
import csv
import gzip
import json
import random
import subprocess
import time
import hashlib
import numpy as np
import pandas as pd
from .data_loading import (ORGANISMS, EVIDENCE_CODES, SEED, session, cached,
                           save_metadata, sha256, load_all, download_all)

def download_go_evidence(root):
    """Stronicowane JSONL: accession i GO wraz z GoEvidenceType."""
    root = Path(root)
    for org in ORGANISMS:
        path = root / f'data/raw/uniprot/{org["go_code"]}-go-evidence.jsonl.gz'
        if cached(path):
            continue
        query = f'(reviewed:true) AND (taxonomy_id:{org["taxon_id"]})'
        params = {'query':query,'format':'json','fields':'accession,go_id','size':500}
        old = json.loads((root/f'data/raw/uniprot/{org["go_code"]}-reviewed.tsv.metadata.json').read_text())
        url = 'https://rest.uniprot.org/uniprotkb/search'
        part = Path(str(path)+'.part')
        pages, count = [], 0
        with session() as s, gzip.open(part,'wt',encoding='utf-8') as f:
            while url:
                r = s.get(url,params=params,timeout=(30,180))
                r.raise_for_status()
                release = r.headers.get('X-UniProt-Release')
                if release != old['uniprot_release']:
                    raise ValueError('Inne wydanie UniProt niż TSV. Wymagane świadome odświeżenie obu źródeł.')
                results = r.json()['results']
                for item in results:
                    f.write(json.dumps(item,ensure_ascii=False)+'\n')
                count += len(results)
                pages.append(r.url)
                total = int(r.headers['X-Total-Results'])
                url = r.links.get('next',{}).get('url')
                params = None
                if url:
                    time.sleep(0.15)
        if count != total or count != old['rows']:
            raise ValueError('Niepełny lub zmieniony eksport UniProt JSON')
        part.replace(path)
        save_metadata(path,{'query':query,'page_urls':pages,'uniprot_release':release,'rows':count})
        print(f'Pobrano dowody GO: {org["species"]}, {count} białek',flush=True)

def parse_go_ids(value):
    return set(x.strip() for x in value.split(';') if x.strip())

def load_evidence(root):
    rows, accessions, raw_go = [], set(), {}
    for org in ORGANISMS:
        with gzip.open(Path(root)/f'data/raw/uniprot/{org["go_code"]}-go-evidence.jsonl.gz','rt') as f:
            for line in f:
                item = json.loads(line)
                pid = item['primaryAccession']
                if pid in accessions:
                    raise ValueError('Powtórzony accession w eksporcie JSON')
                accessions.add(pid)
                terms = set()
                for ref in item.get('uniProtKBCrossReferences',[]):
                    if ref['database'] != 'GO':
                        continue
                    terms.add(ref['id'])
                    props = {x['key']:x['value'] for x in ref.get('properties',[])}
                    code, _, provider = props.get('GoEvidenceType','').partition(':')
                    rows.append({'protein_id':pid,'go_id':ref['id'],'evidence_code':code,'provider':provider,
                                 'aspect':props.get('GoTerm','').partition(':')[0],'source':'UniProt'})
                raw_go[pid] = terms
    return pd.DataFrame(rows), accessions, raw_go

def prepare_master(root):
    root = Path(root)
    proteins, gaf = load_all(root)
    evidence, json_ids, json_go = load_evidence(root)
    assert set(proteins.protein_id) == json_ids, 'Niezgodne białka TSV i JSON'
    uni = {p.protein_id:parse_go_ids(p.go_ids_uniprot) for p in proteins.itertuples()}
    assert uni == json_go, 'Niezgodne GO TSV i JSON: nie łączymy różnych migawek'
    taxon = dict(zip(proteins.protein_id,proteins.taxon_id))
    gaf['is_not'] = gaf.relation.str.split('|',regex=False).map(lambda x:'NOT' in x)
    # Powiązanie po accession; adnotacje specyficzne dla innej izoformy nie są
    # bezwarunkowo przenoszone na sekwencję kanoniczną.
    gaf['is_reviewed_id'] = gaf.protein_id.isin(uni)
    gaf['taxon_matches'] = gaf.apply(lambda r:r['taxon'].split('|')[0] == 'taxon:'+taxon.get(r['protein_id'],''),axis=1)
    gaf['canonical_form'] = gaf.form_id.eq('') | gaf.form_id.eq(gaf.protein_id) | gaf.form_id.eq('UniProtKB:'+gaf.protein_id)
    valid = gaf[gaf.db.eq('UniProtKB') & gaf.is_reviewed_id & gaf.taxon_matches & gaf.canonical_form].copy()
    positive = valid[~valid.is_not]
    positive_sets = positive.groupby('protein_id').go_id.agg(set).to_dict()
    negative = set(zip(valid.loc[valid.is_not,'protein_id'],valid.loc[valid.is_not,'go_id']))
    comparisons = []
    for p in proteins.itertuples():
        u, g = uni[p.protein_id], positive_sets.get(p.protein_id,set())
        comparisons.append({'protein_id':p.protein_id,'species':p.species,
                            'go_ids_gaf':';'.join(sorted(g)),'go_ids_only_uniprot':';'.join(sorted(u-g)),
                            'go_ids_only_gaf':';'.join(sorted(g-u)),
                            'n_go_uniprot':len(u),'n_go_gaf':len(g),'n_shared_go':len(u&g),
                            'same_go_sets':u==g,'has_matching_gaf':p.protein_id in positive_sets})
    comparison = pd.DataFrame(comparisons)
    comparison.to_csv(root/'reports/go_source_comparison.tsv',sep='\t',index=False)
    # UniProt wyznacza dopuszczalny zbiór GO; świadectwo może pochodzić z UniProt
    # lub pasującego GAF. Samo istnienie dowodu dla innego GO nie wystarcza.
    gaf_support = positive[['protein_id','go_id','evidence_code','assigned_by','aspect']].rename(columns={'assigned_by':'provider'})
    gaf_support['source'] = 'GAF'
    support = pd.concat([evidence,gaf_support],ignore_index=True).drop_duplicates()
    support['in_uniprot'] = [go in uni[pid] for pid,go in zip(support.protein_id,support.go_id)]
    support['allowed_evidence'] = support.evidence_code.isin(EVIDENCE_CODES)
    support['not_conflict'] = [(pid,go) in negative for pid,go in zip(support.protein_id,support.go_id)]
    support['selected'] = support.in_uniprot & support.allowed_evidence & ~support.not_conflict
    support.to_csv(root/'data/interim/annotation_evidence_audit.tsv',sep='\t',index=False)
    selected = support[support.selected].copy()
    selected.to_csv(root/'data/processed/label_evidence.tsv',sep='\t',index=False)
    labels = selected.groupby('protein_id').go_id.agg(lambda x:';'.join(sorted(set(x)))).rename('go_ids')
    master = proteins.merge(comparison.drop(columns='species'),on='protein_id',validate='one_to_one').merge(labels,on='protein_id',how='left',validate='one_to_one')
    master['go_ids'] = master.go_ids.fillna('')
    master['n_labels'] = master.go_ids.map(lambda x:len(parse_go_ids(x)))
    master['sequence_sha256'] = master.sequence.map(lambda x:hashlib.sha256(x.encode()).hexdigest())
    master['eligible'] = master.n_labels.gt(0)
    master['exclusion_reason'] = np.where(master.eligible,'','no_UniProt_GO_with_allowed_positive_evidence')
    master.to_csv(root/'data/processed/master_all_reviewed.tsv',sep='\t',index=False)
    eligible = master[master.eligible].copy().sort_values('protein_id').reset_index(drop=True)
    if eligible.empty:
        raise ValueError('Brak białek po filtrowaniu')
    eligible.to_csv(root/'data/processed/master_filtered.tsv',sep='\t',index=False)
    write_fasta(eligible,root/'data/processed/proteins.fasta')
    report = master.groupby('species').agg(reviewed=('protein_id','size'),with_matching_positive_GAF=('has_matching_gaf','sum'),identical_GO_sets=('same_go_sets','sum'),retained=('eligible','sum'))
    report['excluded'] = report.reviewed-report.retained
    report.to_csv(root/'reports/filtering_summary.csv')
    diagnostics = {'label_policy':'GO from UniProt; allowed evidence per protein-GO in UniProt or matching canonical/taxon GAF; exact GAF NOT conflicts excluded',
        'evidence_codes':EVIDENCE_CODES,'gaf_rows':len(gaf),'matching_canonical_taxon_gaf_rows':len(valid),
        'gaf_isoform_rows_excluded':int((gaf.is_reviewed_id & ~gaf.canonical_form).sum()),
        'gaf_taxon_mismatch_rows_excluded':int((gaf.is_reviewed_id & ~gaf.taxon_matches).sum()),
        'not_conflicting_pairs':len(set(zip(support.loc[support.not_conflict,'protein_id'],support.loc[support.not_conflict,'go_id']))),
        'retained_proteins':len(eligible),'unique_GO_terms':len(set().union(*(parse_go_ids(x) for x in eligible.go_ids))),
        'selected_label_pairs':len(selected[['protein_id','go_id']].drop_duplicates())}
    (root/'reports/annotation_policy.json').write_text(json.dumps(diagnostics,indent=2))
    return master, eligible, report

def write_fasta(df,path):
    with Path(path).open('w') as f:
        for p in df.itertuples():
            f.write(f'>{p.protein_id}\n')
            for i in range(0,len(p.sequence),80):
                f.write(p.sequence[i:i+80]+'\n')

def run_mmseqs(root, fasta, prefix, mode='easy-cluster'):
    root = Path(root)
    binary = root/'tools/mmseqs/bin/mmseqs'
    version = subprocess.check_output([str(binary),'version'],text=True).strip()
    temp = root/'data/interim/mmseqs_tmp'/prefix
    temp.parent.mkdir(parents=True,exist_ok=True)
    output = root/'data/interim'/prefix
    args = [str(binary),mode,str(fasta)]
    if mode == 'easy-search':
        args += [str(fasta)]
        output = output.with_suffix('.tsv')
    args += [str(output),str(temp),'--min-seq-id','0.5','-c','0.8','--cov-mode','0',
             '--alignment-mode','3','--threads','2','--split-memory-limit','2G']
    if mode == 'easy-cluster':
        args += ['--cluster-reassign','1']
        expected = Path(str(output)+'_cluster.tsv')
    else:
        args += ['--format-output','query,target,fident,qcov,tcov','--max-seqs','20000','-s','7.5']
        expected = output
    stamp = root/'reports'/f'{prefix}_run.json'
    fingerprint = {'mmseqs_version':version,'command':args,'input_sha256':sha256(fasta)}
    if expected.exists() and stamp.exists():
        previous = json.loads(stamp.read_text())
        if all(previous.get(k)==v for k,v in fingerprint.items()) and previous.get('output_sha256')==sha256(expected):
            return expected
    log = root/'reports'/f'{prefix}.log'
    start = time.time()
    print('Uruchamiam:',mode,flush=True)
    # /usr/bin/time mierzy maksimum RSS procesu i jego zakończonych dzieci.
    rss = root/'reports'/f'{prefix}_resources.txt'
    with log.open('w') as f:
        result = subprocess.run(['/usr/bin/time','-v','-o',str(rss),*args],stdout=f,stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'MMseqs2 zakończony błędem. Zobacz {log}:\n'+log.read_text()[-4000:])
    assert expected.exists()
    fingerprint.update({'elapsed_seconds':time.time()-start,'output_sha256':sha256(expected)})
    stamp.write_text(json.dumps(fingerprint,indent=2))
    print(f'{mode}: {time.time()-start:.1f} s',flush=True)
    return expected

def build_groups(root,eligible):
    root = Path(root)
    fasta = root/'data/processed/proteins.fasta'
    clusters_path = run_mmseqs(root,fasta,'clusters')
    clusters = pd.read_csv(clusters_path,sep='\t',names=['mmseqs_cluster_id','protein_id'],dtype=str)
    assert not clusters.protein_id.duplicated().any()
    assert set(clusters.protein_id)==set(eligible.protein_id)
    # MMseqs2 jest heurystyczny: klaster z reprezentantem sam w sobie nie
    # gwarantuje braku podobnych par między klastrami. Szukamy dodatkowych krawędzi.
    hits_path = run_mmseqs(root,fasta,'homology_audit',mode='easy-search')
    parent = {p:p for p in eligible.protein_id}
    def find(x):
        while parent[x]!=x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b:
            parent[max(a,b)] = min(a,b)
    for rep,pid in clusters.itertuples(index=False,name=None):
        union(rep,pid)
    # Duplikaty sekwencji nie mogą trafić do różnych zbiorów.
    for _, group in eligible.groupby('sequence_sha256'):
        ids = sorted(group.protein_id)
        for pid in ids[1:]:
            union(ids[0],pid)
    hit_count = 0
    with hits_path.open() as f:
        for row in csv.reader(f,delimiter='\t'):
            a,b,identity,qcov,tcov = row
            if float(identity)+1e-6 >= .5 and min(float(qcov),float(tcov))+1e-6 >= .8:
                union(a,b)
                if a!=b:
                    hit_count+=1
    clusters['cluster_id'] = clusters.protein_id.map(lambda x:find(x))
    # MMseqs2 may emit identical cluster members in a different row order.
    clusters = clusters.sort_values('protein_id', kind='stable').reset_index(drop=True)
    clusters.to_csv(root/'data/processed/clusters.tsv',sep='\t',index=False)
    (root/'reports/clustering_summary.json').write_text(json.dumps({'mmseqs_clusters':int(clusters.mmseqs_cluster_id.nunique()),
        'final_groups':int(clusters.cluster_id.nunique()),'nonself_qualifying_directed_hits':hit_count,
        'policy':'Whole MMseqs clusters merged by exact duplicates and qualifying all-vs-all search hits. Search is heuristic, not exhaustive proof.'},indent=2))
    return eligible.merge(clusters,on='protein_id',validate='one_to_one')

def assign_splits(root,df,seed=SEED,trials=128):
    """Wiele losowych kolejności + zachłanny dobór całych grup; test służy tylko
    projektowaniu podziału przed treningiem. Brak doboru modelu na teście."""
    root = Path(root)
    names = ['train','eval','test']
    ratios = np.array([.8,.1,.1])
    species = sorted(df.species.unique())
    terms = sorted(set().union(*(parse_go_ids(x) for x in df.go_ids)))
    term_index = {term:i for i,term in enumerate(terms)}
    groups = []
    for gid, group in df.groupby('cluster_id',sort=True):
        sc = np.array([(group.species==s).sum() for s in species],dtype=float)
        lc = np.zeros(len(terms))
        for value in group.go_ids:
            for term in parse_go_ids(value):
                lc[term_index[term]]+=1
        groups.append((gid,len(group),sc,lc))
    sizes = np.array([g[1] for g in groups],float)
    all_sc = np.sum([g[2] for g in groups],axis=0)
    all_lc = np.sum([g[3] for g in groups],axis=0)
    # Przy ocenie rozkładu etykiet uwzględniamy tylko etykiety obecne w >=10 białkach.
    common = all_lc>=10
    targets = ratios*len(df)
    starget = ratios[:,None]*all_sc
    ltarget = ratios[:,None]*all_lc
    rng = np.random.default_rng(seed)
    best = None
    for trial in range(trials):
        # Duże grupy najpierw; jitter daje alternatywne kolejności w każdej próbie.
        order = np.argsort(-(sizes*rng.uniform(.7,1.3,len(groups))),kind='stable')
        counts=np.zeros(3); scounts=np.zeros((3,len(species))); lcounts=np.zeros((3,len(terms)))
        assignments={}
        for index in order:
            gid,size,sc,lc=groups[index]
            # Przyrost znormalizowanego błędu kwadratowego.
            delta=((counts+size-targets)**2-(counts-targets)**2)/np.maximum(targets,1)
            delta+=.25*np.sum(((scounts+sc-starget)**2-(scounts-starget)**2)/np.maximum(starget,1),axis=1)
            if common.any():
                delta+=.05*np.mean(((lcounts[:,common]+lc[common]-ltarget[:,common])**2-(lcounts[:,common]-ltarget[:,common])**2)/np.maximum(ltarget[:,common],1),axis=1)
            chosen=int(np.argmin(delta))
            counts[chosen]+=size; scounts[chosen]+=sc; lcounts[chosen]+=lc
            assignments[gid]=names[chosen]
        score=float(np.sum((counts-targets)**2/np.maximum(targets,1))+.25*np.sum((scounts-starget)**2/np.maximum(starget,1)))
        if min(counts)>0 and (best is None or score<best[0]):
            best=(score,assignments)
    if best is None:
        raise ValueError('Za mało niezależnych grup do trzech niepustych zbiorów')
    result=df.copy()
    result['split']=result.cluster_id.map(best[1])
    result=result.sort_values(['split','protein_id']).reset_index(drop=True)
    result.to_csv(root/'data/processed/master_split.tsv',sep='\t',index=False)
    (root/'data/processed/splits').mkdir(exist_ok=True)
    for name in names:
        subset=result[result.split.eq(name)]
        subset.to_csv(root/f'data/processed/splits/{name}.tsv',sep='\t',index=False)
        write_fasta(subset,root/f'data/processed/splits/{name}.fasta')
    # Słownik klas modelu wyznaczamy wyłącznie na train. Pozostałe etykiety
    # zachowujemy w plikach, a ich brak w train jawnie raportujemy.
    label_rows=[]
    for p in result.itertuples():
        for term in sorted(parse_go_ids(p.go_ids)):
            label_rows.append((p.protein_id,term,p.split))
    labels=pd.DataFrame(label_rows,columns=['protein_id','go_id','split'])
    labels.to_csv(root/'data/processed/labels.tsv',sep='\t',index=False)
    train_terms=set(labels.loc[labels.split.eq('train'),'go_id'])
    pd.DataFrame({'go_id':sorted(train_terms),'label_index':range(len(train_terms))}).to_csv(root/'data/processed/train_label_vocabulary.tsv',sep='\t',index=False)
    coverage=labels.groupby(['go_id','split']).size().unstack(fill_value=0).reindex(columns=names,fill_value=0)
    coverage['seen_in_train']=coverage.train.gt(0)
    coverage.to_csv(root/'reports/label_coverage.tsv',sep='\t')
    stats=result.groupby('split').agg(proteins=('protein_id','size'),clusters=('cluster_id','nunique'),label_pairs=('n_labels','sum')).reindex(names)
    stats['fraction']=stats.proteins/len(result)
    stats.to_csv(root/'reports/split_summary.csv')
    pd.crosstab(result.species,result.split).reindex(columns=names,fill_value=0).to_csv(root/'reports/split_by_species.csv')
    config={'seed':seed,'trials':trials,'target_ratios':dict(zip(names,ratios.tolist())),
            'objective':'size + species balance; common-label balance within greedy trials; whole groups only',
            'actual_fractions':stats.fraction.to_dict(),'unseen_GO_in_eval':len(set(labels.loc[labels.split.eq('eval'),'go_id'])-train_terms),
            'unseen_GO_in_test':len(set(labels.loc[labels.split.eq('test'),'go_id'])-train_terms)}
    (root/'reports/split_config.json').write_text(json.dumps(config,indent=2))
    return result,stats

def validate_outputs(root,df):
    root=Path(root)
    assert not df.protein_id.duplicated().any()
    assert set(df.split)=={'train','eval','test'}
    assert df.groupby('cluster_id').split.nunique().max()==1
    assert df.groupby('mmseqs_cluster_id').split.nunique().max()==1
    assert df.groupby('sequence_sha256').split.nunique().max()==1
    assert df.n_labels.gt(0).all()
    assert all(parse_go_ids(p.go_ids)<=parse_go_ids(p.go_ids_uniprot) for p in df.itertuples())
    ev=pd.read_csv(root/'data/processed/label_evidence.tsv',sep='\t',keep_default_na=False)
    assert ev.evidence_code.isin(EVIDENCE_CODES).all()
    assert not ev.not_conflict.any()
    splits=df.set_index('protein_id').split.to_dict()
    cross_hits=0
    with (root/'data/interim/homology_audit.tsv').open() as f:
        for a,b,identity,qcov,tcov in csv.reader(f,delimiter='\t'):
            if float(identity)>=.5 and min(float(qcov),float(tcov))>=.8 and splits[a]!=splits[b]:
                cross_hits+=1
    assert cross_hits==0
    for split in ['train','eval','test']:
        expected=set(df.loc[df.split.eq(split),'protein_id'])
        fasta=root/f'data/processed/splits/{split}.fasta'
        actual={line[1:].strip() for line in fasta.read_text().splitlines() if line.startswith('>')}
        assert expected==actual
    metadata = {'status':'passed','proteins':len(df),'clusters':int(df.cluster_id.nunique()),'cross_split_cluster_count':0,
                'cross_split_exact_duplicate_count':0,'cross_split_qualifying_search_hits':cross_hits,
                'limitation':'MMseqs search is heuristic; zero detected cross-split hits is not an exhaustive pairwise guarantee.'}
    (root/'reports/validation.json').write_text(json.dumps(metadata,indent=2))
    manifest={str(p.relative_to(root)):sha256(p) for base in ['data/processed','data/raw'] for p in sorted((root/base).rglob('*')) if p.is_file() and p.name!='.gitkeep'}
    (root/'reports/dataset_manifest.json').write_text(json.dumps(manifest,indent=2))
    return metadata
