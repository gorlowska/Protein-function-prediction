"""Pobieranie i wczytywanie danych Team 2. Bez treningu i podziału danych."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
import os
import platform
import time

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

SEED = 42
ORGANISMS = [
    {"species": "Candida albicans", "taxon_id": 5476, "go_code": "CANAL"},
    {"species": "Candida parapsilosis", "taxon_id": 5480, "go_code": "CANPC"},
    {"species": "Saccharomyces cerevisiae", "taxon_id": 4932, "go_code": "YEAST"},
]
EVIDENCE_CODES = ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP", "TAS", "IC"]
GAF_COLUMNS = ["db", "protein_id", "symbol", "relation", "go_id", "reference", "evidence_code", "with_from", "aspect", "name", "synonym", "object_type", "taxon", "annotation_date", "assigned_by", "extension", "form_id"]

def session():
    s = requests.Session()
    retry = Retry(total=5, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
    s.mount("https://", HTTPAdapter(max_retries=retry))
    s.headers["User-Agent"] = "BioInformatics-Team2-course-project/0.1"
    return s

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def save_metadata(path, meta):
    meta.update({"downloaded_at_utc": datetime.now(timezone.utc).isoformat(), "sha256": sha256(path), "bytes": path.stat().st_size, "python": platform.python_version()})
    sidecar = Path(str(path) + ".metadata.json")
    temporary = Path(str(sidecar) + ".part")
    temporary.write_text(json.dumps(meta, indent=2, ensure_ascii=False), encoding="utf-8")
    os.replace(temporary, sidecar)

def cached(path):
    sidecar = Path(str(path) + ".metadata.json")
    if not path.exists():
        return False
    if not sidecar.exists():
        raise ValueError(f"Plik bez metadanych: {path}. Przenieś go poza data/raw i ponów pobranie.")
    meta = json.loads(sidecar.read_text(encoding="utf-8"))
    if meta["sha256"] != sha256(path):
        raise ValueError(f"Niezgodna suma SHA256: {path}")
    return True

def download_gaf(root, organism):
    path = Path(root) / "data/raw/gaf" / (organism["go_code"] + "-uniprot.gaf.gz")
    if cached(path):
        return path
    url = "https://current.geneontology.org/annotations/gaf/" + path.name
    temporary = Path(str(path) + ".part")
    with session() as s, s.get(url, stream=True, timeout=(30, 180)) as response:
        response.raise_for_status()
        with temporary.open("wb") as f:
            for chunk in response.iter_content(1024 * 1024):
                f.write(chunk)
        header = []
        annotations = 0
        with gzip.open(temporary, "rt", encoding="utf-8") as f:
            for line in f:
                if line.startswith("!"):
                    header.append(line.rstrip())
                elif line.strip():
                    if len(line.rstrip("\n\r").split("\t")) != 17:
                        raise ValueError(f"Nieprawidłowy format GAF: {url}")
                    annotations += 1
        if not annotations:
            raise ValueError(f"Pusty GAF: {url}")
        os.replace(temporary, path)
        save_metadata(path, {"source_url": url, "resolved_url": response.url, "last_modified": response.headers.get("Last-Modified"), "gaf_header": header, "annotation_rows": annotations, "organism": organism})
    return path

def download_uniprot(root, organism):
    path = Path(root) / "data/raw/uniprot" / (organism["go_code"] + "-reviewed.tsv")
    if cached(path):
        return path
    query = f'(reviewed:true) AND (taxonomy_id:{organism["taxon_id"]})'
    params = {"query": query, "format": "tsv", "fields": "accession,organism_name,organism_id,go_id,sequence", "size": 500}
    temporary = Path(str(path) + ".part")
    url = "https://rest.uniprot.org/uniprotkb/search"
    pages = []
    rows = 0
    release = None
    columns = None
    with session() as s, temporary.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        while url:
            response = s.get(url, params=params, timeout=(30, 180))
            response.raise_for_status()
            current_release = response.headers.get("X-UniProt-Release")
            if pages and current_release != release:
                raise ValueError("Wydanie UniProt zmieniło się podczas pobierania. Ponów pobranie.")
            release = current_release
            reader = csv.reader(response.text.splitlines(), delimiter="\t")
            header = next(reader)
            if columns is None:
                columns = header
                writer.writerow(header)
            elif header != columns:
                raise ValueError("Niezgodne kolumny kolejnych stron UniProt")
            for row in reader:
                if len(row) != len(columns):
                    raise ValueError("Nieprawidłowy rekord UniProt")
                writer.writerow(row)
                rows += 1
            pages.append(response.url)
            total = int(response.headers["X-Total-Results"])
            url = response.links.get("next", {}).get("url")
            params = None
            if url:
                time.sleep(0.15)
    if rows != total or rows == 0:
        raise ValueError(f"Niepełne lub puste dane UniProt: {rows}/{total} dla {organism['species']}")
    os.replace(temporary, path)
    save_metadata(path, {"query": query, "page_urls": pages, "uniprot_release": release, "uniprot_release_date": response.headers.get("X-UniProt-Release-Date"), "rows": rows, "organism": organism})
    return path

def download_all(root):
    result = []
    for organism in ORGANISMS:
        print(organism["species"], flush=True)
        gaf = download_gaf(root, organism)
        uniprot = download_uniprot(root, organism)
        result.append({"species": organism["species"], "gaf_file": str(gaf), "uniprot_file": str(uniprot)})
    return pd.DataFrame(result)

def load_all(root):
    proteins, annotations = [], []
    for organism in ORGANISMS:
        code = organism["go_code"]
        p = pd.read_csv(Path(root) / f"data/raw/uniprot/{code}-reviewed.tsv", sep="\t", dtype=str, keep_default_na=False)
        p = p.rename(columns={"Entry": "protein_id", "Organism": "organism", "Organism (ID)": "taxon_id", "Gene Ontology IDs": "go_ids_uniprot", "Sequence": "sequence"})
        required = {"protein_id", "organism", "taxon_id", "go_ids_uniprot", "sequence"}
        if not required.issubset(p.columns):
            raise ValueError(f"Nieoczekiwane kolumny UniProt: {p.columns.tolist()}")
        p["species"] = organism["species"]
        p["length"] = p.sequence.str.len()
        g = pd.read_csv(Path(root) / f"data/raw/gaf/{code}-uniprot.gaf.gz", sep="\t", comment="!", names=GAF_COLUMNS, dtype=str, keep_default_na=False)
        g["species"] = organism["species"]
        proteins.append(p)
        annotations.append(g)
    return pd.concat(proteins, ignore_index=True), pd.concat(annotations, ignore_index=True)

def validate_and_summarize(proteins, annotations):
    if proteins.protein_id.duplicated().any():
        raise ValueError("Powtarzające się identyfikatory UniProt wymagają analizy")
    if proteins.sequence.eq("").any() or not proteins.sequence.str.fullmatch("[A-Z]+").all():
        raise ValueError("Brak sekwencji lub nieprawidłowe znaki")
    if not annotations.go_id.str.fullmatch(r"GO:\d{7}").all():
        raise ValueError("Nieprawidłowy identyfikator GO w GAF")
    rows = []
    for organism in ORGANISMS:
        p = proteins[proteins.species.eq(organism["species"])]
        g = annotations[annotations.species.eq(organism["species"])]
        positive = g[~g.relation.str.split("|").apply(lambda values: "NOT" in values)]
        supported = positive[positive.evidence_code.isin(EVIDENCE_CODES)]
        rows.append({"species": organism["species"], "reviewed_proteins": len(p), "with_GO_in_UniProt": int(p.go_ids_uniprot.ne("").sum()), "GAF_rows": len(g), "reviewed_with_any_GAF": int(p.protein_id.isin(g.protein_id).sum()), "reviewed_with_allowed_positive_GAF": int(p.protein_id.isin(supported.protein_id).sum()), "median_length": p.length.median()})
    return pd.DataFrame(rows)
