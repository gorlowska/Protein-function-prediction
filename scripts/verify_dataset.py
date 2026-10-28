"""Verify the committed dataset snapshot without downloading data or running ML."""
from pathlib import Path
from collections import Counter
import csv
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with path.open(newline='', encoding='utf-8') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def fasta(path):
    records = {}
    key = None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            key = line[1:].strip()
            if key in records:
                raise ValueError(f'Duplicate FASTA identifier: {key}')
            records[key] = ''
        elif line.strip():
            if key is None:
                raise ValueError('FASTA sequence before identifier')
            records[key] += line.strip()
    return records


def verify():
    manifest = json.loads((ROOT / 'reports/dataset_manifest.json').read_text())
    for name, expected in manifest.items():
        path = ROOT / name
        if not path.is_file():
            raise ValueError(f'Missing dataset file: {name}')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f'Dataset checksum mismatch: {name}')
    table = rows(ROOT / 'data/processed/master_split.tsv')
    ids = {r['protein_id'] for r in table}
    if len(ids) != len(table):
        raise ValueError('Duplicate protein IDs in master table')
    counts = Counter(r['split'] for r in table)
    expected_counts = {'train': 5005, 'eval': 625, 'test': 626}
    if counts != expected_counts:
        raise ValueError(f'The meeting-2 snapshot changed: {counts}')
    for column in ['cluster_id', 'mmseqs_cluster_id', 'sequence_sha256']:
        groups = {}
        for record in table:
            old = groups.setdefault(record[column], record['split'])
            if old != record['split']:
                raise ValueError(f'{column} crosses split boundaries')
    for split in expected_counts:
        actual = rows(ROOT / f'data/processed/splits/{split}.tsv')
        expected = [r for r in table if r['split'] == split]
        if actual != expected:
            raise ValueError(f'TSV split differs from master: {split}')
        if fasta(ROOT / f'data/processed/splits/{split}.fasta') != {
            r['protein_id']: r['sequence'] for r in expected
        }:
            raise ValueError(f'FASTA sequences differ from TSV: {split}')
    labels = {(r['protein_id'], go) for r in table for go in r['go_ids'].split(';')}
    exported = rows(ROOT / 'data/processed/labels.tsv')
    if labels != {(r['protein_id'], r['go_id']) for r in exported} or len(labels) != len(exported):
        raise ValueError('Exported label pairs differ from master')
    split_by_id = {r['protein_id']: r['split'] for r in table}
    if any(r['split'] != split_by_id[r['protein_id']] for r in exported):
        raise ValueError('A label is assigned to the wrong split')
    training_terms = {r['go_id'] for r in exported if r['split'] == 'train'}
    vocabulary = rows(ROOT / 'data/processed/train_label_vocabulary.tsv')
    if {r['go_id'] for r in vocabulary} != training_terms:
        raise ValueError('Vocabulary does not match training labels')
    notebook = json.loads((ROOT / 'meeting_02_data_preparation.ipynb').read_text())
    code_cells = [c for c in notebook['cells'] if c['cell_type'] == 'code']
    if any(c['execution_count'] is None for c in code_cells):
        raise ValueError('The shared notebook is missing executed cells')
    if any(o['output_type'] == 'error' for c in code_cells for o in c.get('outputs', [])):
        raise ValueError('The shared notebook contains an execution error')
    return {'status': 'passed', 'verified_files': len(manifest), 'proteins': len(table),
            'label_pairs': len(labels), 'split_counts': dict(counts)}


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2))
