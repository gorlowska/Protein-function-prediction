"""Run the notebook end to end, preserving the committed dataset fingerprints."""
from pathlib import Path
from datetime import datetime, timezone
import json
import os
import sys
import nbformat
from nbclient import NotebookClient
from verify_dataset import verify

ROOT = Path(__file__).resolve().parents[1]
path = ROOT / 'meeting_02_data_preparation.ipynb'
before = verify()
original_manifest = json.loads((ROOT / 'reports/dataset_manifest.json').read_text())
notebook = nbformat.read(path, as_version=4)
nbformat.validate(notebook)
os.environ['JUPYTER_PATH'] = str(Path(sys.prefix) / 'share/jupyter')
NotebookClient(notebook, timeout=1800, kernel_name='bioinformatics',
               resources={'metadata': {'path': str(ROOT)}}).execute()
part = path.with_suffix('.ipynb.part')
nbformat.write(notebook, part)
part.replace(path)
after = verify()
if json.loads((ROOT / 'reports/dataset_manifest.json').read_text()) != original_manifest:
    raise ValueError('Dataset fingerprints changed during notebook execution')
run = {'status': 'passed', 'checked_at_utc': datetime.now(timezone.utc).isoformat(),
       'python': sys.version.split()[0], 'executed_code_cells': sum(c.cell_type == 'code' for c in notebook.cells),
       'dataset_unchanged': True, 'dataset_verification': after}
(ROOT / 'reports/repository_verification.json').write_text(json.dumps(run, indent=2))
print(json.dumps(run, indent=2))
