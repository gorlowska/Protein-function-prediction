"""Prepare project folders and install the pinned Linux x86_64 MMseqs2 build."""
from pathlib import Path
import argparse
import json
import platform
import subprocess
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.data_loading import session, cached, save_metadata, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-mmseqs', action='store_true',
                        help='Prepare folders only; sufficient for viewing saved outputs.')
    args = parser.parse_args()
    for folder in ['data/raw/gaf', 'data/raw/uniprot', 'data/interim',
                   'data/processed', 'reports', 'tools', 'docs']:
        (ROOT / folder).mkdir(parents=True, exist_ok=True)
    if args.skip_mmseqs:
        print('Project directories ready:', ROOT)
        return
    if platform.system() != 'Linux' or platform.machine().lower() not in {'x86_64', 'amd64'}:
        raise SystemExit('Full clustering setup requires Linux x86_64 (Ubuntu/WSL). '
                         'Use --skip-mmseqs to view results or use prepared data.')
    if sys.version_info < (3, 12):
        raise SystemExit('Python 3.12+ is required; the tested environment is Python 3.14.4.')
    if 'sse4_1' not in Path('/proc/cpuinfo').read_text():
        raise SystemExit('The bundled installation target requires the SSE4.1 CPU feature.')
    manifest = json.loads((ROOT / 'docs/source_manifest.json').read_text())
    expected = manifest['tools/mmseqs-linux-sse41.tar.gz.metadata.json']
    archive = ROOT / 'tools/mmseqs-linux-sse41.tar.gz'
    if not cached(archive):
        url = expected['source_url']
        part = Path(str(archive) + '.part')
        with session() as client, client.get(url, stream=True, timeout=(30, 180)) as response:
            response.raise_for_status()
            with part.open('wb') as stream:
                for block in response.iter_content(1024 * 1024):
                    stream.write(block)
            if sha256(part) != expected['sha256']:
                raise SystemExit('The upstream MMseqs2 archive differs from the recorded build. '
                                 'Obtain the original archive from the team; do not silently change versions.')
            part.replace(archive)
            save_metadata(archive, {'source_url': url, 'resolved_url': response.url})
    if sha256(archive) != expected['sha256']:
        raise SystemExit('The cached MMseqs2 archive does not match the recorded build.')
    binary = ROOT / 'tools/mmseqs/bin/mmseqs'
    if not binary.exists():
        with tarfile.open(archive) as package:
            package.extractall(ROOT / 'tools', filter='data')
    version = subprocess.check_output([str(binary), 'version'], text=True).strip()
    recorded = json.loads((ROOT / 'reports/clusters_run.json').read_text())['mmseqs_version']
    if version != recorded:
        raise SystemExit(f'MMseqs2 version mismatch: expected {recorded}, received {version}.')
    print('MMseqs2 ready:', version)
    print('Project:', ROOT)


if __name__ == '__main__':
    main()
