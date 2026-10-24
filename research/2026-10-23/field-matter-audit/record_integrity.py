"""Record source, preserved input and actual output checksums; never synthesize results."""
import argparse
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

root = Path(__file__).resolve().parent
repo = root.parents[2]
parser = argparse.ArgumentParser()
parser.add_argument('--junit', type=Path, required=True)
args = parser.parse_args()
test_run = ET.parse(args.junit).getroot()
suites = test_run.findall('testsuite')
files = list(root.glob('*.py'))+list((root/'results').glob('*'))
files += [repo/'research/2026-10-19/vortex-extended/vacuum_spherical.py',
          repo/'research/2026-10-19/vortex-extended/vacuum_spherical_h0p1_L80.npz']
result = {'junit_summary': [{key: suite.attrib[key] for key in ('tests', 'failures', 'errors', 'skipped', 'time')} for suite in suites],
          'sha256': {str(p.relative_to(repo)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
                     for p in sorted(files) if p.is_file() and p.name != 'integrity.json'}}
(root/'integrity.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
