# Central paths for every notebook in this project.
# Keep this file in the project root (the folder that holds data, notebooks and project_Notes). Nothing here is tied to one computer:
# the project root is simply the folder this file sits in, so moving the whole folder to a new computer needs no path edits.
# Data folder: by default <project root>/data. To keep it somewhere else (e.g. another drive) do ONE of these on that computer:
#   set the environment variable RCB_DATA_DIR, or
#   make a file local_paths.json next to this file with   {"data_dir": "D:/rcb_data"}   (forward slashes; keep it out of git)
# Run  python rcb_paths.py  in this folder to print what the notebooks will use.

import os
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent


def _data_dir():
    env = os.environ.get('RCB_DATA_DIR')
    if env:
        return Path(env).expanduser()
    f = PROJECT_ROOT / 'local_paths.json'
    if f.exists():
        try:
            d = json.loads(f.read_text(encoding='utf-8')).get('data_dir')
            if d:
                return Path(d).expanduser()
        except Exception as e:
            print(f'rcb_paths: could not read {f} ({type(e).__name__}: {e}); using <project root>/data instead')
    return PROJECT_ROOT / 'data'


DATA_DIR = _data_dir()
NOTES_DIR = PROJECT_ROOT / 'project_Notes'
PHASE_1_FILE = NOTES_DIR / 'phase_1.txt'
PHASE_2_FILE = NOTES_DIR / 'phase_2.txt'
SCRIPTS_DIR = PROJECT_ROOT / 'test_CNN' / 'scripts'   # this folder no longer exists; only 03_B still points at it


if __name__ == '__main__':
    for name in ('PROJECT_ROOT', 'DATA_DIR', 'NOTES_DIR', 'PHASE_1_FILE', 'PHASE_2_FILE'):
        p = globals()[name]
        print(f'{name:<13}: {p}   ({"found" if p.exists() else "NOT FOUND"})')
