#!/usr/bin/env python3
"""Explicitly bounded 6B.3 campaigns; never launches N=128 evolution."""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['baseline', 'defect', 'long'])
    p.add_argument('--executable', required=True, type=Path)
    p.add_argument('--output-dir', required=True, type=Path)
    p.add_argument('--jobs', type=int, default=2)
    a = p.parse_args()
    out = a.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    manifest = dict(mode=a.mode, executable_sha256=hashlib.sha256(a.executable.read_bytes()).hexdigest())
    manifest_path = out/'campaign_manifest.json'
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise ValueError('Campaign identity changed; choose a new output directory')
    manifest_path.write_text(json.dumps(manifest, indent=2)+'\n')
    sizes = [8, 16, 32, 64, 128] if a.mode == 'defect' else [8, 16, 32, 64] if a.mode == 'baseline' else [8, 16, 32]
    families = ['regular', 'smooth', 'jitter'] if a.mode != 'long' else ['regular', 'jitter']
    tasks = [(f, n, d) for f in families for n in sizes
             for d in ([0] if a.mode == 'defect' else [1, 2, 4])]

    def run(task):
        f, n, d = task
        path = out/f'{f}_{n}_{d}.csv'
        command = [str(a.executable.resolve()), f, str(n)]
        if d:
            command.append(str(d))
        command.append(str(path))
        start = time.monotonic()
        with path.with_suffix('.log').open('w') as log:
            subprocess.run(['/usr/bin/time', '-l', *command], check=True, stdout=log, stderr=subprocess.STDOUT,
                           env={**os.environ, 'OMP_NUM_THREADS': '1'})
        return f'{task}: {time.monotonic()-start:.2f}s'

    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for message in pool.map(run, tasks):
            print(message, flush=True)


if __name__ == '__main__':
    main()
