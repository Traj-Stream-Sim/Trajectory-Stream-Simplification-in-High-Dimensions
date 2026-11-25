#!/usr/bin/env python3
"""Run the full benchmark over the five datasets and store the raw output.

Each dataset is processed by run_bisect in retained-size mode, which sweeps the six retention
targets and reports, for every algorithm, the achieved retention, the output size, the
continuous Frechet error and the wall-clock time per input point. The same runs produce the
accuracy table and the runtime table.
"""
import os
import subprocess
import sys

BIN = os.path.join('bin', 'run_bisect')
OUT = 'results'

DATASETS = [
    # name,      directory,   layout,   max files
    ('geolife', 'data/geolife', 'tlast', 120),
    ('rome',    'data/rome',    'tlast', 120),
    ('ibtracs', 'data/ibtracs', 'tlast', 120),
    ('gas',     'data/gas',     'tfirst', 60),
    ('smd',     'data/smd',     'tfirst', 120),
]


def main():
    if not os.path.isfile(BIN) and not os.path.isfile(BIN + '.exe'):
        sys.exit('binary not found: run "make" first')
    os.makedirs(OUT, exist_ok=True)
    exe = BIN + ('.exe' if os.name == 'nt' else '')
    for name, directory, layout, n in DATASETS:
        if not os.path.isdir(directory):
            print('skip %-8s (missing %s)' % (name, directory))
            continue
        cmd = [exe, directory, '--layout=' + layout,
               '--max-files=%d' % n, '--mode=size']
        print('run  %-8s %s' % (name, ' '.join(cmd)))
        with open(os.path.join(OUT, 'bisect_%s.txt' % name), 'w') as fh:
            rc = subprocess.call(cmd, stdout=fh)
        if rc != 0:
            print('  exited with code %d' % rc)


if __name__ == '__main__':
    main()
