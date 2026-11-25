#!/usr/bin/env python3
"""End-to-end check: simplify, time, bisect and evaluate Frechet on two datasets."""
import os
import subprocess
import sys

BIN = 'bin/run_bisect' + ('.exe' if os.name == 'nt' else '')

CASES = [
    ('ibtracs', 'data/ibtracs', 'tlast', 8),
    ('smd', 'data/smd', 'tfirst', 4),
]


def main():
    if not os.path.isfile(BIN):
        sys.exit('binary not found: run "make" first')
    ok = True
    for name, directory, layout, n in CASES:
        if not os.path.isdir(directory):
            print('skip %s: missing %s' % (name, directory))
            continue
        cmd = [BIN, directory, '--layout=' + layout, '--max-files=%d' % n, '--mode=size']
        print('== %s ==' % name)
        print('   ' + ' '.join(cmd))
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            print('   FAILED with code %d' % res.returncode)
            ok = False
            continue
        rows = [l for l in res.stdout.splitlines()
                if l.split() and l.split()[-1] == 'ok']
        print('   %d completed bisections' % len(rows))
        for line in rows[:4]:
            print('   ' + line.strip())
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
