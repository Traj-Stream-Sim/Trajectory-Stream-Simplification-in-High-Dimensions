#!/usr/bin/env python3
"""Build the accuracy and runtime tables from the stored benchmark output.

Output is written to results/tables.csv and printed in Markdown.
"""
import os
import re
import sys

ALGOS = ['stream', 'svaried', 'squishe', 'dots']
SETS = ['geolife', 'rome', 'ibtracs', 'gas', 'smd']


def read_lines(path):
    """Read a results file, tolerating UTF-16 as written by PowerShell redirection."""
    raw = open(path, 'rb').read()
    if raw[:2] in (b'\xff\xfe', b'\xfe\xff'):
        return raw.decode('utf-16', errors='replace').splitlines()
    return raw.decode('utf-8', errors='replace').splitlines()


def parse(path):
    """Return {algo: {target_percent: (achieved, out_size, frechet, us_per_point)}}."""
    data = {}
    for line in read_lines(path):
        parts = line.split()
        if len(parts) < 8 or parts[0] not in ALGOS:
            continue
        if parts[-1] != 'ok':
            continue
        nums = []
        for tok in parts[1:]:
            try:
                nums.append(float(tok.rstrip('%')))
            except ValueError:
                pass
        if len(nums) < 5:
            continue
        target, achieved, out, frechet, us = nums[0], nums[1], nums[2], nums[3], nums[4]
        data.setdefault(parts[0], {})[int(round(target))] = (
            achieved, out, frechet, us)
    return data


def main():
    rows = []
    for name in SETS:
        path = os.path.join('results', 'bisect_%s.txt' % name)
        if not os.path.isfile(path):
            continue
        rows.append((name, parse(path)))

    header = 'dataset,table,retention,' + ','.join(ALGOS)
    lines = [header]
    for name, data in rows:
        for target in (60, 40, 30, 20, 10, 5):
            err = [name, 'frechet_error', target]
            tim = [name, 'runtime_us_per_point', target]
            for a in ALGOS:
                cell = data.get(a, {}).get(target)
                err.append('%.6e' % cell[2] if cell else '')
                tim.append('%.4f' % cell[3] if cell else '')
            lines.append(','.join(str(x) for x in err))
            lines.append(','.join(str(x) for x in tim))

    os.makedirs('results', exist_ok=True)
    with open(os.path.join('results', 'tables.csv'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('wrote results/tables.csv')

    for name, data in rows:
        print()
        print(name)
        print('%-6s %-42s %s' % ('ret', 'frechet error (stream/svaried/squishe/dots)',
                                 'runtime us/point'))
        for target in (60, 40, 30, 20, 10, 5):
            cells = [data.get(a, {}).get(target) for a in ALGOS]
            e = ' '.join('%.2e' % c[2] if c else '   --   ' for c in cells)
            t = ' '.join('%.3f' % c[3] if c else '  --  ' for c in cells)
            print('%-6d %-42s %s' % (target, e, t))


if __name__ == '__main__':
    main()
