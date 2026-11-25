# Supplemental code and data

This package contains the C++ implementations, the trajectory data and the driver programs
needed to reproduce the accuracy and runtime tables of the paper. It is self-contained: no file
outside this directory is read.

## Requirements

- A C++17 compiler. The reported numbers were produced with g++ 13.2.0 and re-checked with
  g++ 14.2.0 on x86-64. Other conforming compilers should work but will give slightly different
  wall-clock times.
- GNU Make.
- Python 3.6 or later for the driver and aggregation scripts. Only the standard library is used.

No external libraries are required. The Fréchet evaluator and the four simplifiers are
implemented in this package.

## Directory layout

```
include/            headers
  common.h            point and trajectory types
  fretchet.h          continuous Frechet distance, 2D reference implementation
  qtcompat.h          minimal container compatibility layer used by the published baselines
  DotsException.h     exception type used by the published DOTS code
  dim/
    metrics_d.h       synchronous Euclidean distance, LISSED, projections
    frechet_d.h       continuous Frechet distance, d-dimensional
    stream_d.h        the streaming delta-simplifier
    stream_varied_d.h the re-anchoring variant
    squish_e_d.h      SQUISH-E baseline
    dots_published_d.h      DOTS baseline
    dots_published_d_impl.h out-of-line part of DOTS
src/
  run_bisect.cpp    bisection driver: reaches a target retention or a target error
  run_unified.cpp   sweep driver: evaluates a fixed parameter grid
  DotsException.cpp
scripts/
  run_experiments.py  runs the benchmark over the five datasets
  make_tables.py      aggregates the output into the tables
  smoke_test.py       short end-to-end check
data/               the five datasets, see below
results/            script output (empty on delivery)
```

## Build

```
make
```

This produces `bin/run_bisect` and `bin/run_unified` with `-O3 -std=c++17`. To use a specific
compiler, pass it on the command line:

```
make CXX=clang++
```

## Data

Each dataset is a directory of plain-text trajectories. Every line is one sample. Columns are
either coordinates followed by a timestamp, or a timestamp followed by coordinates, selected at
run time with `--layout=tlast` or `--layout=tfirst`.

| directory | trajectories | dimension | layout | files | source |
|---|---|---|---|---|---|
| `data/geolife` | 120 | 2 | `tlast` | 130 | GeoLife GPS trajectories, Microsoft Research |
| `data/rome` | 120 | 2 | `tlast` | 130 | CRAWDAD roma/taxi |
| `data/ibtracs` | 120 | 2 | `tlast` | 121 | IBTrACS best-track archive |
| `data/gas` | 60 | 16 | `tfirst` | 70 | UCI gas sensor array under dynamic mixtures |
| `data/smd` | 120 | 38 | `tfirst` | 130 | Server Machine Dataset (OmniAnomaly) |

The file counts exceed the number of trajectories used because the drivers sort the directory
and keep the first *N* files that parse, so a small margin is included to make the selection
deterministic. The gas and server recordings are continuous streams and were cut into
consecutive blocks of 2000 samples. Timestamps are seconds. The server telemetry carries no
timestamp in its source, so a clock of one sample per minute was assigned.

If a dataset is absent the scripts skip it and continue.

## Normalisation

Every trajectory is normalised once, before any algorithm runs, by dividing all coordinates by
the largest single-axis span of that trajectory:

```
s    = max over axes k of ( max_i c[i][k] - min_i c[i][k] )
c'   = ( c - min ) / s
```

One divisor is shared by all axes, so aspect ratios are preserved. Timestamps are not scaled.

## Running the experiments

```
make
python3 scripts/run_experiments.py
```

For each dataset this runs `bin/run_bisect` in retained-size mode and writes
`results/bisect_<dataset>.txt`. Each line reports, for one algorithm and one retention target,
the achieved retention, the output size, the continuous Fréchet error, and the wall-clock time
per input point.

The driver bisects on each algorithm's own parameter until the retained fraction is within 0.5%
of the target. If that tolerance cannot be met it accepts the closest result within 2%; it also
stops when the parameter interval falls below 0.05% or after 40 iterations. The search reports
the closest evaluation it found, and a target that no parameter can reach is reported with status
`UNREACHABLE` and must not be quoted.

To reproduce the runtime tables, run the driver on an otherwise idle machine.

## Producing the tables

```
python3 scripts/make_tables.py
```

This writes `results/tables.csv` and prints one block per dataset. Two rows are produced for
each retention target:

- `frechet_error` corresponds to the accuracy table.
- `runtime_us_per_point` corresponds to the runtime table.

Both come from the same runs, so the two tables always describe the same operating points.

## Checking the package

```
make test
```

runs a short end-to-end check on one two-dimensional and one high-dimensional dataset,
exercising simplification, timing, bisection and Fréchet evaluation.

## Notes

**Timing.** Wall-clock times depend on the machine, the compiler and the load. The accuracy
columns are reproducible exactly; the runtime columns are not, and should be read as ratios
rather than absolute values. Times exclude Fréchet evaluation, which is measured separately and
not charged to any algorithm.

**Fréchet tolerance.** The continuous Fréchet distance is computed by the decision procedure
of Baldus and Bringmann combined with a binary search on the distance, at an absolute tolerance
of 1e-12 and a relative tolerance of 1e-9. Reported distances are exact for the output curve they
describe, to better than nine significant figures.

**Third-party code.** `include/dim/squish_e_d.h` and `include/dim/dots_published_d.h` are
d-dimensional adaptations of the published SQUISH-E and DOTS algorithms. The eviction rule of
SQUISH-E and the admissibility test of DOTS are reproduced as specified in the original
publications and are not modified. `include/qtcompat.h` supplies a small subset of container
types so that the published sources compile without an external framework.

## Licence

The trajectory data are redistributed under the terms of their original sources, cited in the
table above. The code in `include/dim/squish_e_d.h` and `include/dim/dots_published_d.h` derives
from published reference implementations and retains their terms.
