# `hpc/` — SLURM harness for the posit32 CG sweep

Reproduces every table in *Where Exact Accumulation Helps Conjugate Gradient*
(CoNGA 2026) by running each experiment on 50 seeds across 13 SuiteSparse
matrices. One task is one (matrix, seed) pair: a single-core process that reads
one matrix and writes one log. Nothing is shared and nothing communicates, so
the campaign scales to whatever the scheduler gives it.

The full set is 27 job arrays, 5,909 tasks and roughly 21,000 core-hours. The
reported results cost 20,807 core-hours, about \$106 at \$0.0051/core-hour.

Written for a Slurm cluster. It was run on the University of Maine's ARCSIM
Penobscot cluster, and the defaults in `config.sh` reflect that; see
[Moving to another cluster](#moving-to-another-cluster).

## Quick start

From the repository root, not from inside `hpc/`:

```bash
export SBATCH_ACCOUNT=<your-account>     # sacctmgr show assoc where user=$USER format=account
bash hpc/go.sh
```

`go.sh` builds the binaries, verifies the quire is the standard 512 bits,
regenerates the job lists and submits every array. It is safe to re-run: a
completed log is never recomputed, so running it again after a failure fills
only the gaps.

Then disconnect. The jobs belong to Slurm, not to your shell, and closing the
terminal does not touch them. To check back:

```bash
bash hpc/progress.sh        # one bar per array class
bash hpc/progress.sh -v     # per matrix, for whatever is incomplete
```

To run the pieces separately:

```bash
bash hpc/build.sh                       # once
python3 hpc/make_joblists.py            # writes hpc/jobs/ and the manifest
DRYRUN=1 bash hpc/submit.sh             # print the sbatch commands, submit nothing
bash hpc/submit.sh                      # all 27 arrays, longest first
bash hpc/submit.sh sweep ladder         # only classes with these name prefixes
```

Before committing a cluster to 21,000 core-hours, ask it for four tasks:

```bash
mkdir -p hpc/slurm results/testjob
sbatch hpc/testjob.sbatch
cat hpc/slurm/pcgtest_*.out
```

`testjob.sbatch` has the same shape as a production array with four two-second
tasks, and documents the real run's resource profile in its header. It is also
the file to send to a cluster administrator who offers to check your request.

## Configuration

Everything site-specific lives in **`config.sh`**, which every other script
sources. It is the only file to edit when moving clusters.

| Variable | Default | What it controls |
|---|---|---|
| `HPC_DIR` | `hpc` | Where this harness sits, relative to the repo root |
| `SLURM_OUT` | `$HPC_DIR/slurm` | Slurm's per-task `.out` files |
| `PARTITION` | `epyc` | Partition to submit to |
| `HPC_ACCOUNT` | `$SBATCH_ACCOUNT` | Account for usage tracking |
| `GCC_MODULE` | `gcc/12.2.0` | Compiler module; empty skips `module load` |
| `CLUSTER` | `2.0` | Per-core speed relative to the machine the cost model was fitted on |
| `RATE` | `0.008` | Dollars per core-hour, for the estimate the manifest prints |
| `MAX_WALL_H` | `168` | Wall-clock ceiling per task; longer tasks are dropped and reported |
| `SEEDS` | `50` | Seeds per matrix |

Anything already exported wins, so a one-off override needs no edit:

```bash
PARTITION=haswell SEEDS=5 bash hpc/submit.sh
```

A permanent account belongs in `~/.bashrc`, not in this file. Slurm reads
`SBATCH_ACCOUNT` and `SBATCH_PARTITION` natively, so exporting those covers
every `sbatch` you run, not only this harness.

### Why `#SBATCH` lines are not configurable

`sbatch` reads `#SBATCH` directives as literal text before any shell runs, so a
variable there does **not** expand: `#SBATCH --output=$SLURM_OUT/%x.out`
creates a directory named `$SLURM_OUT`. The `.sbatch` files therefore carry
hardcoded defaults, and `submit.sh` overrides `--partition`, `--output`,
`--time`, `--mem`, `--array` and `--account` on the command line, which is what
takes effect for the production arrays. If you `sbatch testjob.sbatch` or
`probe.sbatch` by hand, the literals in those two files apply.

## What each file does

| File | Role |
|---|---|
| `config.sh` | Site settings. Sourced by every script |
| `go.sh` | Build, generate, submit. The one-command path |
| `build.sh` | Compiles 11 binaries into `bin/`, then asserts the quire is 512 bits |
| `make_joblists.py` | Cost model and manifest. Writes `jobs/*.txt` and `jobs/manifest.tsv` |
| `submit.sh` | One `sbatch` per manifest row. Refuses to run while `pcg_*` jobs are live |
| `sweep.sbatch` | The array runner. Reads one `matrix seed [maxiter]` line per task |
| `progress.sh` | Completion per class, counted against the job lists |
| `audit.sh` | Post-run `sacct` report: timeouts, near-limit finishes, peak `MaxRSS`, log integrity |
| `requeue.py` | Manifest covering only the missing seeds, with scaled `--time` |
| `fit_costs.sh` | Refits the cost constants from what the run actually did |
| `testjob.sbatch` | Four-task smoke test, and the file to show an administrator |
| `probe.sbatch` | Reports what a partition actually gives you: node, clock, limits, a benchmark |
| `universal-headers.tgz` | Pinned Universal headers; see [Dependencies](#dependencies) |

## Which program produces which table

`make_joblists.py` prints this map when it runs, and is the authority if the two
ever disagree.

| Class | Binary | Output | Tables |
|---|---|---|---|
| `sweep` | `cg_compare_seeded` | `results/seed_sweep_50` | 5, 6, 8, 9 |
| `ladder` | `generic_ladder_seeded` | `results/ladder_50` | 4 |
| `quiregain` | `ladder_quiregain` | `results/quiregain_50` | 4 |
| `qmv` | `cg_compare_seeded_q512_qmv` | `results/qmv_50/qmv` | 10, quire matvec arm |
| `qmvplain` | `qmv_plain` | `results/qmv_50/plain` | 10, plain matvec arm |
| `refine11` | `cg_refinement_seeded_v2` | `results/refine11_50` | 11 |
| `ref12quire` | `cg_refinement_seeded_v3_q512` | `results/refine12_quire_50` | 12, quire arm |
| `ref12plain` | `cg_refinement_noquire` | `results/refine12_plain_50` | 12, plain arm |
| `refexactb` | `cg_refinement_exactb_q512` | `results/refine_exactb_50` | 13 |
| `maxn` | `cg_compare_seeded` | `results/maxn_50` | iteration-budget check |
| `budget2k` | `cg_compare_seeded` | `results/budget2k_50` | iteration-budget check |

Each class is split into `_xl`, `_lg` and `_sm` by matrix size, which is what
lets a 30-hour task and a 30-second task carry different `--time` requests.

The last two reuse `cg_compare_seeded` with a third argument that raises the
iteration budget, `n` for `maxn` and 4000 for `budget2k`, and they run 5 and 3
seeds rather than 50. `build.sh` also compiles `cg_compare_maxn`, which the
manifest does not use; it predates the `maxiter` argument.

## Dependencies

A compiler with `-std=c++20`, meaning GCC 11 or newer, and Python 3 for the two
generators. No MPI, no threads, no GPU, no Python packages beyond the standard
library.

Posit and quire arithmetic come from the header-only
[Universal](https://github.com/stillwater-sc/universal) library, shipped here as
`universal-headers.tgz` rather than cloned. **Do not clone upstream.** These
programs use the pre-generalized quire API, `quire<nbits, es, capacity>`;
current upstream `main` replaced it with `quire<NumberType, capacity, LimbType>`,
so a fresh clone fails with

```
type/value mismatch at argument 1 ... expected a type, got 32
```

on every `quire<32,2,31>`. The tarball is a pin, not a fork: the headers are
unmodified, taken from the release installed at `/usr/local/include/universal`
on the machine that produced the reported results, which also removes "which
Universal version" as a variable. Universal is MIT licensed. `build.sh` probes
the API first and exits with that explanation rather than several hundred lines
of template errors.

To refresh the pin:

```bash
tar czf hpc/universal-headers.tgz -C /usr/local/include universal
```

## Results and logs

Logs land in `results/<program>_50/<matrix>_seed<seed>.log`, one per task, with
an `_it<N>` suffix when a non-default iteration budget is used so that a
4000-iteration run never overwrites the 2000-iteration one it is compared
against.

The full set is **977 MB in 5,909 files** and is **not in git**; `.gitignore`
excludes `*.log`. The tables are built from the aggregates in `results/csv/`,
which are tracked. Pull the logs from the cluster with `rsync`:

```bash
rsync -avPz <you>@<host>:<repo>/results/ ./results/
```

Writes are atomic. Each task writes `<log>.tmp` and moves it into place, so a
task killed mid-write leaves a `.tmp` that `progress.sh` reports rather than a
truncated log that looks finished.

One caution on mixing runs: logs from the archived earlier campaign have 33
columns and the current binaries write 19, having dropped the posit8/16/64 and
guard columns and added `solerr_ffma`. Position-based parsing will silently
produce wrong numbers if the two are combined. Regenerate rather than top up.

## When tasks die

Not with `go.sh`, which rebuilds the whole manifest and resubmits all 27 arrays.
Use the requeue path, which submits only the gaps:

```bash
bash hpc/audit.sh                   # why they died; MaxRSS against ReqMem
python3 hpc/requeue.py              # dry run
python3 hpc/requeue.py --write
bash hpc/submit.sh --requeue
```

`requeue.py` scales `--time` by each class's observed timeout rate and
reallocates the concurrency slots to minimize the slowest array's finish.

## Moving to another cluster

1. Edit `PARTITION`, `HPC_ACCOUNT` and `GCC_MODULE` in `config.sh`.
2. Check the real wall-clock ceiling with `sinfo -o '%P %l'` and set `MAX_WALL_H`.
3. Run `sbatch hpc/probe.sbatch` on each candidate partition to see what you
   actually get, then set `CLUSTER` from the benchmark it prints.
4. Run `sbatch hpc/testjob.sbatch` and confirm four logs appear.
5. `DRYRUN=1 bash hpc/submit.sh` and read the `sbatch` lines.
6. `bash hpc/go.sh`.

Two findings worth carrying over. **Memory, not cores, set the throughput.**
Peak `MaxRSS` was 207 MB, so a 16 GB request was up to 80x over and capped a
503 GB node at 31 tasks whatever its core count; 2 GB is the right request.
That is also why `epyc` beat `haswell` here despite haswell being 1.6x faster
per core, since haswell's 62 GB fits far fewer tasks. And `getrusage` inside the
binary reports about half of Slurm's figure, because the cgroup is charged for
the whole job step, so size against Slurm's number.

**Do not scale local timings by `CLUSTER` when they were measured under
contention.** A 12-way-contended laptop run is slower than one idle cluster
core, not faster. Scaling the wrong way projected 23 hours for jobs that took 5
and nearly prompted a needless limit increase.

Recalibrate the cost model after any run:

```bash
bash hpc/fit_costs.sh
```
