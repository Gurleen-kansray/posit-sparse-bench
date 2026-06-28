# Alpha Step Size Analysis — Quire vs Naive Accumulation

Per suggestion from Prof. James Quinlan, we measured the CG step size
alpha = r·z / p·Ap in full posit32 arithmetic (both numerator and denominator),
in addition to the primary pAp metric.

## Measurement Definition

- `pAp_d` = p^T A p in double64 (ground truth)
- `p32q` = p^T A p in posit32+quire
- `p32n` = p^T A p in posit32 naive
- `rz32q` = r·z in posit32+quire
- `rz32n` = r·z in posit32 naive
- `alpha_d` = rz_d / pAp_d (double64 reference)
- `alpha_full_q` = rz32q / p32q (full posit32+quire alpha)
- `alpha_full_n` = rz32n / p32n (full posit32 naive alpha)

All errors are max relative error over 300 CG iterations.
Filtered: iterations where |pAp_d| < 1e-30 (near-convergence regime).

## Results

| matrix | pAp gain | rz gain | full alpha gain | filtered |
|---|---|---|---|---|
| bcsstk03 | 7x | 10x | 4x | 98 |
| bcsstk14 | 39x | 38x | 31x | 0 |
| bcsstk36 | 111x | 800x | 653x | 0 |
| bcsstk37 | 93x | 484x | 280x | 0 |
| bcsstk38 | 56x | 362x | 224x | 0 |
| nasasrb | 79x | 1,319x | 432x | 0 |
| mhd4800b | 23x | 14x | 12x | 221 |
| s3dkt3m2 | 1,535x | 3,049x | 1,271x | 0 |
| s3dkq4m2 | 4,531x | 2,046x | 2,351x | 0 |
| sts4098 | 39x | 90x | 45x | 0 |
| nasa4704 | 8x | 38x | 11x | 0 |
| nos2 | 10x | 114x | 13x | 0 |
| bodyy4 | 128x | 148x | 67x | 0 |

## Key Observations

1. quire beats naive on ALL three metrics across all 13 matrices — consistent result.
2. rz gain is often higher than pAp gain (bcsstk36/37/38, nasasrb, s3dkt3m2, nos2),
   meaning r·z also benefits significantly from quire accumulation.
3. full alpha gain is generally between pAp and rz gains — the combined effect.
4. bcsstk14 shows near-identical pAp gain (39x) and rz gain (38x) — both inner
   products contribute equally.
5. mhd4800b and bcsstk03 are weaker results with high filtered counts — near-convergence
   regime where posit32 hits precision floor regardless of accumulation method.

## Implementation Note

rz and pAp are computed independently in posit32+quire and posit32 naive.
The CG solver itself runs in double64 throughout — posit computations are
instrumentation only and do not affect solver trajectory.
Raw per-iteration logs: results/ladder_logs/*_ladder.log
