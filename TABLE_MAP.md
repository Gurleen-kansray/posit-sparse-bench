# Table -> source map (CoNGA'26 revision, paper numbering)
All C++ quire sources use quire<32,2,31> (512-bit). Original 482-bit logs: git tag pre-q512.

| Paper item | C++ source | Logs | Aggregation | Verified |
|---|---|---|---|---|
| Table 4 (quire gain) | none (James: ladder_seeded.jl in james_scripts/) | james_scripts/ladder_seeded_summary.csv | n/a | not reproduced from this repo |
| Table 5 (non-transfer) | src/cg_compare_seeded.cpp | results/ladder_logs/seed_sweep | python3 compute_gain_ratios.py [dir] | yes, exact |
| Table 6 (pre-conv gain) | src/cg_compare_seeded.cpp | results/ladder_logs/seed_sweep | python3 compute_preconv_gain3.py (path hardcoded) | yes, exact |
| Table 7 (TOST) | none (James's pipeline) | james_scripts/stats_equivalence.csv | n/a | not from this repo |
| Table 8 (iterations) | James's runs | unclear | unclear | open question to James |
| Table 9 (stay-below) | src/cg_compare_seeded.cpp | results/ladder_logs/seed_sweep | python3 stay_table.py (needs ./bnorm_tool, data/) | yes, exact |
| Table 10 (quire matvec) | src/cg_compare_seeded_q512_qmv.cpp | results/table10_q512 | python3 agg_table10.py | new 20-seed rerun |
| Table 11 (refinement) | src/cg_refinement_seeded_v2.cpp | results/refinement_v2 | python3 compute_refinement_gains.py | yes, exact (232 runs) |
| Table 12 (refinement baseline) | src/cg_refinement_noquire.cpp, src/cg_refinement_seeded_v3_q512.cpp | results/refinement_noquire, results/refinement_q512 | revision/ | to confirm |
| Table 13 (exact b) | src/cg_refinement_exactb_q512.cpp | results/exactb_q512 | revision/ | to confirm |
| Fig. 2 (nasa4704, n iters) | src/cg_compare_maxn.cpp | results/maxn_all | revision/ | to confirm |

Table 10 build:
    g++ -std=c++20 -O2 -I../universal/include -o cg_q512_plainmv -DPLAIN_MV src/cg_compare_seeded_q512_qmv.cpp -lm
    g++ -std=c++20 -O2 -I../universal/include -o cg_q512_qmv src/cg_compare_seeded_q512_qmv.cpp -lm
Other builds: see run_seeded_sweep.sh, run_q512.sh, run_exactb_q512.sh, run_newmat_q512.sh, run_table10.sh.
