import csv

CONGA_MATRICES = {"bcsstk03","bcsstk14","bcsstk36","bcsstk37","bcsstk38",
                   "mhd4800b","nasasrb","nasa4704","nos2","bodyy4",
                   "s3dkt3m2","s3dkq4m2","sts4098"}

HEADER = ["matrix","bitwidth","lambda_max","condest_cmsw","nnz_static",
          "saturated_count","total_entries","sat_fraction"]

rows = []
with open("results/csv/static_conditioning.csv") as f:
    reader = csv.reader(f)
    next(reader)
    for r in reader:
        if len(r) != 8:
            continue
        matrix, bitwidth, lambda_max, condest_cmsw, nnz_static, sat_count, total, sat_frac = r
        if matrix not in CONGA_MATRICES:
            continue
        try:
            int(nnz_static); int(sat_count); int(total); float(sat_frac)
        except ValueError:
            continue
        if condest_cmsw != "CHOL_FAIL":
            try:
                float(condest_cmsw)
            except ValueError:
                continue
        rows.append(r)

latest = {}
for r in rows:
    latest[(r[0], r[1])] = r

clean_rows = sorted(latest.values(), key=lambda r: (r[0], int(r[1])))

with open("results/csv/static_conditioning_clean.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(HEADER)
    w.writerows(clean_rows)

print(f"Kept {len(clean_rows)} rows across {len(set(r[0] for r in clean_rows))} matrices")
for m in sorted(CONGA_MATRICES):
    bw = sorted(int(r[1]) for r in clean_rows if r[0]==m)
    print(f"  {m}: bitwidths {bw}" + ("" if bw==[8,16,32,64] else "  <-- INCOMPLETE"))
