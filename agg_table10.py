import os, statistics as st
B = 'results/table10_q512'
mats = ['bcsstk03','bcsstk14','bodyy4','mhd4800b','nasa4704','nos2','sts4098']

def final_err(f):
    L = [l.split() for l in open(f) if l.strip()]
    h = next(i for i, l in enumerate(L) if l[0] == 'iter')
    H = L[h]
    rows = [l for l in L[h+1:] if len(l) == len(H) and l[0].isdigit()]
    return float(rows[-1][H.index('solerr_p32q')])

print(f"{'Matrix':10s} {'Seeds':>5s} {'Naive(plain mv)':>16s} {'Quire mv':>12s} {'Ratio':>8s}")
for m in mats:
    p, q = [], []
    for s in range(1, 51):
        fp = f'{B}/plainmv/{m}_seed{s}.log'
        fq = f'{B}/qmv/{m}_seed{s}.log'
        if os.path.exists(fp) and os.path.exists(fq):
            p.append(final_err(fp)); q.append(final_err(fq))
    if p:
        print(f"{m:10s} {len(p):5d} {st.mean(p):16.3e} {st.mean(q):12.3e} {st.mean(p)/st.mean(q):8.2f}")
    else:
        print(f"{m:10s}     0  (no paired seeds)")
