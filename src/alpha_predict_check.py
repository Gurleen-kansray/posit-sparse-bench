import re, sys

def parse_ladder(path, max_iter=30):
    with open(path) as f:
        text = f.read()
    entries = re.findall(
        r'iter=(\d+)\s+pAp_d=([-\d.e+]+)\s+p8q=[-\d.e+]+\s+p8n=[-\d.e+]+\s+p8q_e2=[-\d.e+]+\s+p8n_e2=[-\d.e+]+\s+p16q=[-\d.e+]+\s+p16n=[-\d.e+]+\s+p16q_e2=[-\d.e+]+\s+p16n_e2=[-\d.e+]+\s+p32q=([-\d.e+]+)\s+p32n=([-\d.e+]+)\s+p64q=[-\d.e+]+\s+p64n=[-\d.e+]+\s+alpha_d=([-\d.e+]+)\s+alpha_p32q=[-\d.e+]+\s+alpha_p32n=[-\d.e+]+\s+rz32q=([-\d.e+]+)\s+rz32n=([-\d.e+]+)\s+alpha_full_q=[-\d.e+]+\s+alpha_full_n=([-\d.e+]+)',
        text
    )
    rows = []
    for e in entries[:max_iter]:
        it, pAp_d, pAp_q, pAp_n, alpha_d, rz_q, rz_n, alpha_full_n = e
        pAp_d, pAp_n, alpha_d, rz_q, rz_n, alpha_full_n = map(float, (pAp_d, pAp_n, alpha_d, rz_q, rz_n, alpha_full_n))
        rel_err_pAp = (pAp_n - pAp_d) / pAp_d
        rel_err_rz = (rz_n - rz_q) / rz_q
        predicted = rel_err_rz - rel_err_pAp
        actual = (alpha_full_n - alpha_d) / alpha_d
        rows.append((int(it), rel_err_rz, rel_err_pAp, predicted, actual))
    return rows

path = sys.argv[1]
print("iter\trel_err_rz\trel_err_pAp\tpredicted\tactual")
for it, rz_err, pAp_err, pred, actual in parse_ladder(path):
    print(f"{it}\t{rz_err:.3e}\t{pAp_err:.3e}\t{pred:.3e}\t{actual:.3e}")
