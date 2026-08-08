import numpy as np
import scipy.io as sio
import glob, os, sys

def characterize(x, y):
    p = x * y
    xy = np.abs(x) @ np.abs(y)
    dot = x @ y
    cond = 2 * xy / abs(dot) if dot != 0 else np.inf
    nz = p[p != 0]
    D = np.log2(np.max(np.abs(nz)) / np.min(np.abs(nz))) if len(nz) else 0.0
    partial = np.cumsum(p)
    peak = np.max(np.abs(partial)) if len(partial) else 0.0
    u = 2**-24
    n_eff = np.sum(np.abs(p) >= u * peak) if peak > 0 else 0
    return cond, D, n_eff, len(p)

if __name__ == "__main__":
    mtx_dir = sys.argv[1] if len(sys.argv) > 1 else "data/matrices"
    for f in sorted(glob.glob(os.path.join(mtx_dir, "*.mtx"))):
        try:
            A = sio.mmread(f).tocsr()
            n = A.shape[0]
            np.random.seed(0)
            p_vec = np.random.randn(n)
            Ap = A @ p_vec
            cond, D, n_eff, n_total = characterize(p_vec, Ap)
            print(f"{os.path.basename(f):20s} cond={cond:.3e} D={D:.1f} n_eff={n_eff}/{n_total}")
        except Exception as e:
            print(f"{os.path.basename(f):20s} ERROR: {e}")
