import re
import matplotlib.pyplot as plt

def parse_ladder(path):
    iters, err_q, err_n = [], [], []
    with open(path) as f:
        for line in f:
            m = re.search(r'iter=(\d+) pAp_d=([\d.eE+-]+).*p32q=([\d.eE+-]+) p32n=([\d.eE+-]+)', line)
            if m:
                it = int(m.group(1))
                d = float(m.group(2))
                q = float(m.group(3))
                n = float(m.group(4))
                if d != 0:
                    iters.append(it)
                    err_q.append(abs(q - d) / abs(d))
                    err_n.append(abs(n - d) / abs(d))
    return iters, err_q, err_n

fig, axes = plt.subplots(1, 2, figsize=(12,5))

for ax, fname, title in zip(axes, ['bcsstk38_ladder.log', 'mhd4800b_ladder.log'],
                             ['bcsstk38', 'mhd4800b']):
    it, eq, en = parse_ladder(fname)
    ax.semilogy(it, eq, 'g--', linewidth=1.5, label='p32+quire relative error')
    ax.semilogy(it, en, 'r:', linewidth=1.5, label='p32 naive relative error')
    ax.set_title(title)
    ax.set_xlabel('CG iteration')
    ax.set_ylabel('relative error vs double64')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('results/figures/p32_relerror_compare.png', dpi=150)
print("Saved: results/figures/p32_relerror_compare.png")
