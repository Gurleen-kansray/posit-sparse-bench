import re
import matplotlib.pyplot as plt

def parse_ladder(path):
    iters, pAp_d, p16q, p16n = [], [], [], []
    with open(path) as f:
        for line in f:
            m = re.search(r'iter=(\d+) pAp_d=([\d.eE+-]+).*p16q=([\d.eE+-]+) p16n=([\d.eE+-]+)', line)
            if m:
                iters.append(int(m.group(1)))
                pAp_d.append(float(m.group(2)))
                p16q.append(float(m.group(3)))
                p16n.append(float(m.group(4)))
    return iters, pAp_d, p16q, p16n

fig, axes = plt.subplots(1, 2, figsize=(12,5))

for ax, fname, title in zip(axes, ['bcsstk38_ladder.log', 'mhd4800b_ladder.log'],
                             ['bcsstk38', 'mhd4800b']):
    it, d, q, n = parse_ladder(fname)
    ax.semilogy(it, d, 'k-', label='pAp_d (double64, ground truth)')
    ax.semilogy(it, q, 'g--', label='p16q (posit16+quire)')
    ax.semilogy(it, n, 'r:', label='p16n (posit16 naive)')
    ax.set_title(title)
    ax.set_xlabel('CG iteration')
    ax.set_ylabel('pAp value')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('results/figures/p16_saturation_compare.png', dpi=150)
print("Saved: results/figures/p16_saturation_compare.png")
