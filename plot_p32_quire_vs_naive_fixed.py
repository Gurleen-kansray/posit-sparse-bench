import matplotlib.pyplot as plt
import numpy as np
import re

correct_gains = {
    'bcsstk03': 6.6,
    'bcsstk14': 39,
    'bcsstk38': 56,
    'nasasrb':  79,
    'mhd4800b': 23,
    's3dkt3m2': 1535,
}

titles = {
    'bcsstk03': 'bcsstk03\n(98 iters filtered: |pAp|<1e-30)',
    'bcsstk14': 'bcsstk14',
    'bcsstk38': 'bcsstk38',
    'nasasrb':  'nasasrb',
    'mhd4800b': 'mhd4800b',
    's3dkt3m2': 's3dkt3m2',
}

matrices = ['bcsstk03', 'bcsstk14', 'bcsstk38', 'nasasrb', 'mhd4800b', 's3dkt3m2']

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.patch.set_facecolor('#1a1a2e')
axes = axes.flatten()

for idx, mat in enumerate(matrices):
    ax = axes[idx]
    ax.set_facecolor('#1a1a2e')

    logfile = f'results/ladder_logs/{mat}_ladder.log'
    
    iters, p32q_errs, p32n_errs = [], [], []
    with open(logfile) as f:
        for line in f:
            if not line.startswith('iter='): continue
            it    = int(re.search(r'iter=(\d+)', line).group(1))
            pAp_d = float(re.search(r'pAp_d=([\de+\-\.]+)', line).group(1))
            p32q  = float(re.search(r'p32q=([\de+\-\.]+)', line).group(1))
            p32n  = float(re.search(r'p32n=([\de+\-\.]+)', line).group(1))
            if abs(pAp_d) < 1e-30: continue
            p32q_errs.append(abs(p32q - pAp_d) / abs(pAp_d))
            p32n_errs.append(abs(p32n - pAp_d) / abs(pAp_d))
            iters.append(it)

    ax.semilogy(iters, p32q_errs, 'b-',  linewidth=0.8, label='posit32+quire')
    ax.semilogy(iters, p32n_errs, 'r--', linewidth=0.8, label='posit32 naive')

    gain = correct_gains[mat]
    ax.text(0.97, 0.95, f'max quire gain: {gain:,}x',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=9, color='black',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#90EE90', alpha=0.9))

    ax.set_title(titles[mat], fontsize=10, color='white')
    ax.set_xlabel('CG iteration', color='white', fontsize=9)
    ax.set_ylabel('relative error vs float64', color='white', fontsize=9)
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_edgecolor('gray')
    ax.grid(True, alpha=0.2, color='gray')
    if idx >= 4:
        ax.legend(fontsize=8, facecolor='#2a2a4e', labelcolor='white')

fig.suptitle('posit32+quire vs posit32 naive: pAp relative error over 300 CG iterations\n(lower is better)',
             fontsize=13, fontweight='bold', color='white')
plt.tight_layout()
plt.savefig('results/figures/p32_quire_vs_naive_6mat.png', dpi=150, bbox_inches='tight',
            facecolor='#1a1a2e')
print("Saved: results/figures/p32_quire_vs_naive_6mat.png")
