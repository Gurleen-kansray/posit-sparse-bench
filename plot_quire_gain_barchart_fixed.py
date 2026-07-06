import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

matrices = ['bcsstk03', 'bcsstk14', 'bcsstk38', 'nasasrb', 'mhd4800b', 's3dkt3m2']
gains    = [7, 39, 56, 79, 23, 1535]
domains  = ['FEM structural', 'FEM structural', 'FEM structural',
            'NASA structural', 'MHD physics', '3D structural FEM']
colors   = {'FEM structural': '#4C9BE8', 'NASA structural': '#3CB371',
            'MHD physics': '#FFA500', '3D structural FEM': '#9370DB'}

fig, ax = plt.subplots(figsize=(10, 6))
fig.patch.set_facecolor('#1a1a2e')
ax.set_facecolor('#1a1a2e')

bar_colors = [colors[d] for d in domains]
bars = ax.bar(matrices, gains, color=bar_colors, edgecolor='white', linewidth=0.5)

for bar, gain, domain in zip(bars, gains, domains):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.15,
            f'{gain:,}x', ha='center', va='bottom', fontsize=11, fontweight='bold', color='white')
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height()/2,
            domain.replace(' ', '\n'), ha='center', va='center', fontsize=8, color='white', fontweight='bold')

ax.set_yscale('log')
ax.set_ylabel('Quire improvement over naive accumulation\n(p32 naive error / p32+quire error)', fontsize=11, color='white')
ax.set_xlabel('Sparse matrix', fontsize=11, color='white')
ax.set_title('posit32+quire accuracy gain over posit32 naive\nacross 6 real-world sparse matrices (300 CG iterations)', fontsize=12, fontweight='bold', color='white')
ax.tick_params(colors='white')
for spine in ax.spines.values():
    spine.set_edgecolor('gray')

legend_elements = [Patch(facecolor=colors[d], label=d) for d in dict.fromkeys(domains)]
ax.legend(handles=legend_elements, fontsize=9, loc='upper left', facecolor='#2a2a4e', labelcolor='white')
ax.grid(axis='y', alpha=0.3, color='gray')
ax.set_ylim(bottom=1)

plt.tight_layout()
plt.savefig('results/figures/quire_gain_barchart.png', dpi=150, bbox_inches='tight', facecolor='#1a1a2e')
print("Saved: results/figures/quire_gain_barchart.png")
