import matplotlib.pyplot as plt
import numpy as np

matrices = ['bcsstk03','bcsstk14','bcsstk36','bcsstk37','bcsstk38',
            'nasasrb','mhd4800b','s3dkt3m2','s3dkq4m2','sts4098',
            'nasa4704','nos2','bodyy4']

diag_ratio = [1.52e6, 8.94e9, 1.74e9, 9.61e8, 9.26e12,
              2.55e5, 3.73e12, 2.52e7, 1.44e7, 6.02e7,
              1.52e5, 1.23e5, 2.45e2]

val_ratio  = [3.78e16, 5.53e36, 1.13e37, 8.81e27, 2.00e38,
              2.32e22, 5.75e20, 1.01e40, 1.62e26, 5.66e54,
              2.60e20, 2.46e5, 1.84e19]

gains      = [7, 39, 111, 93, 56, 79, 23, 1535, 4531, 39, 8, 10, 128]

colors = plt.cm.tab20(np.linspace(0, 1, len(matrices)))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor('#1a1a2e')

for ax in [ax1, ax2]:
    ax.set_facecolor('#1a1a2e')
    ax.tick_params(colors='white')
    ax.yaxis.label.set_color('white')
    ax.xaxis.label.set_color('white')
    ax.title.set_color('white')
    for spine in ax.spines.values():
        spine.set_edgecolor('gray')
    ax.grid(True, alpha=0.2, color='gray')

for i, mat in enumerate(matrices):
    ax1.scatter(diag_ratio[i], gains[i], color=colors[i], s=80, zorder=5)
    ax1.annotate(mat, (diag_ratio[i], gains[i]), fontsize=7, color='white',
                 xytext=(5,3), textcoords='offset points')

    ax2.scatter(val_ratio[i], gains[i], color=colors[i], s=80, zorder=5)
    ax2.annotate(mat, (val_ratio[i], gains[i]), fontsize=7, color='white',
                 xytext=(5,3), textcoords='offset points')

ax1.set_xscale('log'); ax1.set_yscale('log')
ax1.set_xlabel('Diagonal ratio (max/min diagonal)', fontsize=10)
ax1.set_ylabel('p32 quire gain over naive', fontsize=10)
ax1.set_title('Diagonal ratio vs quire gain', fontsize=11)

ax2.set_xscale('log'); ax2.set_yscale('log')
ax2.set_xlabel('Value range ratio (max/min |value|)', fontsize=10)
ax2.set_ylabel('p32 quire gain over naive', fontsize=10)
ax2.set_title('Value range vs quire gain', fontsize=11)

fig.suptitle('Matrix properties vs posit32 quire accuracy gain\n(no single property predicts quire benefit)',
             fontsize=12, fontweight='bold', color='white')
plt.tight_layout()
plt.savefig('results/figures/property_vs_quire_gain.png', dpi=150, bbox_inches='tight',
            facecolor='#1a1a2e')
print("Saved: results/figures/property_vs_quire_gain.png")
