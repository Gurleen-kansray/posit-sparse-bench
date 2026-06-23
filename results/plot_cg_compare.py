import matplotlib.pyplot as plt
import numpy as np

matrices = {
    'bcsstk03': 'FEM structural (n=112)',
    'bcsstk14': 'FEM structural (n=1806)',
    'bcsstk38': 'FEM structural (n=8032)',
    'mhd4800b': 'MHD physics (n=4800)',
    'sts4098':  'Structural (n=4098)',
}

fig, axes = plt.subplots(2, 3, figsize=(15, 9))
axes = axes.flatten()

for idx, (mat, title) in enumerate(matrices.items()):
    ax = axes[idx]
    data = np.loadtxt(f'results/cg_compare/{mat}_cg.log', skiprows=1)
    iters = data[:, 0]
    res_d = data[:, 1]
    res_f = data[:, 2]
    res_p = data[:, 3]

    ax.semilogy(iters, res_d, 'b-',  linewidth=1.5, label='double64')
    ax.semilogy(iters, res_f, 'r--', linewidth=1.5, label='float32')
    ax.semilogy(iters, res_p, 'g-',  linewidth=1.5, label='posit32+quire')

    ax.set_title(f'{mat}\n{title}', fontsize=10)
    ax.set_xlabel('CG iteration')
    ax.set_ylabel('Residual norm')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

# hide unused subplot
axes[5].set_visible(False)

fig.suptitle('Full CG Solver: double64 vs float32 vs posit32+quire\nResidual convergence across domains', 
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('results/figures/cg_convergence_compare.png', dpi=150, bbox_inches='tight')
print("Saved: results/figures/cg_convergence_compare.png")
