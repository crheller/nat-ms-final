"""
Cartoon schematic - two neurons, 3 stimuli
High vs. low arousal changes
    gain
    FF
    noise corr.
"""
import sys
sys.path.append("/auto/users/hellerc/code/projects/nat-ms-final/")
from path_settings import PY_FIGURES_DIR
import numpy as np
import charlieTools.plotting as cplt
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['font.size'] = 8

fig_fn = f"{PY_FIGURES_DIR}figure1.svg"

# generate "baseline" sample data
cov = np.array([[1, 0.4], [0.4, 1]])
u1 = np.array([1, 2])
u2 = np.array([2, 1])
u3 = np.array([2.4, 3.1])
k = 5000
d1 = np.random.multivariate_normal(u1, cov, k)
d2 = np.random.multivariate_normal(u2, cov, k)
d3 = np.random.multivariate_normal(u3, cov, k)
e1 = cplt.compute_ellipse(d1[:, 0], d1[:, 1])
e2 = cplt.compute_ellipse(d2[:, 0], d2[:, 1])
e3 = cplt.compute_ellipse(d3[:, 0], d3[:, 1])

# plot 4 subplots: baseline, change gain, change variance, change noise corr.
f, ax = plt.subplots(1, 4, figsize=(8, 2))

# baseline
ax[0].set_title("Low Arousal")
ax[0].plot(e1[0], e1[1], lw=2, label="Stimulus 1")
ax[0].plot(e2[0], e2[1], lw=2, label="Stimulus 2")
ax[0].plot(e3[0], e3[1], lw=2, label="Stimulus 3")
ax[0].axis("square")
ax[0].set_xticks([])
ax[0].set_yticks([])
ax[0].set_xlabel("Dim. 1")
ax[0].set_ylabel("Dim. 2")
ax[0].legend()

# change gain
g = 1.3
_u1 = u1 * g
_u2 = u2 * g
_u3 = u3 * g
_d1 = np.random.multivariate_normal(_u1, cov, k)
_d2 = np.random.multivariate_normal(_u2, cov, k)
_d3 = np.random.multivariate_normal(_u3, cov, k)
e1 = cplt.compute_ellipse(_d1[:, 0], _d1[:, 1])
e2 = cplt.compute_ellipse(_d2[:, 0], _d2[:, 1])
e3 = cplt.compute_ellipse(_d3[:, 0], _d3[:, 1])
ax[1].set_title("Modulate gain")
ax[1].plot(e1[0], e1[1], lw=2)
ax[1].plot(e2[0], e2[1], lw=2)
ax[1].plot(e3[0], e3[1], lw=2)
ax[1].axis("square")
ax[1].set_xticks([])
ax[1].set_yticks([])
ax[1].set_xlabel("Dim. 1")
ax[1].set_ylabel("Dim. 2")

# change variance
_cov = cov * 0.6
_d1 = np.random.multivariate_normal(u1, _cov, k)
_d2 = np.random.multivariate_normal(u2, _cov, k)
_d3 = np.random.multivariate_normal(u3, _cov, k)
e1 = cplt.compute_ellipse(_d1[:, 0], _d1[:, 1])
e2 = cplt.compute_ellipse(_d2[:, 0], _d2[:, 1])
e3 = cplt.compute_ellipse(_d3[:, 0], _d3[:, 1])
ax[2].set_title("Modulate variance")
ax[2].plot(e1[0], e1[1], lw=2)
ax[2].plot(e2[0], e2[1], lw=2)
ax[2].plot(e3[0], e3[1], lw=2)
ax[2].axis("square")
ax[2].set_xticks([])
ax[2].set_yticks([])
ax[2].set_xlabel("Dim. 1")
ax[2].set_ylabel("Dim. 2")

# change correlations
_d1 = d1; _d1[:, 0] = _d1[np.random.choice(range(d1.shape[0]), d1.shape[0]), 0]
_d2 = d2; _d2[:, 0] = _d2[np.random.choice(range(d1.shape[0]), d1.shape[0]), 0]
_d3 = d3; _d3[:, 0] = _d3[np.random.choice(range(d1.shape[0]), d1.shape[0]), 0]
e1 = cplt.compute_ellipse(_d1[:, 0], _d1[:, 1])
e2 = cplt.compute_ellipse(_d2[:, 0], _d2[:, 1])
e3 = cplt.compute_ellipse(_d3[:, 0], _d3[:, 1])
ax[3].set_title("Modulate correlation")
ax[3].plot(e1[0], e1[1], lw=2)
ax[3].plot(e2[0], e2[1], lw=2)
ax[3].plot(e3[0], e3[1], lw=2)
ax[3].axis("square")
ax[3].set_xticks([])
ax[3].set_yticks([])
ax[3].set_xlabel("Dim. 1")
ax[3].set_ylabel("Dim. 2")

for a in ax:
    a.set_xlim((-1, 6))
    a.set_ylim((-1, 6))

f.tight_layout()

f.savefig(fig_fn)
f.savefig(fig_fn.replace(".svg", ".png"))