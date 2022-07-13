"""
Population subspaces -- Hypothesis: We have a (relatively) higher dimensional stimulus space than noise space.
Figures:
    * Stimulus space is low-D relative to N neurons, but dependent on number of neurons / stimuli
    * FA shared space is low-D relative to N neurons, but does not depend on neurons / stimuli (lrg / small pupil separate to test hyp. / follow Yu paper)
    * FA shared space is smaller than stimulus space for most recording sites
    * Do / do not overlap - more / less likely to overlap depending on stim dimensionality? More / less overlap in big / small pupil?

Supplementals:
    * Stimulus space dimensionality ~consistent between large / small pupil so just use overall dim.
    * FA pop metric changes are a bit variable, but all tend to correspond to decreased noise correlation -- error bar plots or scatter?
        * Actually, possibly this is main text figure -- significant change in loading similarity across both 322 / 331.
    - 
""" 
import scipy.stats as ss
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['font.size'] = 12

import sys
sys.path.append("/auto/users/hellerc/code/projects/nat-ms-final/")
from path_settings import PY_FIGURES_DIR
from global_settings import CPN_SITES, HIGHR_SITES
import colors
from FactorAnalysis.loader import load_pop_metrics

from nems_lbhb.analysis.pop_models import subspace_overlap

savefig = False
fig_fn = PY_FIGURES_DIR + "figure5.svg"
fig_fnS1 = PY_FIGURES_DIR + "figure5_S1.svg"

modelname = "factor_analysis_pca_evoked"
#modelname = "factor_analysis_pca"
sites = CPN_SITES + HIGHR_SITES
batches = [331]*len(CPN_SITES) + [322 if s.startswith("BOL")==False else 294 for s in HIGHR_SITES]
columns = ["bp_sv", "sp_sv", "bp_dim", "sp_dim", "bp_ls", "sp_ls", "nCells", "nStim", "bp_stim_dim", "sp_stim_dim", "stim_dim", "site", "batch"]
metrics = pd.DataFrame(columns=columns, index=range(len(sites)))
for i, (s, b) in enumerate(zip(sites, batches)):
    r = load_pop_metrics(site=s, batch=b, modelname=modelname)
    metrics.iloc[i] = [r["bp_sv"], r["sp_sv"], 
                    r["bp_dim95"], r["sp_dim95"],
                    r["bp_loading_sim"], r["sp_loading_sim"], 
                    r["nCells"], r["nStim"],
                    r["pca_results"]["bp_ta_dim"], r["pca_results"]["sp_ta_dim"],
                    r["pca_results"]["ta_dim"],
                    s, b]

# ===================== POPULATION STIM. DIMENSIONALITY =========================
# depends on interaction between n cells and n stimuli
f, ax = plt.subplots(1, 2, figsize=(9, 4))
offset = np.random.normal(0, 2, metrics.shape[0])
im = ax[0].scatter(metrics["nCells"], metrics["nStim"]+offset, c=metrics["sp_stim_dim"].astype(float), s=100, cmap="viridis", vmin=0, vmax=13, edgecolor="white")
f.colorbar(im, ax=ax[0])
for nstim in metrics["nStim"].unique():
    ax[0].axhline(nstim, color="k", linestyle="--", zorder=-1)
ax[0].set_title("Number of Stimulus \n Dimensions (small pupil)")
ax[0].set_xlabel("Number of Cells")
ax[0].set_ylabel("Number of Stimuli")

im = ax[1].scatter(metrics["nCells"], metrics["nStim"]+offset, c=metrics["bp_stim_dim"].astype(float), s=100, cmap="viridis", vmin=0, vmax=13, edgecolor="white")
f.colorbar(im, ax=ax[1])
for nstim in metrics["nStim"].unique():
    ax[1].axhline(nstim, color="k", linestyle="--", zorder=-1)
ax[1].set_title("Number of Stimulus \n Dimensions (large pupil)")
ax[1].set_xlabel("Number of Cells")
ax[1].set_ylabel("Number of Stimuli")

f.tight_layout()



# ===================== POPULATION FA DIMENSIONALITY =========================
f, ax = plt.subplots(1, 2, figsize=(9, 4))
offset = np.random.normal(0, 2, metrics.shape[0])
im = ax[0].scatter(metrics["nCells"], metrics["nStim"]+offset, c=metrics["sp_dim"].astype(float), s=100, cmap="viridis", vmin=0, vmax=13, edgecolor="white")
f.colorbar(im, ax=ax[0])
for nstim in metrics["nStim"].unique():
    ax[0].axhline(nstim, color="k", linestyle="--", zorder=-1)
ax[0].set_title("Number of Stimulus-Independent \n Dimensions (small pupil)")
ax[0].set_xlabel("Number of Cells")
ax[0].set_ylabel("Number of Stimuli")

im = ax[1].scatter(metrics["nCells"], metrics["nStim"]+offset, c=metrics["bp_dim"].astype(float), s=100, cmap="viridis", vmin=0, vmax=13, edgecolor="white")
f.colorbar(im, ax=ax[1])
for nstim in metrics["nStim"].unique():
    ax[1].axhline(nstim, color="k", linestyle="--", zorder=-1)
ax[1].set_title("Number of Stimulus-Independent \n Dimensions (large pupil)")
ax[1].set_xlabel("Number of Cells")
ax[1].set_ylabel("Number of Stimuli")

f.tight_layout()

# ===================== Stimulus space vs. FA space =========================
mask =  [True]*len(metrics) #
f, ax = plt.subplots(1, 2, figsize=(8, 4))

ax[0].scatter(metrics[mask]["sp_stim_dim"], metrics[mask]["sp_dim"], s=35, color=colors.SMALL, edgecolor="white")
ax[0].set_ylabel("Number of shared factors")
ax[0].set_xlabel("Number of stimulus dimensions")
ax[0].set_title("Small pupil")

ax[1].scatter(metrics[mask]["bp_stim_dim"], metrics[mask]["bp_dim"], s=35, color=colors.LARGE, edgecolor="white")
ax[1].set_ylabel("Number of shared factors")
ax[1].set_xlabel("Number of stimulus dimensions")
ax[1].set_title("Large pupil")

ll, ul = (0, 15) #ax[0].get_xlim()
ll = 0
ax[0].plot([ll, ul], [ll, ul], "k--")
ax[1].plot([ll, ul], [ll, ul], "k--")

f.tight_layout()

if savefig:
    f.savefig(fig_fn)
    f.savefig(fig_fn.replace(".svg", ".png"))

# ================== OVERLAP BETWEEN SUBSPACES ======================



# ================== SUPPLEMENTAL FIGURES ======================

# ================== STIM SPACE LRG VS. SMALL ======================
mask =  [True]*len(metrics) #
f, ax = plt.subplots(1, 1, figsize=(5, 5))

ax.scatter(metrics[mask]["bp_stim_dim"], metrics[mask]["sp_stim_dim"], edgecolor="white", color="k", s=75)
ax.plot([0, 15], [0, 15], "k--")
ax.set_xlabel("Large pupil")
ax.set_ylabel("Small pupil")
ax.set_title("Stimulus space dimensionality")


# ================== FA POPULATION METRICS ======================
mask = [True]*len(metrics) #
f, ax = plt.subplots(2, 3, figsize=(9, 6))

ax[0, 0].scatter(metrics[mask]["bp_dim"], metrics[mask]["sp_dim"], color="k", edgecolor="white")
ax[0, 0].plot([0, 15], [0, 15], "k--")
ax[0, 0].set_xlabel("Large pupil")
ax[0, 0].set_ylabel("Small pupil")
ax[0, 0].set_title("Dimenionsality")

ax[1, 0].errorbar([1], metrics[mask]["bp_dim"].mean(), yerr=metrics[mask]["bp_dim"].sem(),
                        marker="o", capsize=3, color=colors.LARGE)
ax[1, 0].errorbar([0], metrics[mask]["sp_dim"].mean(), yerr=metrics[mask]["sp_dim"].sem(),
                        marker="o", capsize=3, color=colors.SMALL)
p = np.round(ss.wilcoxon(metrics[mask]["bp_dim"], metrics[mask]["sp_dim"]).pvalue, 3)
uu = ax[1, 0].get_ylim()[-1]
ax[1, 0].plot([0, 1], [uu, uu], "k-")
ax[1, 0].text(0.3, uu+(uu/100), f"p={p}")
ax[1, 0].set_xlim((-0.5, 1.5))

ax[0, 1].scatter(metrics[mask]["bp_sv"], metrics[mask]["sp_sv"], color="k", edgecolor="white")
ax[0, 1].plot([0, 1], [0, 1], "k--")
ax[0, 1].set_xlabel("Large pupil")
ax[0, 1].set_ylabel("Small pupil")
ax[0, 1].set_title("% Shared Variance")

ax[1, 1].errorbar([1], metrics[mask]["bp_sv"].mean(), yerr=metrics[mask]["bp_sv"].sem(),
                        marker="o", capsize=3, color=colors.LARGE)
ax[1, 1].errorbar([0], metrics[mask]["sp_sv"].mean(), yerr=metrics[mask]["sp_sv"].sem(),
                        marker="o", capsize=3, color=colors.SMALL)
p = np.round(ss.wilcoxon(metrics[mask]["bp_sv"], metrics[mask]["sp_sv"]).pvalue, 3)
uu = ax[1, 1].get_ylim()[-1]
ax[1, 1].plot([0, 1], [uu, uu], "k-")
ax[1, 1].text(0.3, uu+(uu/100), f"p={p}")
ax[1, 1].set_xlim((-0.5, 1.5))

ax[0, 2].scatter(metrics[mask]["bp_ls"], metrics[mask]["sp_ls"], color="k", edgecolor="white")
ax[0, 2].plot([0, 1], [0, 1], "k--")
ax[0, 2].set_xlabel("Large pupil")
ax[0, 2].set_ylabel("Small pupil")
ax[0, 2].set_title("Loading similarity")

ax[1, 2].errorbar([1], metrics[mask]["bp_ls"].mean(), yerr=metrics[mask]["bp_ls"].sem(),
                        marker="o", capsize=3, color=colors.LARGE)
ax[1, 2].errorbar([0], metrics[mask]["sp_ls"].mean(), yerr=metrics[mask]["sp_ls"].sem(),
                        marker="o", capsize=3, color=colors.SMALL)
p = np.round(ss.wilcoxon(metrics[mask]["bp_ls"], metrics[mask]["sp_ls"]).pvalue, 3)
uu = ax[1, 2].get_ylim()[-1]
ax[1, 2].plot([0, 1], [uu, uu], "k-")
ax[1, 2].text(0.3, uu+(uu/100), f"p={p}")
ax[1, 2].set_xlim((-0.5, 1.5))

f.tight_layout()