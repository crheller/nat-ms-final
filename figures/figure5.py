"""
Factor analysis model fitting
    - activity is low dimensional
Supplementals
    - population level metric changes are variable, but all tend to correspond
        to decreased noise correlation
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
sys.path.append("/auto/users/hellerc/code/projects/nat_pupil_ms/")
from path_settings import PY_FIGURES_DIR
from global_settings import CPN_SITES, HIGHR_SITES
import colors
import sys
sys.path.append("/auto/users/hellerc/code/projects/nat-ms-final/FactorAnalysis")
from FactorAnalysis.loader import load_pop_metrics

savefig = True
fig_fn = PY_FIGURES_DIR + "figure5.svg"
fig_fnS1 = PY_FIGURES_DIR + "S1_figure5.svg"

modelname = "factor_analysis"

sites = CPN_SITES + HIGHR_SITES
batches = [331]*len(CPN_SITES) + [322 if s.startswith("BOL")==False else 294 for s in HIGHR_SITES]
columns = ["bp_sv", "sp_sv", "bp_dim", "sp_dim", "bp_dim_sem", "sp_dim_sem", "bp_ls", "sp_ls", "all_dim", "all_dim_sem", "site", "batch", "nCells"]
metrics = pd.DataFrame(columns=columns, index=range(len(sites)))
metrics_all = pd.DataFrame(columns=columns, index=range(len(sites)))
for i, (s, b) in enumerate(zip(sites, batches)):
    r = load_pop_metrics(site=s, batch=b, modelname=modelname)
    metrics.iloc[i] = [r["bp_sv"], r["sp_sv"], r["bp_dim95"], r["sp_dim95"], r["bp_dim_sem"], r["sp_dim_sem"], r["bp_loading_sim"], r["sp_loading_sim"], r["all_dim95"], r["all_dim95_sd"], s, b, r["nCells"]]
    metrics_all.iloc[i] = [r["final_fit"]["bp_sv_all"], r["final_fit"]["sp_sv_all"], r["final_fit"]["bp_dim95_all"], r["final_fit"]["sp_dim95_all"], r["bp_dim_sem"], r["sp_dim_sem"], r["final_fit"]["bp_ls_all"], r["final_fit"]["sp_ls_all"], r["final_fit"]["all_dim95_all"], r["all_dim_sem"], s, b, r["nCells"]]


# ===================== POPULATION DIMENSIONALITY =========================
mask =  [True]*len(metrics) #
f, ax = plt.subplots(1, 3, figsize=(8, 2.66))

ax[0].scatter(metrics[mask]["nCells"], metrics[mask]["bp_dim"], s=35, color=colors.LARGE, edgecolor="white")
ax[0].set_ylabel("Number of shared factors")
ax[0].set_xlabel("Number of recorded units")
ax[0].set_title("Large pupil")

ax[1].scatter(metrics[mask]["nCells"], metrics[mask]["sp_dim"], s=35, color=colors.SMALL, edgecolor="white")
ax[1].set_ylabel("Number of shared factors")
ax[1].set_xlabel("Number of recorded units")
ax[1].set_title("Small pupil")

ll, ul = ax[0].get_xlim()
ll = 0
ax[0].plot([ll, ul], [ll, ul], "k--")
ax[1].plot([ll, ul], [ll, ul], "k--")

ax[2].errorbar([0], 
            [metrics[mask]["bp_dim"].mean()],
            yerr=[metrics[mask]["bp_dim"].sem()],
            capsize=2, lw=2, marker="o", color=colors.LARGE)
ax[2].errorbar([1], 
            [metrics[mask]["sp_dim"].mean()],
            yerr=[metrics[mask]["sp_dim"].sem()],
            capsize=3, lw=2, marker="o", color=colors.SMALL)
ax[2].set_xticks([0, 1])
ax[2].set_xticklabels(["Large", "Small"])
ax[2].set_xlim((-2, 3))
ax[2].set_ylabel("Number of shared factors")
top = ax[2].get_ylim()[-1]
ax[2].plot([0, 1], [top, top], lw=2, color="k")
pval = ss.wilcoxon(metrics[mask]["bp_dim"], metrics[mask]["sp_dim"]).pvalue
ax[2].text(-0.1, top+0.1, r"$p = %s$" % (np.round(pval, 3)), fontsize=8)


f.tight_layout()

if savefig:
    f.savefig(fig_fn)
    f.savefig(fig_fn.replace(".svg", ".png"))

# ================== SUPPLEMENTAL POPULATION METRICS ======================
# WORK IN PROGRESS -- WANT TO REFINE THIS IF INCLUDING IT
mask =  [True]*len(metrics) #
f, ax = plt.subplots(1, 3, figsize=(16, 5))

ax[0].scatter(metrics[mask]["bp_sv"], metrics[mask]["sp_sv"], edgecolor="white")
ax[0].scatter(metrics[mask]["bp_sv"].mean(), metrics[mask]["sp_sv"].mean(), s=100, edgecolor="white")
ax[0].plot([0, 1], [0, 1], "k--")
ax[0].set_xlabel("Large pupil")
ax[0].set_ylabel("Small pupil")
ax[0].set_title("% sv")

ax[1].scatter(metrics[mask]["bp_ls"], metrics[mask]["sp_ls"], edgecolor="white")
ax[1].scatter(metrics[mask]["bp_ls"].mean(), metrics[mask]["sp_ls"].mean(), s=100, edgecolor="white")
ax[1].plot([0, 1], [0, 1], "k--")
ax[1].set_xlabel("Large pupil")
ax[1].set_ylabel("Small pupil")
ax[1].set_title("Loading similarity")

ax[2].scatter(metrics[mask]["bp_dim"], metrics[mask]["sp_dim"], edgecolor="white")
ax[2].scatter(metrics[mask]["bp_dim"].mean(), metrics[mask]["sp_dim"].mean(), s=100, edgecolor="white")
ax[2].plot([0, 20], [0, 20], "k--")
ax[2].set_xlabel("Large pupil")
ax[2].set_ylabel("Small pupil")
ax[2].set_title("nDim")

if savefig:
    f.savefig(fig_fnS1)
    f.savefig(fig_fnS1.replace(".svg", ".png"))