"""
Exploratory summary of difference in large vs. small pupil noise / signal 
subspaces
""" 
import nems.db as nd
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

#CPN_SITES, _ = nd.get_batch_sites(331)

modelname = "factor_analysis_pca_evoked"
sites = CPN_SITES + HIGHR_SITES
batches = [331]*len(CPN_SITES) + [322 if s.startswith("BOL")==False else 294 for s in HIGHR_SITES]
columns = ["bp_overlap", "sp_overlap", "stim_overlap", "noise_overlap",
                "noise_cos_sim", "stim_cos_sim",
                "bp_sv", "sp_sv", "bp_dim", "sp_dim", 
                "bp_ls", "sp_ls", "nCells", "nStim", 
                "bp_stim_dim", "sp_stim_dim", "stim_dim", 
                "site", "batch"]
metrics = pd.DataFrame(columns=columns, index=range(len(sites)))
for i, (s, b) in enumerate(zip(sites, batches)):
    try:
        r = load_pop_metrics(site=s, batch=b, modelname=modelname)

        # compute subspace overlap (stim vs. noise space)
        stim_big = r["pca_results"]["pca_ta_large.components_"]
        stim_small = r["pca_results"]["pca_ta_small.components_"]
        noise_big = r["final_fit"]["fa_big.components_"]
        noise_small = r["final_fit"]["fa_small.components_"]

        # calc big overlap btwn stim vs. noise
        if noise_big.shape[0] > stim_big.shape[0]:
            ee = stim_big.shape[0]
            big_overlap = subspace_overlap(stim_big[:ee, :], noise_big[:ee, :])
        elif noise_big.shape[0] < stim_big.shape[0]:
            ee = noise_big.shape[0]
            big_overlap = subspace_overlap(stim_big[:ee, :], noise_big[:ee, :])
        else:
            big_overlap = subspace_overlap(stim_big, noise_big)
        
        # calc small overlap between stim vs. noise
        if noise_small.shape[0] > stim_small.shape[0]:
            ee = stim_small.shape[0]
            small_overlap = subspace_overlap(stim_small[:ee, :], noise_small[:ee, :])
        elif noise_small.shape[0] < stim_small.shape[0]:
            ee = noise_small.shape[0]
            small_overlap = subspace_overlap(stim_small[:ee, :], noise_small[:ee, :])
        else:
            small_overlap = subspace_overlap(stim_small, noise_small)

        # calc noise overlap
        if noise_small.shape[0] > noise_big.shape[0]:
            ee = noise_big.shape[0]
            noise_overlap = subspace_overlap(noise_small[:ee, :], noise_small[:ee, :])
        elif noise_small.shape[0] < noise_big.shape[0]:
            ee = noise_small.shape[0]
            noise_overlap = subspace_overlap(noise_big[:ee, :], noise_small[:ee, :])
        else:
            noise_overlap = subspace_overlap(noise_big, noise_small)

        # calc stim overlap
        if stim_small.shape[0] > stim_big.shape[0]:
            ee = stim_big.shape[0]
            stim_overlap = subspace_overlap(stim_small[:ee, :], stim_small[:ee, :])
        elif stim_small.shape[0] < stim_big.shape[0]:
            ee = stim_small.shape[0]
            stim_overlap = subspace_overlap(stim_big[:ee, :], stim_small[:ee, :])
        else:
            stim_overlap = subspace_overlap(stim_big, stim_small)

        noise_cos_sim = abs(noise_small[0, :].dot(noise_big[0, :].T))
        stim_cos_sim = abs(stim_small[0, :].dot(stim_big[0, :].T))

        metrics.iloc[i] = [big_overlap, small_overlap, 
                        stim_overlap, noise_overlap, 
                        noise_cos_sim, stim_cos_sim,
                        r["bp_sv"], r["sp_sv"], 
                        r["bp_dim95"], r["sp_dim95"],
                        r["bp_loading_sim"], r["sp_loading_sim"], 
                        r["nCells"], r["nStim"],
                        r["pca_results"]["bp_ta_dim"], r["pca_results"]["sp_ta_dim"],
                        r["pca_results"]["ta_dim"],
                        s, b]
    except: 
        print(f"site: {s} not found")

metrics0 = pd.DataFrame(columns=columns, index=range(len(sites)))
for i, (s, b) in enumerate(zip(sites, batches)):
    try:
        r = load_pop_metrics(site=s, batch=b, modelname=modelname+"_shuff")

        # compute subspace overlap (stim vs. noise space)
        stim_big = r["pca_results"]["pca_ta_large.components_"]
        stim_small = r["pca_results"]["pca_ta_small.components_"]
        noise_big = r["final_fit"]["fa_big.components_"]
        noise_small = r["final_fit"]["fa_small.components_"]

        # calc big overlap btwn stim vs. noise
        if noise_big.shape[0] > stim_big.shape[0]:
            ee = stim_big.shape[0]
            big_overlap = subspace_overlap(stim_big[:ee, :], noise_big[:ee, :])
        elif noise_big.shape[0] < stim_big.shape[0]:
            ee = noise_big.shape[0]
            big_overlap = subspace_overlap(stim_big[:ee, :], noise_big[:ee, :])
        else:
            big_overlap = subspace_overlap(stim_big, noise_big)
        
        # calc small overlap between stim vs. noise
        if noise_small.shape[0] > stim_small.shape[0]:
            ee = stim_small.shape[0]
            small_overlap = subspace_overlap(stim_small[:ee, :], noise_small[:ee, :])
        elif noise_small.shape[0] < stim_small.shape[0]:
            ee = noise_small.shape[0]
            small_overlap = subspace_overlap(stim_small[:ee, :], noise_small[:ee, :])
        else:
            small_overlap = subspace_overlap(stim_small, noise_small)

        # calc noise overlap
        if noise_small.shape[0] > noise_big.shape[0]:
            ee = noise_big.shape[0]
            noise_overlap = subspace_overlap(noise_small[:ee, :], noise_small[:ee, :])
        elif noise_small.shape[0] < noise_big.shape[0]:
            ee = noise_small.shape[0]
            noise_overlap = subspace_overlap(noise_big[:ee, :], noise_small[:ee, :])
        else:
            noise_overlap = subspace_overlap(noise_big, noise_small)

        # calc stim overlap
        if stim_small.shape[0] > stim_big.shape[0]:
            ee = stim_big.shape[0]
            stim_overlap = subspace_overlap(stim_small[:ee, :], stim_small[:ee, :])
        elif stim_small.shape[0] < stim_big.shape[0]:
            ee = stim_small.shape[0]
            stim_overlap = subspace_overlap(stim_big[:ee, :], stim_small[:ee, :])
        else:
            stim_overlap = subspace_overlap(stim_big, stim_small)

        noise_cos_sim = abs(noise_small[0, :].dot(noise_big[0, :].T))
        stim_cos_sim = abs(stim_small[0, :].dot(stim_big[0, :].T))

        metrics0.iloc[i] = [big_overlap, small_overlap, 
                        stim_overlap, noise_overlap, 
                        noise_cos_sim, stim_cos_sim,
                        r["bp_sv"], r["sp_sv"], 
                        r["bp_dim95"], r["sp_dim95"],
                        r["bp_loading_sim"], r["sp_loading_sim"], 
                        r["nCells"], r["nStim"],
                        r["pca_results"]["bp_ta_dim"], r["pca_results"]["sp_ta_dim"],
                        r["pca_results"]["ta_dim"],
                        s, b]
    except: 
        print(f"site: {s} not found")



f, ax = plt.subplots(2, 2, figsize=(10, 10))

ax[0, 0].scatter(metrics[metrics.batch==322]["bp_overlap"],
                metrics[metrics.batch==322]["sp_overlap"])
ax[0, 0].scatter(metrics[metrics.batch==331]["bp_overlap"],
                metrics[metrics.batch==331]["sp_overlap"])
ax[0, 0].plot([0, 1], [0, 1], "k--")
ax[0, 0].set_xlabel("Big pupil")
ax[0, 0].set_ylabel("Small pupil")
ax[0, 0].set_title("Noise vs. stim overlap")

ax[0, 1].scatter(metrics[metrics.batch==322]["noise_overlap"],
                metrics[metrics.batch==322]["stim_overlap"])
ax[0, 1].scatter(metrics[metrics.batch==331]["noise_overlap"],
                metrics[metrics.batch==331]["stim_overlap"])
ax[0, 1].plot([0, 1], [0, 1], "k--")
ax[0, 1].set_xlabel("Noise overlap")
ax[0, 1].set_ylabel("Stim overlap")
ax[0, 1].set_title("Noise vs. stim reliability")

ax[1, 0].scatter(metrics0[metrics0.batch==322]["bp_overlap"],
                metrics0[metrics0.batch==322]["sp_overlap"])
ax[1, 0].scatter(metrics0[metrics0.batch==331]["bp_overlap"],
                metrics0[metrics0.batch==331]["sp_overlap"])
ax[1, 0].plot([0, 1], [0, 1], "k--")
ax[1, 0].set_xlabel("Big pupil")
ax[1, 0].set_ylabel("Small pupil")
ax[1, 0].set_title("Noise vs. stim overlap")

ax[1, 1].scatter(metrics0[metrics0.batch==322]["noise_overlap"],
                metrics0[metrics0.batch==322]["stim_overlap"])
ax[1, 1].scatter(metrics0[metrics0.batch==331]["noise_overlap"],
                metrics0[metrics0.batch==331]["stim_overlap"])
ax[1, 1].plot([0, 1], [0, 1], "k--")
ax[1, 1].set_xlabel("Noise overlap")
ax[1, 1].set_ylabel("Stim overlap")
ax[1, 1].set_title("Noise vs. stim reliability")

f.tight_layout()

# direct compare shuffled vs. not
# full space overlaps
f, ax = plt.subplots(1, 2, figsize=(10, 5))

ax[0].scatter(metrics[metrics.batch==322]["noise_overlap"],
                metrics0[metrics0.batch==322]["noise_overlap"])
ax[0].scatter(metrics[metrics.batch==331]["noise_overlap"],
                metrics0[metrics0.batch==331]["noise_overlap"])
ax[0].plot([0, 1], [0, 1], "k--")
ax[0].set_xlabel("raw")
ax[0].set_ylabel("shuffled pupil")
ax[0].set_title("noise overlap")

ax[1].scatter(metrics[metrics.batch==322]["stim_overlap"],
                metrics0[metrics0.batch==322]["stim_overlap"])
ax[1].scatter(metrics[metrics.batch==331]["stim_overlap"],
                metrics0[metrics0.batch==331]["stim_overlap"])
ax[1].plot([0.6, 1], [0.6, 1], "k--")
ax[1].set_xlabel("raw")
ax[1].set_ylabel("shuffled pupil")
ax[1].set_title("stim overlap")

f.tight_layout()

# Factor 1 cosine similarity
f, ax = plt.subplots(1, 2, figsize=(10, 5))

ax[0].scatter(metrics[metrics.batch==322]["noise_cos_sim"],
                metrics0[metrics0.batch==322]["noise_cos_sim"])
ax[0].scatter(metrics[metrics.batch==331]["noise_cos_sim"],
                metrics0[metrics0.batch==331]["noise_cos_sim"])
ax[0].plot([0, 1], [0, 1], "k--")
ax[0].set_xlabel("raw")
ax[0].set_ylabel("shuffled pupil")
ax[0].set_title("noise cos sim")

ax[1].scatter(metrics[metrics.batch==322]["stim_cos_sim"],
                metrics0[metrics0.batch==322]["stim_cos_sim"])
ax[1].scatter(metrics[metrics.batch==331]["stim_cos_sim"],
                metrics0[metrics0.batch==331]["stim_cos_sim"])
ax[1].plot([0.6, 1], [0.6, 1], "k--")
ax[1].set_xlabel("raw")
ax[1].set_ylabel("shuffled pupil")
ax[1].set_title("stim cos sim")

f.tight_layout()

# shared variance / dimensionality / loading sim
f, ax = plt.subplots(2, 3, figsize=(12, 8))

for a, k, t in zip(ax[0, :], ["_sv", "_dim", "_ls"], ["shared var.", "dimensionality", "loading sim"]):
    a.scatter(
        metrics[metrics.batch==322]["bp"+k],
        metrics[metrics.batch==322]["sp"+k]
    )
    a.scatter(
        metrics[metrics.batch==331]["bp"+k],
        metrics[metrics.batch==331]["sp"+k]
    )
    mi = np.min(a.get_xlim() + a.get_ylim())
    ma = np.max(a.get_xlim() + a.get_ylim())
    a.plot([mi, ma], [mi, ma], "k--")
    a.set_title(t)
    a.set_xlabel("big")
    a.set_ylabel("small")

for a, k, t in zip(ax[1, :], ["_sv", "_dim", "_ls"], ["shared var.", "dimensionality", "loading sim"]):
    a.scatter(
        metrics0[metrics0.batch==322]["bp"+k],
        metrics0[metrics0.batch==322]["sp"+k]
    )
    a.scatter(
        metrics0[metrics0.batch==331]["bp"+k],
        metrics0[metrics0.batch==331]["sp"+k]
    )
    mi = np.min(a.get_xlim() + a.get_ylim())
    ma = np.max(a.get_xlim() + a.get_ylim())
    a.plot([mi, ma], [mi, ma], "k--")
    a.set_title(t)
    a.set_xlabel("big")
    a.set_ylabel("small")

f.tight_layout()