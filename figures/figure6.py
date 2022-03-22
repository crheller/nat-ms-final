"""
Decoding results for factor analysis.
    * Example site (scatter plot - more intuitive for explaining things)
    * Summary (report R2 and absolute error) for each model
"""
import sys
sys.path.append("/auto/users/hellerc/code/projects/nat-ms-final/")
from path_settings import DPRIME_DIR, PY_FIGURES_DIR
from global_settings import HIGHR_SITES, CPN_SITES
import confusion_matrix_helper as chelp

import charlieTools.nat_sounds_ms.decoding as decoding
from charlieTools.statistics import get_direct_prob, get_bootstrapped_sample

import os
import pandas as pd
import scipy.stats as ss
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['font.size'] = 8

np.random.seed(123)
savefig = False
fig_fn = PY_FIGURES_DIR + "figure6.svg"

sites = CPN_SITES + HIGHR_SITES
batches = [331]*len(CPN_SITES) + [322 if s.startswith("BOL")==False else 294 for s in HIGHR_SITES]
site = 'ARM033a' #'DRX007a.e65:128' #'DRX008b.e65:128' #'DRX007a.e65:128'
batch = 331

noise = ""
noise331 = "_noiseDim-6"
loader = decoding.DecodingResults()

raw_df = []
fa_null_df = []
fa_ind_df = []
fa_full_df = []
for i, (s, b) in enumerate(zip(sites, batches)):
    if b == 331:
        _noise = noise331
    else:
        _noise = noise
    rmodel = f"dprime_jk10_zscore_fixtdr2-fa{_noise}"                    # reaw data
    famodel_null = f"dprime_faModel.ind-null_jk10_zscore_fixtdr2-fa{_noise}" # fixed cov matrix between lrg / small
    famodel_ind = f"dprime_faModel.ind_jk10_zscore_fixtdr2-fa{_noise}"   # only change ind. variance. (diag cov matrix)
    #famodel = f"dprime_faModel.rr1_jk10_zscore_fixtdr2-fa{_noise}"           # full (reduced) rank cov matrix
    famodel = f"dprime_faModel_jk10_zscore_fixtdr2-fa{_noise}"  

    fn = os.path.join(DPRIME_DIR, str(b), s, rmodel+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df = results.numeric_results.loc[results.evoked_stimulus_pairs]

    fn = os.path.join(DPRIME_DIR, str(b), s, famodel_null+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df_null_fa = results.numeric_results.loc[results.evoked_stimulus_pairs]

    fn = os.path.join(DPRIME_DIR, str(b), s, famodel_ind+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df_ind_fa = results.numeric_results.loc[results.evoked_stimulus_pairs]

    fn = os.path.join(DPRIME_DIR, str(b), s, famodel+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df_fa = results.numeric_results.loc[results.evoked_stimulus_pairs]

    df["delta"] = (df["bp_dp"]-df["sp_dp"]) / (df["bp_dp"]+df["sp_dp"])
    df["site"] = s
    df["batch"] = b

    df_null_fa["delta"] = (df_null_fa["bp_dp"]-df_null_fa["sp_dp"])/(df_null_fa["bp_dp"]+df_null_fa["sp_dp"])
    df_null_fa["site"] = s
    df_null_fa["batch"] = b

    df_ind_fa["delta"] = (df_ind_fa["bp_dp"]-df_ind_fa["sp_dp"])/(df_ind_fa["bp_dp"]+df_ind_fa["sp_dp"])
    df_ind_fa["site"] = s
    df_ind_fa["batch"] = b

    df_fa["delta"] = (df_fa["bp_dp"]-df_fa["sp_dp"])/(df_fa["bp_dp"]+df_fa["sp_dp"])
    df_fa["site"] = s
    df_fa["batch"] = b

    raw_df.append(df)
    fa_null_df.append(df_null_fa)
    fa_ind_df.append(df_ind_fa)
    fa_full_df.append(df_fa)

df = pd.concat(raw_df)
df_null_fa = pd.concat(fa_null_df)
df_ind_fa = pd.concat(fa_ind_df)
df_full_fa = pd.concat(fa_full_df)
na = df["delta"].isna()
df = df[na==False]
df_null_fa = df_null_fa[na==False]
df_ind_fa = df_ind_fa[na==False]
df_full_fa = df_full_fa[na==False]

# compute r2 / err for each model / site / batch
df_null_fa["err"] = abs(df["delta"] - df_null_fa["delta"])
df_ind_fa["err"] = abs(df["delta"] - df_ind_fa["delta"])
df_full_fa["err"] = abs(df["delta"] - df_full_fa["delta"])


for i, (s, b) in enumerate(zip(sites, batches)):
    m = (df.site==s) & (df.batch==b)
    # fit models
    X = df[m][["delta"]]
    sm.add_constant(X)
    y = df_null_fa[m]["delta"]
    model = sm.OLS(y, X).fit()
    r2_1 = model.rsquared
    df_null_fa.loc[(df_null_fa.site==s) & (df_null_fa.batch==b), "r2"] = r2_1
    y = df_ind_fa[m]["delta"]
    model = sm.OLS(y, X).fit()
    r2_2 = model.rsquared
    df_ind_fa.loc[(df_null_fa.site==s) & (df_null_fa.batch==b), "r2"] = r2_2
    y = df_full_fa[m]["delta"]
    model = sm.OLS(y, X).fit()
    r2_3 = model.rsquared
    df_full_fa.loc[(df_null_fa.site==s) & (df_null_fa.batch==b), "r2"] = r2_3

# ================================== BUILD FIGURE =======================================
f = plt.figure(figsize=(10, 5))

null_ax = plt.subplot2grid((2, 8), (0, 0), colspan=2)
ind_ax = plt.subplot2grid((2, 8), (0, 2), colspan=2)
full_ax = plt.subplot2grid((2, 8), (0, 4), colspan=2)
err_ax = plt.subplot2grid((2, 8), (0, 6), colspan=1)
r2_ax = plt.subplot2grid((2, 8), (0, 7), colspan=1)

# example site
mask = (df.site==site) & (df.batch==batch)

# null model
null_ax.scatter(df[mask]["delta"], df_null_fa[mask]["delta"], s=15, color="tab:blue", edgecolor="white")
null_ax.set_title(r"$\Sigma_{small}=\Sigma_{large}=\Psi_{large}$")

# ind model
ind_ax.scatter(df[mask]["delta"], df_ind_fa[mask]["delta"], s=15, color="tab:orange", edgecolor="white")
ind_ax.set_title(r"$\Sigma_{small}=\Psi_{small}$, $\Sigma_{large}=\Psi_{large}$")

# full model
full_ax.scatter(df[mask]["delta"], df_full_fa[mask]["delta"], s=15, color="tab:green", edgecolor="white")
full_ax.set_title(r"$\Sigma_{small}=\Sigma_{shared, small}+\Psi_{small}$,"+"\n"+"$\Sigma_{large}=\Sigma_{shared, large}+\Psi_{large}$")

for a in [null_ax, ind_ax, full_ax]:
    a.plot([-0.6, 0.6], [-0.6, 0.6], "grey", linestyle="--", zorder=-1)
    a.axvline(0, linestyle="--", color="grey", zorder=-1)
    a.axhline(0, linestyle="--", color="grey", zorder=-1)
    a.set_xlabel(r"Raw $\Delta d'^2$")
    a.set_ylabel(r"Predicted $\Delta d'^2$")

# summary across sites

# r2 
r2_ax.errorbar(0, df_null_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_null_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:blue", marker="o")
r2_ax.errorbar(1, df_ind_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_ind_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:orange", marker="o")
r2_ax.errorbar(2, df_full_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_full_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:green", marker="o")
r2_ax.set_xlim((-0.5, 2.5))
r2_ax.set_ylabel(r"$R^2$")

# err
err_ax.errorbar(0, df_null_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_null_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:blue", marker="o")
err_ax.errorbar(1, df_ind_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_ind_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:orange", marker="o")
err_ax.errorbar(2, df_full_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_full_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:green", marker="o")
err_ax.set_xlim((-0.5, 2.5))
err_ax.set_ylabel(r"$|\Delta d'^2 - \Delta d'^2_{pred}|$")

f.tight_layout()

if savefig:
    f.savefig(fig_fn)
    f.savefig(fig_fn.replace(".svg", ".png"))