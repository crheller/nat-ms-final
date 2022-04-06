"""
Decoding results for factor analysis: Show that change in low-D shared factor space is enough to account for d-prime diversity
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
fa_s1_df = []
fa_s2_df = []
fa_s3_df = []
fa_s4_df = []
for i, (s, b) in enumerate(zip(sites, batches)):
    if b == 331:
        _noise = noise331
    else:
        _noise = noise
    rmodel = f"dprime_jk10_zscore_fixtdr2-fa{_noise}"                         # reaw data
    famodel_s1 = f"dprime_faModel.sim1_jk10_zscore_fixtdr2-fa{_noise}"      # fixed cov matrix between lrg / small
    famodel_s2 = f"dprime_faModel.sim2_jk10_zscore_fixtdr2-fa{_noise}"       # only change ind. variance. (fix abs covrirance)
    famodel_s3 = f"dprime_faModel.sim3_jk10_zscore_fixtdr2-fa{_noise}"       # only change ind. variance. (fix relative covariance (correlations))
    famodel_s4 = f"dprime_faModel.sim4_jk10_zscore_fixtdr2-fa{_noise}"  

    fn = os.path.join(DPRIME_DIR, str(b), s, rmodel+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df = results.numeric_results.loc[results.evoked_stimulus_pairs]

    fn = os.path.join(DPRIME_DIR, str(b), s, famodel_s1+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df_s1_fa = results.numeric_results.loc[results.evoked_stimulus_pairs]

    fn = os.path.join(DPRIME_DIR, str(b), s, famodel_s2+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df_s2_fa = results.numeric_results.loc[results.evoked_stimulus_pairs]

    fn = os.path.join(DPRIME_DIR, str(b), s, famodel_s3+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df_s3_fa = results.numeric_results.loc[results.evoked_stimulus_pairs]

    fn = os.path.join(DPRIME_DIR, str(b), s, famodel_s4+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df_s4_fa = results.numeric_results.loc[results.evoked_stimulus_pairs]

    df["delta"] = (df["bp_dp"]-df["sp_dp"]) / (df["bp_dp"]+df["sp_dp"])
    df["site"] = s
    df["batch"] = b

    df_s1_fa["delta"] = (df_s1_fa["bp_dp"]-df_s1_fa["sp_dp"])/(df_s1_fa["bp_dp"]+df_s1_fa["sp_dp"])
    df_s1_fa["site"] = s
    df_s1_fa["batch"] = b

    df_s2_fa["delta"] = (df_s2_fa["bp_dp"]-df_s2_fa["sp_dp"])/(df_s2_fa["bp_dp"]+df_s2_fa["sp_dp"])
    df_s2_fa["site"] = s
    df_s2_fa["batch"] = b

    df_s3_fa["delta"] = (df_s3_fa["bp_dp"]-df_s3_fa["sp_dp"])/(df_s3_fa["bp_dp"]+df_s3_fa["sp_dp"])
    df_s3_fa["site"] = s
    df_s3_fa["batch"] = b

    df_s4_fa["delta"] = (df_s4_fa["bp_dp"]-df_s4_fa["sp_dp"])/(df_s4_fa["bp_dp"]+df_s4_fa["sp_dp"])
    df_s4_fa["site"] = s
    df_s4_fa["batch"] = b

    raw_df.append(df)
    fa_s1_df.append(df_s1_fa)
    fa_s2_df.append(df_s2_fa)
    fa_s3_df.append(df_s3_fa)
    fa_s4_df.append(df_s4_fa)

df = pd.concat(raw_df)
df_s1_fa = pd.concat(fa_s1_df)
df_s2_fa = pd.concat(fa_s2_df)
df_s3_fa = pd.concat(fa_s3_df)
df_s4_fa = pd.concat(fa_s4_df)
na = df["delta"].isna()
df = df[na==False]
df_s1_fa = df_s1_fa[na==False]
df_s2_fa = df_s2_fa[na==False]
df_s3_fa = df_s3_fa[na==False]
df_s4_fa = df_s4_fa[na==False]

# compute r2 / err for each model / site / batch
df_s1_fa["err"] = abs(df["delta"] - df_s1_fa["delta"])
df_s2_fa["err"] = abs(df["delta"] - df_s2_fa["delta"])
df_s3_fa["err"] = abs(df["delta"] - df_s3_fa["delta"])
df_s4_fa["err"] = abs(df["delta"] - df_s4_fa["delta"])


for i, (s, b) in enumerate(zip(sites, batches)):
    m = (df.site==s) & (df.batch==b)
    # fit models
    X = df[m][["delta"]]
    sm.add_constant(X)
    y = df_s1_fa[m]["delta"]
    model = sm.OLS(y, X).fit()
    r2_1 = model.rsquared
    df_s1_fa.loc[(df_s1_fa.site==s) & (df_s1_fa.batch==b), "r2"] = r2_1
    y = df_s2_fa[m]["delta"]
    model = sm.OLS(y, X).fit()
    r2_2 = model.rsquared
    df_s2_fa.loc[(df_s1_fa.site==s) & (df_s1_fa.batch==b), "r2"] = r2_2
    y = df_s3_fa[m]["delta"]
    model = sm.OLS(y, X).fit()
    r2_3 = model.rsquared
    df_s3_fa.loc[(df_s1_fa.site==s) & (df_s1_fa.batch==b), "r2"] = r2_3
    y = df_s4_fa[m]["delta"]
    model = sm.OLS(y, X).fit()
    r2_4 = model.rsquared
    df_s4_fa.loc[(df_s1_fa.site==s) & (df_s1_fa.batch==b), "r2"] = r2_4

# ================================== BUILD FIGURE =======================================
f = plt.figure(figsize=(9, 8))
ms = 50
s1_ax = plt.subplot2grid((2, 5), (0, 0), colspan=2)
s2_ax = plt.subplot2grid((2, 5), (0, 2), colspan=2)
s3_ax = plt.subplot2grid((2, 5), (1, 0), colspan=2)
s4_ax = plt.subplot2grid((2, 5), (1, 2), colspan=2)
err_ax = plt.subplot2grid((2, 5), (0, 4), colspan=1)
r2_ax = plt.subplot2grid((2, 5), (1, 4), colspan=1)

# example site
mask = (df.site==site) & (df.batch==batch)

# null model
s1_ax.scatter(df[mask]["delta"], df_s1_fa[mask]["delta"], s=ms, color="tab:blue", edgecolor="white")
s1_ax.set_title(r"$\Sigma_{small}=\Sigma_{large}=\Psi_{large}$")

# ind model
s2_ax.scatter(df[mask]["delta"], df_s2_fa[mask]["delta"], s=ms, color="tab:orange", edgecolor="white")
s2_ax.set_title(r"$\Sigma_{small}=\Psi_{small}$, $\Sigma_{large}=\Psi_{large}$")

s3_ax.scatter(df[mask]["delta"], df_s3_fa[mask]["delta"], s=ms, color="tab:orange", edgecolor="white")
s3_ax.set_title(r"$\Sigma_{small}=\Psi_{small}$, $\Sigma_{large}=\Psi_{large}$")

# full model
s4_ax.scatter(df[mask]["delta"], df_s4_fa[mask]["delta"], s=ms, color="tab:green", edgecolor="white")
s4_ax.set_title(r"$\Sigma_{small}=\Sigma_{shared, small}+\Psi_{small}$,"+"\n"+"$\Sigma_{large}=\Sigma_{shared, large}+\Psi_{large}$")

for a in [s1_ax, s2_ax, s3_ax, s4_ax]:
    a.plot([-0.6, 0.6], [-0.6, 0.6], "grey", linestyle="--", zorder=-1)
    a.axvline(0, linestyle="--", color="grey", zorder=-1)
    a.axhline(0, linestyle="--", color="grey", zorder=-1)
    a.set_xlabel(r"Raw $\Delta d'^2$")
    a.set_ylabel(r"Predicted $\Delta d'^2$")

# summary across sites

# r2 
r2_ax.errorbar(0, df_s1_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_s1_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:blue", marker="o")
r2_ax.errorbar(1, df_s2_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_s2_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:orange", marker="o")
r2_ax.errorbar(2, df_s3_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_s3_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:orange", marker="o")
r2_ax.errorbar(3, df_s4_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_s4_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:green", marker="o")
r2_ax.set_xlim((-0.5, 3.5))
r2_ax.set_ylabel(r"$R^2$")

# err
err_ax.errorbar(0, df_s1_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_s1_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:blue", marker="o")
err_ax.errorbar(1, df_s2_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_s2_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:orange", marker="o")
err_ax.errorbar(2, df_s3_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_s3_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:orange", marker="o")
err_ax.errorbar(3, df_s4_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_s4_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:green", marker="o")
err_ax.set_xlim((-0.5, 3.5))
err_ax.set_ylabel(r"$|\Delta d'^2 - \Delta d'^2_{pred}|$")

f.tight_layout()

if savefig:
    f.savefig(fig_fn)
    f.savefig(fig_fn.replace(".svg", ".png"))