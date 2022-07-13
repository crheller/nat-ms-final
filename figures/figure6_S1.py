"""
Decoding results for factor analysis: 
    Show results across all 6 possible simulations. Probably don't show all of these in the main 
    text.
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
metric = "raw_delta" # "delta", "raw_delta"
noise = ""
noise331 = "_noiseDim-6"
loader = decoding.DecodingResults()

raw_df = []
fa_s1_df = []
fa_s2_df = []
fa_s3_df = []
fa_s4_df = []
fa_s5_df = []
fa_s6_df = []
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
    famodel_s5 = f"dprime_faModel.ind_jk10_zscore_fixtdr2-fa{_noise}"    # old model, diag cov.
    famodel_s6 = f"dprime_faModel.ind-null_jk10_zscore_fixtdr2-fa{_noise}"  # old model, diag cov.

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

    fn = os.path.join(DPRIME_DIR, str(b), s, famodel_s5+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df_s5_fa = results.numeric_results.loc[results.evoked_stimulus_pairs]

    fn = os.path.join(DPRIME_DIR, str(b), s, famodel_s6+'_TDR.pickle')
    results = loader.load_results(fn, cache_path=None, recache=False)
    df_s6_fa = results.numeric_results.loc[results.evoked_stimulus_pairs]

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

    df_s5_fa["delta"] = (df_s5_fa["bp_dp"]-df_s5_fa["sp_dp"])/(df_s5_fa["bp_dp"]+df_s5_fa["sp_dp"])
    df_s5_fa["site"] = s
    df_s5_fa["batch"] = b
    
    df_s6_fa["delta"] = (df_s6_fa["bp_dp"]-df_s6_fa["sp_dp"])/(df_s6_fa["bp_dp"]+df_s6_fa["sp_dp"])
    df_s6_fa["site"] = s
    df_s6_fa["batch"] = b

    raw_df.append(df)
    fa_s1_df.append(df_s1_fa)
    fa_s2_df.append(df_s2_fa)
    fa_s3_df.append(df_s3_fa)
    fa_s4_df.append(df_s4_fa)
    fa_s5_df.append(df_s5_fa)
    fa_s6_df.append(df_s6_fa)

df = pd.concat(raw_df)
df_s1_fa = pd.concat(fa_s1_df)
df_s2_fa = pd.concat(fa_s2_df)
df_s3_fa = pd.concat(fa_s3_df)
df_s4_fa = pd.concat(fa_s4_df)
df_s5_fa = pd.concat(fa_s5_df)
df_s6_fa = pd.concat(fa_s6_df)
for d in [df, df_s1_fa, df_s2_fa, df_s3_fa, df_s4_fa, df_s5_fa, df_s6_fa]:
    d["raw_delta"] = d["delta"] * (d["bp_dp"] + d["sp_dp"])

# remove nans
na = (df[metric].isna() | \
        df_s1_fa[metric].isna() | \
        df_s2_fa[metric].isna() | \
        df_s3_fa[metric].isna() | \
        df_s4_fa[metric].isna() | \
        df_s5_fa[metric].isna() | \
        df_s6_fa[metric].isna())
df = df[na==False]
df_s1_fa = df_s1_fa[na==False]
df_s2_fa = df_s2_fa[na==False]
df_s3_fa = df_s3_fa[na==False]
df_s4_fa = df_s4_fa[na==False]
df_s5_fa = df_s5_fa[na==False]
df_s6_fa = df_s6_fa[na==False]

# compute r2 / err for each model / site / batch
df_s1_fa["err"] = abs(df[metric] - df_s1_fa[metric])
df_s2_fa["err"] = abs(df[metric] - df_s2_fa[metric])
df_s3_fa["err"] = abs(df[metric] - df_s3_fa[metric])
df_s4_fa["err"] = abs(df[metric] - df_s4_fa[metric])
df_s5_fa["err"] = abs(df[metric] - df_s5_fa[metric])
df_s6_fa["err"] = abs(df[metric] - df_s6_fa[metric])


for i, (s, b) in enumerate(zip(sites, batches)):
    m = (df.site==s) & (df.batch==b)
    # fit models
    X = df[m][[metric]]
    sm.add_constant(X)
    y = df_s1_fa[m][metric]
    model = sm.OLS(y, X).fit()
    r2_1 = model.rsquared
    df_s1_fa.loc[(df_s1_fa.site==s) & (df_s1_fa.batch==b), "r2"] = r2_1
    y = df_s2_fa[m][metric]
    model = sm.OLS(y, X).fit()
    r2_2 = model.rsquared
    df_s2_fa.loc[(df_s1_fa.site==s) & (df_s1_fa.batch==b), "r2"] = r2_2
    y = df_s3_fa[m][metric]
    model = sm.OLS(y, X).fit()
    r2_3 = model.rsquared
    df_s3_fa.loc[(df_s1_fa.site==s) & (df_s1_fa.batch==b), "r2"] = r2_3
    y = df_s4_fa[m][metric]
    model = sm.OLS(y, X).fit()
    r2_4 = model.rsquared
    df_s4_fa.loc[(df_s1_fa.site==s) & (df_s1_fa.batch==b), "r2"] = r2_4
    y = df_s5_fa[m][metric]
    model = sm.OLS(y, X).fit()
    r2_5 = model.rsquared
    df_s5_fa.loc[(df_s1_fa.site==s) & (df_s1_fa.batch==b), "r2"] = r2_5
    y = df_s6_fa[m][metric]
    model = sm.OLS(y, X).fit()
    r2_6 = model.rsquared
    df_s6_fa.loc[(df_s1_fa.site==s) & (df_s1_fa.batch==b), "r2"] = r2_6

# ================================== BUILD FIGURE =======================================
# for this, don't bother with example site. Just summary, I think.
f = plt.figure(figsize=(8, 4))
err_ax = plt.subplot2grid((1, 2), (0, 0))
r2_ax = plt.subplot2grid((1, 2), (0, 1))


# summary across sites

# r2 
r2_ax.errorbar(0, df_s6_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_s6_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:blue", marker="o")
r2_ax.errorbar(1, df_s1_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_s1_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:orange", marker="o")
r2_ax.errorbar(2, df_s5_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_s5_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:green", marker="o")
r2_ax.errorbar(3, df_s3_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_s3_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:purple", marker="o")
r2_ax.errorbar(4, df_s2_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_s2_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:red", marker="o")
r2_ax.errorbar(5, df_s4_fa.groupby(by=["site", "batch"]).mean()["r2"].mean(), 
                    df_s4_fa.groupby(by=["site", "batch"]).mean()["r2"].sem(), 
                            capsize=2, color="tab:grey", marker="o")
r2_ax.set_xlim((-0.5, 5.5))
r2_ax.set_ylabel(r"$R^2$")

# err
err_ax.errorbar(0, df_s6_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_s6_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:blue", marker="o")
err_ax.errorbar(1, df_s1_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_s1_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:orange", marker="o")
err_ax.errorbar(2, df_s5_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_s5_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:green", marker="o")
err_ax.errorbar(3, df_s3_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_s3_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:purple", marker="o")
err_ax.errorbar(4, df_s2_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_s2_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:red", marker="o")
err_ax.errorbar(5, df_s4_fa.groupby(by=["site", "batch"]).mean()["err"].mean(), 
                    df_s4_fa.groupby(by=["site", "batch"]).mean()["err"].sem(), 
                            capsize=2, color="tab:grey", marker="o")
err_ax.set_xlim((-0.5, 5.5))
err_ax.set_ylabel(r"$|\Delta d'^2 - \Delta d'^2_{pred}|$")

f.tight_layout()


if savefig:
    f.savefig(fig_fn)
    f.savefig(fig_fn.replace(".svg", ".png"))