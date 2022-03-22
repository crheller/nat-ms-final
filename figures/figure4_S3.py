"""
Third supplement to figure 4.

Look at delta dprime as a function of the response magnitude (z-scored) for each stimulus.

Idea is that this is a negative result showing that the "goodness" of the stimulus (its ability to drive neurons)
doesn't predict the diversity in decoding changes.
"""
import sys
sys.path.append("/auto/users/hellerc/code/projects/nat-ms-final")
from path_settings import DPRIME_DIR, PY_FIGURES_DIR
from global_settings import HIGHR_SITES, CPN_SITES

import charlieTools.nat_sounds_ms.decoding as decoding
from regression_helper import fit_OLS_model

import statsmodels.api as sm
import pickle
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['font.size'] = 8

modelname = 'dprime_jk10_zscore_fixtdr2-fa'
nComponents = 2
modelname331 = 'dprime_jk10_zscore_fixtdr2-fa_noiseDim-6'
nComponents331 = 8
sig_pairs_only = False
recache = True
savefig = False
fig_fn = PY_FIGURES_DIR + 'S3_figure4.svg'

np.random.seed(123)

# ============================= LOAD DPRIME =========================================
path = DPRIME_DIR
loader = decoding.DecodingResults()
df_all = []
sites = CPN_SITES + HIGHR_SITES
batches = [331]*len(CPN_SITES) + [322] * len(HIGHR_SITES)
batches = [b if s.startswith("BOL")==False else 294 for s, b in zip(sites, batches)]
for batch, site in zip(batches, sites):
    if batch == 331:
        mn = modelname331
        n_components = nComponents331
    else:
        mn = modelname
        n_components = nComponents
    try:
        fn = os.path.join(path, str(batch), site, mn+'_TDR.pickle')
        results = loader.load_results(fn, cache_path=None, recache=recache)
        _df = results.numeric_results
    except:
        print(f"WARNING!! NOT LOADING SITE {site}")

    stim = results.evoked_stimulus_pairs
    _df = _df.loc[pd.IndexSlice[stim, n_components], :]
    _df['cos_dU'] = results.slice_array_results('cos_dU_evec_test', stim, n_components, idx=(0,0))[0] #.apply(lambda x: np.arccos(x)*180/np.pi)
    _df['site'] = site

    # zscore the resp mag within site to prevent bias from number of cells
    alldata = pd.concat([_df['r1mag_test'], _df['r2mag_test']])
    m = alldata.mean()
    sd = alldata.std()
    _df['r1mag_test'] = (_df['r1mag_test'] - m) / sd
    _df['r2mag_test'] = (_df['r2mag_test'] - m) / sd

    df_all.append(_df)

df_all = pd.concat(df_all)

df_all['delta'] = (df_all['bp_dp'] - df_all['sp_dp']) / (df_all['bp_dp'] + df_all['sp_dp'])
df_all['delta_dU'] = df_all['bp_dU_mag'] - df_all['sp_dU_mag']
df_all = df_all[df_all["delta"].isna()==False]

c1 = []
c2 = []
c12 = []
c1_ci = []
c2_ci = []
c12_ci = []
r2 = []
r2_ci = []
for s in sites:
    print(f"Running regression for site {s}")
    df_dup = df_all.copy()
    df_dup = df_dup[df_dup.site==s]
    df_dup = pd.concat([df_all, df_all])
    df_dup['r1mag_test'] = pd.concat([df_all['r1mag_test'], df_all['r2mag_test']])
    df_dup['r2mag_test'] = pd.concat([df_all['r2mag_test'], df_all['r1mag_test']])

    X = df_dup[['r1mag_test', 'r2mag_test']]
    #X = pd.DataFrame()
    X['interaction'] = df_dup['r1mag_test'] * df_dup['r2mag_test']
    X -= X.mean()
    X /= X.std()
    X = sm.add_constant(X)

    y = df_dup['delta']

    res = fit_OLS_model(X, y, replace=True, nboot=10, njacks=2)
    r2.append(res["r2"]["full"])
    r2_ci.append(res["ci"]["full"])
    c1.append(res['coef']['r1mag_test'])
    c1_ci.append(res['ci_coef']['r1mag_test'])
    c2.append(res['coef']['r2mag_test'])
    c2_ci.append(res['ci_coef']['r2mag_test'])
    c12.append(res['coef']['interaction'])
    c12_ci.append(res['ci_coef']['interaction'])

bins = 40
mm = 0.75
f, ax = plt.subplots(1, 2, figsize=(10, 4))

df_all.plot.hexbin(x='r1mag_test',
                   y='r2mag_test',
                   C='delta',
                   gridsize=bins,
                   vmin=-mm, vmax=mm, cmap='bwr', ax=ax[0])
ax[0].set_xlabel(r"$|\mathbf{r}_a|$ (normalized)")
ax[0].set_ylabel(r"$|\mathbf{r}_b|$ (normalized)")
ax[0].set_title(r"$\Delta d'^2$")

# crop for visualization (z-scored, so including 2 standard deviations)
ax[0].set_xlim((-2, 2))
ax[0].set_ylim((-2, 2))

for i in range(len(r2)):
    ax[1].plot(i, r2[i], "o", color="k")
    ax[1].plot([i, i], [r2_ci[i][0], r2_ci[i][1]], color="k")

ax[1].axhline(0, linestyle="--", color="grey")
ax[1].set_ylim((-0.001, 0.01))
ax[1].set_ylabel(r"$R^2$")
ax[1].set_xlabel("Site")

f.tight_layout()

if savefig:
    f.savefig(fig_fn)
    f.savefig(fig_fn.replace(".svg", ".png"))

plt.show()