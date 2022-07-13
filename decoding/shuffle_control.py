"""
Compare delta dprime for shuffled pupil with real pupil
"""
import sys
sys.path.append("/auto/users/hellerc/code/projects/nat-ms-final")
from global_settings import CPN_SITES
from path_settings import DPRIME_DIR
import charlieTools.plotting as cplt
import charlieTools.nat_sounds_ms.decoding as decoding
import numpy as np
import matplotlib.pyplot as plt
import nems.db as nd
import pandas as pd
import os
import matplotlib as mpl
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['font.size'] = 12

decoder = "dprime_jk10_zscore_fixtdr2-fa"
decoder0 = "dprime_shufpup_jk10_zscore_fixtdr2-fa"
decoder00 = "dprime_shuffle_jk10_zscore_fixtdr2-fa"
batch = 331
sites, _ = nd.get_batch_sites(batch)

results = []
results0 = []
results00 = []
for site in sites:
    loader = decoding.DecodingResults()
    # raw data
    fn = os.path.join(DPRIME_DIR, str(batch), site, decoder+'_TDR.pickle')
    res = loader.load_results(fn, cache_path=None)
    df = res.numeric_results.loc[pd.IndexSlice[res.evoked_stimulus_pairs, 2], :]
    df['delta_dprime'] = (df['bp_dp'] - df['sp_dp']) / (df['bp_dp'] + df['sp_dp'])
    df['raw_delta'] = (df['bp_dp'] - df['sp_dp'])
    df['batch'] = batch
    df['site'] = site
    results.append(df)
    # shuffle pupil
    fn = os.path.join(DPRIME_DIR, str(batch), site, decoder0+'_TDR.pickle')
    res = loader.load_results(fn, cache_path=None)
    df = res.numeric_results.loc[pd.IndexSlice[res.evoked_stimulus_pairs, 2], :]
    df['delta_dprime'] = (df['bp_dp'] - df['sp_dp']) / (df['bp_dp'] + df['sp_dp'])
    df['raw_delta'] = (df['bp_dp'] - df['sp_dp'])
    df['batch'] = batch
    df['site'] = site
    results0.append(df)
    # shuffle correlations
    fn = os.path.join(DPRIME_DIR, str(batch), site, decoder00+'_TDR.pickle')
    res = loader.load_results(fn, cache_path=None)
    df = res.numeric_results.loc[pd.IndexSlice[res.evoked_stimulus_pairs, 2], :]
    df['delta_dprime'] = (df['bp_dp'] - df['sp_dp']) / (df['bp_dp'] + df['sp_dp'])
    df['raw_delta'] = (df['bp_dp'] - df['sp_dp'])
    df['batch'] = batch
    df['site'] = site
    results00.append(df)

results = pd.concat(results)
results0 = pd.concat(results0)
results00 = pd.concat(results00)

# histogram of delta dprime
bins = np.arange(-1, 1, 0.02)
f, ax = plt.subplots(1, 1, figsize=(5, 5))

ax.hist(results["delta_dprime"], bins=bins,
            histtype="step", lw=2, label="raw")
ax.hist(results0["delta_dprime"], bins=bins,
            histtype="step", lw=2, label="pup shuff")
ax.hist(results00["delta_dprime"], bins=bins,
            histtype="step", lw=2, label="trial shuff")
ax.set_xlabel("Delta d-prime", fontsize=12)
ax.legend(frameon=False, fontsize=12)

f.tight_layout()