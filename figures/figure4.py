"""
Example "confusion" matrix plotting dprime / population responses for all
stimuli at a site. Highlight that changes are not explained by first order 
response properties (e.g. PC_1 response variance)

Also, show big /small pupil d-pdrime scatter plot
"""
import sys
sys.path.append("/auto/users/hellerc/code/projects/nat-ms-final/")
from path_settings import DPRIME_DIR, PY_FIGURES_DIR
from global_settings import LOWR_SITES, HIGHR_SITES, CPN_SITES
import confusion_matrix_helper as chelp

import charlieTools.nat_sounds_ms.decoding as decoding
from charlieTools.statistics import get_direct_prob, get_bootstrapped_sample

import scipy.stats as ss
import pickle
import os
import pandas as pd
from scipy.stats import gaussian_kde
from itertools import combinations
from sklearn.decomposition import PCA
from scipy.io import wavfile
import nems.analysis.gammatone.gtgram as gt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['font.size'] = 8

np.random.seed(123)
savefig = True
fig_fn = PY_FIGURES_DIR + 'figure4.svg'
figS1 = PY_FIGURES_DIR + 'figure4_S1.svg'
figS2 = PY_FIGURES_DIR + 'figure4_S2.svg'

modelname = 'dprime_jk10_zscore_fixtdr2-fa'
nComponents = 2
modelname331 = 'dprime_jk10_zscore_fixtdr2-fa_noiseDim-6'
nComponents331 = 8
site = 'DRX008b.e65:128' #'DRX007a.e65:128' #'DRX008b.e65:128' #'DRX007a.e65:128'
batch = 289
duration = 1   # length of stimuli (for making spectrograms)
prestim = 0.25 # legnth of stimuli (for making spectrograms)
pc_div  = 8 # how much space to leave around conf. matrix edge for PC (bigger this is, the less space. Default = 16)
soundpath = '/auto/users/hellerc/code/baphy/Config/lbhb/SoundObjects/@NaturalSounds/sounds_set4/'

# get decoding results
loader = decoding.DecodingResults()
fn = os.path.join(DPRIME_DIR, str(batch), site, modelname+'_TDR.pickle')
results = loader.load_results(fn, cache_path=None)
df = results.numeric_results.loc[results.evoked_stimulus_pairs]

X, sp_bins, X_pup, pup_mask, epochs = decoding.load_site(site=site, batch=batch, return_epoch_list=True)
ncells = X.shape[0]
nreps = X.shape[1]
nstim = X.shape[2]
nbins = X.shape[3]
sp_bins = sp_bins.reshape(1, sp_bins.shape[1], nstim * nbins)
nstim = nstim * nbins

# ============================= LOAD DPRIME =========================================
path = DPRIME_DIR
loader = decoding.DecodingResults()
recache = False
df_all = []
sites = CPN_SITES + HIGHR_SITES
batches = [331] * len(CPN_SITES) + [322] * len(HIGHR_SITES)
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
        results = loader.load_results(fn, cache_path=None)
        _df = results.numeric_results
    except:
        print(f"WARNING!! NOT LOADING SITE {site}")

    stim = results.evoked_stimulus_pairs
    _df = _df.loc[pd.IndexSlice[stim, n_components], :]
    _df['site'] = site
    _df["batch"] = batch
    df_all.append(_df)

df_all = pd.concat(df_all)

# =========================== generate a list of stim pairs ==========================
# these are the indices of the decoding results dataframes
all_combos = list(combinations(range(nstim), 2))
spont_bins = np.argwhere(sp_bins[0, 0, :])
spont_combos = [c for c in all_combos if (c[0] in spont_bins) & (c[1] in spont_bins)]
ev_ev_combos = [c for c in all_combos if (c[0] not in spont_bins) & (c[1] not in spont_bins)]
spont_ev_combos = [c for c in all_combos if (c not in ev_ev_combos) & (c not in spont_combos)]

X = X.reshape(ncells, nreps, nstim)
pup_mask = pup_mask.reshape(1, nreps, nstim)
ev_bins = list(set(range(X.shape[-1])).difference(set(spont_bins.squeeze())))
Xev = X[:, :, ev_bins]

# ============================= DO PCA ================================
Xu = Xev.mean(axis=1)
spont = X[:, :, spont_bins.squeeze()].mean(axis=1).mean(axis=-1, keepdims=True)
Xu_center = Xu - spont # subtract spont
pca = PCA()
pca.fit(Xu_center.T)

spont = spont[:, :, np.newaxis] # for subtracting from single trial data
X_spont = X - spont
proj = (X_spont).T.dot(pca.components_.T)

# ================== BUILD SPECTROGRAM FOR FIGURE ======================
# (high res)
stimulus = []
for i, epoch in enumerate(epochs):
    print(f"Building spectrogram {i} / {len(epochs)} for epoch: {epoch}")
    soundfile = soundpath + epoch.strip('STIM_')
    fs, data = wavfile.read(soundfile)
    # pad / crop data
    data = data[:int(duration * fs)]
    spbins = int(prestim * fs)
    data = np.concatenate((np.zeros(spbins), data))
    spec = gt.gtgram(data, fs, 0.004, 0.002, 100, 0)
    stimulus.append(spec)
stimulus = np.concatenate(stimulus, axis=-1)
stim_fs = 500

# ================================== BUILD FIGURE =======================================
f = plt.figure(figsize=(13, 6))

gs = mpl.gridspec.GridSpec(2, 12, width_ratios=np.ones(12), height_ratios=[1, 0.05],
         wspace=0.0, hspace=0.0, top=0.9, bottom=0.1, left=0.0, right=1.0)
bp = f.add_subplot(gs[0, 0:3])
sp = f.add_subplot(gs[0, 3:6])
diff = f.add_subplot(gs[0, 6:9])
scax = f.add_subplot(gs[0, 9:])

# get big pupil / small pupil projected response, scale the same way to put between 0 / 1
bpsp_proj = proj[:, :, :2].copy()
ma = bpsp_proj.max()
mi = bpsp_proj[:, :, :2].min()
ran = ma - mi
bpsp_proj += abs(mi)
bpsp_proj /= ran
baseline = abs(mi)
baseline /= ran

nbins = int(((stimulus.shape[1]/stim_fs) * 4) / len(epochs))
df['delta'] = (df['bp_dp'] - df['sp_dp']) / (df['bp_dp'] + df['sp_dp'])
# plot big pupil
bp_proj = np.stack([bpsp_proj[i, pup_mask[0, :, i], :] for i in range(pup_mask.shape[-1])])
im = chelp.plot_confusion_matrix(df, 
                    metric='bp_dp',
                    spectrogram=np.sqrt(stimulus)**(1/2),
                    sortby=('delta', nbins),
                    sort_method='1D',
                    resp_fs=4,
                    stim_fs=stim_fs,
                    pcs = bp_proj,
                    baseline=baseline,
                    vmin=0,
                    vmax=100,
                    pc_div=8,
                    ax=bp
                    )
bp.set_title(r"Large pupil $d'^2$")
cax = plt.subplot(gs[1, 1])
f.colorbar(im, ax=bp, cax=cax, orientation='horizontal', ticks=[0, 50, 100])
# plot small pupil
sp_proj = np.stack([bpsp_proj[i, ~pup_mask[0, :, i], :] for i in range(pup_mask.shape[-1])])
im = chelp.plot_confusion_matrix(df, 
                    metric='sp_dp',
                    spectrogram=np.sqrt(stimulus)**(1/2),
                    sortby=('delta', nbins),
                    sort_method='1D',
                    resp_fs=4,
                    stim_fs=stim_fs,
                    pcs = sp_proj,
                    baseline=baseline,
                    vmin=0,
                    vmax=100,
                    pc_div=pc_div,
                    ax=sp
                    )
sp.set_title(r"Small pupil $d'^2$")
cax = plt.subplot(gs[1, 4])
f.colorbar(im, ax=sp, cax=cax, orientation='horizontal', ticks=[0, 50, 100])

# plot difference
im = chelp.plot_confusion_matrix(df, 
                    metric='delta',
                    spectrogram=np.sqrt(stimulus)**(1/2),
                    sortby=('delta', nbins),
                    sort_method='1D',
                    resp_fs=4,
                    stim_fs=stim_fs,
                    vmin=-1,
                    vmax=1,
                    pc_div=pc_div,
                    ax=diff
                    )
diff.set_title(r"$\Delta d'^2$")
#cax = f.add_axes([0.1, 0.1, 0.1, 0.05])
cax = plt.subplot(gs[1, 7])
cbar = f.colorbar(im, ax=diff, cax=cax, orientation='horizontal', ticks=[-1, 0, 1])

# plot scatter plot of delta dprime results
# plot dprime results
nSamples = 3000
idx = df_all[['bp_dp', 'sp_dp']].max(axis=1) < 100
if idx.sum()<nSamples:
    nSamples = idx.sum()
sidx = np.random.choice(range(idx.sum()), nSamples, replace=False)
bp = df_all['bp_dp'].values[idx][sidx]
sp = df_all['sp_dp'].values[idx][sidx]
s = 5
xy = np.vstack([bp, sp])
z = gaussian_kde(xy)(xy)
scax.scatter(sp, bp, s=s, c=z, cmap='inferno')
scax.plot([0, 100], [0, 100], 'k--')
scax.set_xlabel("Small pupil")
scax.set_ylabel("Large pupil")
scax.set_title(r"Stimulus discriminability ($d'^2$)")
scax.axis('square')

# get statistics for all data
df_all['delta'] = (df_all['bp_dp'] - df_all['sp_dp']) #/ (df_all['bp_dp'] + df_all['sp_dp'])
na_mask = np.isnan(df_all["delta"])==False
df_all = df_all[na_mask]
df_all["delta_norm"] = df_all["delta"] / (df_all['bp_dp'] + df_all['sp_dp'])
d = {s: df_all[df_all.site==s]['delta'].values for s in df_all.site.unique()}
bs = get_bootstrapped_sample(d, even_sample=False, nboot=1000)
p = get_direct_prob(bs, np.zeros(bs.shape[0]))[0]

print(f"mean large pupil d': {df_all['bp_dp'].mean()}, {df_all['bp_dp'].sem()}")
print(f"mean small pupil d': {df_all['sp_dp'].mean()}, {df_all['sp_dp'].sem()}")
print(f"pval (bootstrapped): {p}")
print(f"Mean n stimulus pairs per session: {np.mean([len(d[s]) for s in d.keys()])}, {np.std([len(d[s]) for s in d.keys()]) / np.sqrt(len(d.keys()))}")


frac = []
for s in d.keys():
    frac.append(np.sum(d[s]<0) / len(d[s]))
print(f"Fraction of stimlulus pairs with decreases per site: {np.mean(frac)}, sem: {np.std(frac)/len(frac)}")

#f.tight_layout()

if savefig:
    f.savefig(fig_fn)
    f.savefig(fig_fn.replace(".svg", ".png"))


################################################# SUPPLEMENTAL FIGS ##########################################################
# mainly control analyses
# S1 - delta dprime does not depend on absolute dprime

lg_pos_fract = []
lg_neg_fract = []
sm_pos_fract = []
sm_neg_fract = []
f, ax = plt.subplots(1, 1, figsize=(2.5, 3))
for s in df_all.site.unique():
    dfr = df_all[df_all.site==s].reset_index()
    if dfr.shape[0] >= 20:
        mask = dfr['dp_opt_test'] > np.median(dfr['dp_opt_test'])
        lg_pos = (dfr[mask]['delta'] > 0).sum() / mask.sum(); lg_pos_fract.append(lg_pos)
        lg_neg = (dfr[mask]['delta'] < 0).sum() / mask.sum(); lg_neg_fract.append(lg_neg)
        sm_pos = (dfr[~mask]['delta'] > 0).sum() / (~mask).sum(); sm_pos_fract.append(sm_pos)
        sm_neg = (dfr[~mask]['delta'] < 0).sum() / (~mask).sum(); sm_neg_fract.append(sm_neg)
        if s in CPN_SITES:
            ax.plot([0, 1], [sm_pos, lg_pos], color='tab:blue')
        else:
            ax.plot([0, 1], [sm_pos, lg_pos], color='tab:orange')
ax.bar([0, 1], [np.mean(sm_pos_fract), np.mean(lg_pos_fract)], lw=2, color='none', edgecolor='k')
ax.set_ylabel("Prop. stim. pairs where "+r"$\Delta d'^2>0$")
ax.set_xticks([0, 1])
ax.set_xticklabels(["Small", "Large"])
ax.set_xlabel(r"Baseline discriminability ($d'^2$)")

stat, p = ss.wilcoxon(lg_pos_fract, sm_pos_fract)
ax.text(0.1, ax.get_ylim()[1], r"p=%s"%round(p, 3))

f.tight_layout()
if savefig:
    f.savefig(figS1)
    f.savefig(figS1.replace(".svg", ".png"))


# S2 - box plots - delta dprime for each site
boxprops = {"color": "k", "linewidth": 2}
f, ax = plt.subplots(1, 1, figsize=(8, 3))
bp = ax.boxplot([df_all[(df_all.site==s) & (df_all.batch==b)]["delta_norm"] for s, b in zip(sites, batches)], 
                showfliers=False,
                boxprops=boxprops,
                medianprops=boxprops,
                whiskerprops=boxprops,
                capprops=boxprops)
ax.axhline(0, linestyle="--", color="grey", zorder=-1)
ax.set_xlabel("Site")
ax.set_ylabel(r"$\Delta d'^2$")

if savefig:
    f.savefig(figS2)
    f.savefig(figS2.replace(".svg", ".png"))
