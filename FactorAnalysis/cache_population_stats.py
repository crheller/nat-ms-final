"""
Inspired by Umakanthan 2021 neuron paper -- bridging neuronal correlations / dim reduction
Cache script --  
for each site, calculate:
    - % shared variance
    - loading similarity
    - dimensionality
    - also compute dimensionality / components of full space, including stimuli, using PCA (trial-averaged or raw data?)
        - think raw data is okay, but we'll do both. The point we want to make is that variability coflucates in a space with dim < total dim
"""
import numpy as np
from sklearn.decomposition import FactorAnalysis, PCA
from scipy.signal import argrelextrema
from charlieTools.nat_sounds_ms.decoding import load_site
import pickle


import os
import sys
import nems.db as nd
import logging 
import nems

log = logging.getLogger(__name__)

if 'QUEUEID' in os.environ:
    queueid = os.environ['QUEUEID']
    nems.utils.progress_fun = nd.update_job_tick

else:
    queueid = 0

if queueid:
    log.info("Starting QUEUEID={}".format(queueid))
    nd.update_job_start(queueid)


# Get sys args
site = sys.argv[1]  
batch = sys.argv[2]
modelname = sys.argv[3]

evoked = False
shuffle = False
for op in modelname.split("_"):
    if op == "evoked":
        evoked = True
    if op == "shuff":
        shuffle = True

# measure change in dimensionality, %sv, loading sim, across jackknifes
def get_dim(LL):
    try:
        return argrelextrema(LL, np.greater)[0][0]+1
    except:
        log.info("No relative LL max, choosing overall maximum")
        return np.argmax(LL)+1

def sigma_shared(model):
    return (model.components_.T @ model.components_)

def sigma_ind(model):
    return np.diag(model.noise_variance_)

def pred_cov(model):
    return sigma_shared(model) + sigma_ind(model)

def get_dim95(model):
    """
    number of dims to explain 95% of shared var
    """
    ss = sigma_shared(model)
    evals, _ = np.linalg.eig(ss)
    evals = evals / sum(evals)
    return np.argwhere(np.cumsum(evals)>=0.95)[0][0]+1

def get_sv(model):
    sig_shared = sigma_shared(model) # rank n_components cov matrix
    full_cov_pred = pred_cov(model)
    # % shared variance
    # per neuron
    pn = np.diag(sig_shared) / np.diag(full_cov_pred)
    # average
    sv = np.mean(pn)
    return sv

def get_loading_similarity(model, dim=0):
    # loading similarity
    loading = model.components_[dim, :]
    loading /= np.linalg.norm(loading)
    load_sim = 1 - (np.var(loading) / (1 / len(loading)))
    return load_sim


# ============ perform analysis ==================

# load data
X, sp_bins, X_pup, pup_mask, epochs = load_site(site=site, batch=int(batch), regress_pupil=False, use_xforms=False, return_epoch_list=True)

if evoked:
    log.info("Extracting only evoked activity")
    X = X[:, :, :, sp_bins[0, 0, 0, :]==False]
    pup_mask = pup_mask[:, :, :, sp_bins[0, 0, 0, :]==False]

# fit all stim together, after subtracting psth
# "special" cross-validation -- fitting individual stims doesn't work, not enough data
# instead, leave-one-stim out fitting to find dims that are shared / stimulus-independent
nstim = X.shape[-1] * X.shape[-2]
nCells = X.shape[0]
Xsub = (X - X.mean(axis=1, keepdims=True))
Xfa = Xsub.reshape(X.shape[0], X.shape[1], nstim)
Xpca = X.reshape(X.shape[0], X.shape[1], nstim)
Xpca_ta = X.mean(axis=1, keepdims=True).reshape(X.shape[0], 1, nstim)

pm = pup_mask.reshape(pup_mask.shape[0], pup_mask.shape[1], nstim)
if shuffle:
    pm = pm[:, np.random.choice(range(pup_mask.shape[1]), pup_mask.shape[1], replace=False), :]

# large / small trial-average matrices
Xpca_ta_large = np.concatenate([Xpca[:, pm[0, :, i]==True, [i]].mean(axis=1, keepdims=True) for i in range(Xpca.shape[-1])], axis=-1)[:, np.newaxis, :]
Xpca_ta_small = np.concatenate([Xpca[:, pm[0, :, i]==False, [i]].mean(axis=1, keepdims=True) for i in range(Xpca.shape[-1])], axis=-1)[:, np.newaxis, :]

nfold = nstim
nComponents = 50
if X.shape[0] < nComponents:
    nComponents = X.shape[0]

if X.shape[0] > nstim:
    smallestDim = nstim
else:
    smallestDim = X.shape[0]

log.info("\nComputing log-likelihood across models / nfolds")
LL = np.zeros((nComponents, nfold))
LL_small = np.zeros((nComponents, nfold))
LL_large = np.zeros((nComponents, nfold))
LL_pca_all = np.zeros((nComponents, nfold))
LL_pca_ta_all = np.zeros((nComponents, nfold))
LL_pca_small = np.zeros((nComponents, nfold))
LL_pca_large = np.zeros((nComponents, nfold))
LL_pca_ta_small = np.zeros((nComponents, nfold))
LL_pca_ta_large = np.zeros((nComponents, nfold))
for ii in np.arange(1, LL.shape[0]+1):
    log.info(f"{ii} / {LL.shape[0]}")
    fa = FactorAnalysis(n_components=ii, random_state=0) # init model
    pca = PCA(n_components=ii)
    for nf in range(nfold):
        fit = [x for x in np.arange(0, nstim) if x != nf]

        # fit all data
        fa.fit(Xfa[:, :, fit].reshape(Xfa.shape[0], -1).T) # fit model
        # Get LL score
        LL[ii-1, nf] = fa.score(Xfa[:, :, nf].T)

        # fit large pupil
        fa.fit(Xfa[:, :, fit][:, pm[0, :, fit].T].reshape(Xfa.shape[0], -1).T) # fit model 
        LL_large[ii-1, nf] = fa.score(Xfa[:, pm[0, :, nf], nf].T)

        # fit small pupil
        fa.fit(Xfa[:, :, fit][:, pm[0, :, fit].T==False].reshape(Xfa.shape[0], -1).T) # fit model 
        LL_small[ii-1, nf] = fa.score(Xfa[:, pm[0, :, nf]==False, nf].T)

        # PCA
        try:
            pca.fit(Xpca_ta[:, 0, fit].T)
            LL_pca_ta_all[ii-1, nf] = pca.score(Xpca_ta[:, 0, [nf]].T)

            pca.fit(Xpca_ta_large[:, 0, fit].T)
            LL_pca_ta_large[ii-1, nf] = pca.score(Xpca_ta_large[:, 0, [nf]].T)

            pca.fit(Xpca_ta_small[:, 0, fit].T)
            LL_pca_ta_small[ii-1, nf] = pca.score(Xpca_ta_small[:, 0, [nf]].T)
        except:
            LL_pca_ta_all[ii-1, nf] = np.nan
            LL_pca_ta_large[ii-1, nf] = np.nan
            LL_pca_ta_small[ii-1, nf] = np.nan

        # CV different for full data pca
        _xpca_flat = Xpca.reshape(nCells, -1)
        _pmask_flat = pm[0, :, :].reshape(1, -1)
        _xpca_lrg = _xpca_flat[:, _pmask_flat[0]==True]
        _xpca_sm = _xpca_flat[:, _pmask_flat[0]==False]

        nval = int(_xpca_lrg.shape[-1] / nfold)
        val = np.arange(nf*nval, (nf+1)*nval)
        fit = np.array([v for v in np.arange(0, _xpca_lrg.shape[-1]) if v not in val])
        pca.fit(_xpca_lrg[:, fit].T)
        LL_pca_large[ii-1, nf] = pca.score(_xpca_lrg[:, val].T)

        nval = int(_xpca_sm.shape[-1] / nfold)
        val = np.arange(nf*nval, (nf+1)*nval)
        fit = np.array([v for v in np.arange(0, _xpca_sm.shape[-1]) if v not in val])
        pca.fit(_xpca_sm[:, fit].T)
        LL_pca_small[ii-1, nf] = pca.score(_xpca_sm[:, val].T)

        nval = int(_xpca_flat.shape[-1] / nfold)
        val = np.arange(nf*nval, (nf+1)*nval)
        fit = np.array([v for v in np.arange(0, _xpca_flat.shape[-1]) if v not in val])
        pca.fit(_xpca_flat[:, fit].T)
        LL_pca_all[ii-1, nf] = pca.score(_xpca_flat[:, val].T)


log.info("Estimating %sv and loading similarity for the 'best' model")
# all data
all_dim_sem = np.std([get_dim(LL[:, i]) for i in range(LL.shape[1])]) / np.sqrt(LL.shape[1])
all_dim = get_dim(LL.mean(axis=-1))
# fit the "best" model over jackknifes
all_sv = np.zeros(nfold)
all_loading_sim = np.zeros(nfold)
all_dim95 = np.zeros(nfold)
for nf in range(nfold):
    fit = [x for x in np.arange(0, nstim) if x != nf]
    fa_all = FactorAnalysis(n_components=all_dim, random_state=0) 
    fa_all.fit(Xfa[:, :, fit].reshape(nCells, -1).T)
    all_sv[nf] = get_sv(fa_all)
    all_loading_sim[nf] = get_loading_similarity(fa_all)
    # get n dims needs to explain 95% of shared variance
    all_dim95[nf] = get_dim95(fa_all)

# small pupil
sp_dim_sem = np.std([get_dim(LL_small[:, i]) for i in range(LL.shape[1])]) / np.sqrt(LL.shape[1])
sp_dim = get_dim(LL_small.mean(axis=-1))
# fit the "best" model over jackknifes
sp_sv = np.zeros(nfold)
sp_loading_sim = np.zeros(nfold)
sp_dim95 = np.zeros(nfold)
for nf in range(nfold):
    fit = [x for x in np.arange(0, nstim) if x != nf]
    fa_small = FactorAnalysis(n_components=sp_dim, random_state=0) 
    fa_small.fit(Xfa[:, :, fit][:, pm[0, :, fit].T==False].T)
    sp_sv[nf] = get_sv(fa_small)
    sp_loading_sim[nf] = get_loading_similarity(fa_small)
    # get n dims needs to explain 95% of shared variance
    sp_dim95[nf] = get_dim95(fa_small)

# large pupil
bp_dim_sem = np.std([get_dim(LL_large[:, i]) for i in range(LL.shape[1])]) / np.sqrt(LL.shape[1])
bp_dim = get_dim(LL_large.mean(axis=-1))
# fit the "best" model over jackknifes
bp_sv = np.zeros(nfold)
bp_loading_sim = np.zeros(nfold)
bp_dim95 = np.zeros(nfold)
for nf in range(nfold):
    fit = [x for x in np.arange(0, nstim) if x != nf]
    fa_big = FactorAnalysis(n_components=bp_dim, random_state=0) 
    fa_big.fit(Xfa[:, :, fit][:, pm[0, :, fit].T==True].T)
    bp_sv[nf] = get_sv(fa_big)
    bp_loading_sim[nf] = get_loading_similarity(fa_big)
    # get n dims needs to explain 95% of shared variance
    bp_dim95[nf] = get_dim95(fa_big)


# final fit with all data to get components
fa_all = FactorAnalysis(n_components=all_dim, random_state=0) 
fa_all.fit(Xfa.reshape(nCells, -1).T)
all_sv_all = get_sv(fa_all)
all_ls_all = get_loading_similarity(fa_all)
all_dim95_all = get_dim95(fa_all)

fa_big = FactorAnalysis(n_components=bp_dim, random_state=0) 
fa_big.fit(Xfa[:, pm[0]==True].T)
bp_sv_all = get_sv(fa_big)
bp_ls_all = get_loading_similarity(fa_big)
bp_dim95_all = get_dim95(fa_big)

fa_small = FactorAnalysis(n_components=sp_dim, random_state=0) 
fa_small.fit(Xfa[:, pm[0]==False].T)
sp_sv_all = get_sv(fa_small)
sp_ls_all = get_loading_similarity(fa_small)
sp_dim95_all = get_dim95(fa_small)



# PCA results
bp_pca_ta_dim_sem = np.std([get_dim(LL_pca_ta_large[:, i]) for i in range(LL.shape[1])]) / np.sqrt(LL.shape[1])
bp_pca_ta_dim = get_dim(LL_pca_ta_large.mean(axis=-1))
sp_pca_ta_dim_sem = np.std([get_dim(LL_pca_ta_small[:, i]) for i in range(LL.shape[1])]) / np.sqrt(LL.shape[1])
sp_pca_ta_dim = get_dim(LL_pca_ta_small.mean(axis=-1))
bp_pca_dim_sem = np.std([get_dim(LL_pca_large[:, i]) for i in range(LL.shape[1])]) / np.sqrt(LL.shape[1])
bp_pca_dim = get_dim(LL_pca_large.mean(axis=-1))
sp_pca_dim_sem = np.std([get_dim(LL_pca_small[:, i]) for i in range(LL.shape[1])]) / np.sqrt(LL.shape[1])
sp_pca_dim = get_dim(LL_pca_small.mean(axis=-1))

pca_ta_dim_sem = np.std([get_dim(LL_pca_ta_all[:, i]) for i in range(LL.shape[1])]) / np.sqrt(LL.shape[1])
pca_ta_dim = get_dim(LL_pca_ta_all.mean(axis=-1))
pca_dim_sem = np.std([get_dim(LL_pca_all[:, i]) for i in range(LL.shape[1])]) / np.sqrt(LL.shape[1])
pca_dim = get_dim(LL_pca_all.mean(axis=-1))

pca_big = PCA(n_components=bp_pca_ta_dim) 
pca_big.fit(Xpca_ta_large[:, 0, :].T)

pca_small = PCA(n_components=sp_pca_ta_dim) 
pca_small.fit(Xpca_ta_small[:, 0, :].T)

# Save results
results = {
    "all_sv": all_sv.mean(),
    "all_sv_sd": all_sv.std(),
    "all_loading_sim": all_loading_sim.mean(),
    "all_loading_sim_sd": all_loading_sim.std(),
    "all_dim95": all_dim95.mean(),
    "all_dim95_sd": all_dim95.std(),
    "all_dim": all_dim,
    "all_dim_sem": all_dim_sem,
    "bp_sv": bp_sv.mean(),
    "sp_sv": sp_sv.mean(),
    "bp_sv_sd": bp_sv.std(),
    "sp_sv_sd": sp_sv.std(),
    "bp_loading_sim": bp_loading_sim.mean(),
    "sp_loading_sim": sp_loading_sim.mean(),
    "bp_loading_sim_sd": bp_loading_sim.std(),
    "sp_loading_sim_sd": sp_loading_sim.std(),
    "bp_dim95": bp_dim95.mean(),
    "sp_dim95": sp_dim95.mean(),
    "bp_dim95_sd": bp_dim95.std(),
    "sp_dim95_sd": sp_dim95.std(),
    "bp_dim": bp_dim,
    "sp_dim": sp_dim,
    "bp_dim_sem": bp_dim_sem,
    "sp_dim_sem": sp_dim_sem,
    "nCells": X.shape[0],
    "nStim": nstim,
    "final_fit": {
        "fa_all.components_": fa_all.components_,
        "fa_all.sigma_shared": sigma_shared(fa_all),
        "fa_all.sigma_ind": sigma_ind(fa_all),
        "fa_all.sigma_full": pred_cov(fa_all),
        "all_sv_all": all_sv_all,
        "all_ls_all": all_ls_all,
        "all_dim95_all": all_dim95_all,
        "fa_big.components_": fa_big.components_,
        "fa_small.components_": fa_small.components_,
        "fa_big.sigma_shared": sigma_shared(fa_big),
        "fa_small.sigma_shared": sigma_shared(fa_small),
        "fa_big.sigma_ind": sigma_ind(fa_big),
        "fa_small.sigma_ind": sigma_ind(fa_small),
        "fa_big.sigma_full": pred_cov(fa_big),
        "fa_small.sigma_full": pred_cov(fa_small),
        "bp_sv_all": bp_sv_all,
        "sp_sv_all": sp_sv_all,
        "bp_ls_all": bp_ls_all,
        "sp_ls_all": sp_ls_all,
        "bp_dim95_all": bp_dim95_all,
        "sp_dim95_all": sp_dim95_all
    },
    "pca_results": {
        "bp_full_dim": bp_pca_dim,
        "bp_full_dim_sem": bp_pca_dim_sem,
        "sp_full_dim": sp_pca_dim,
        "sp_full_dim_sem": sp_pca_dim_sem,
        "bp_ta_dim": bp_pca_ta_dim,
        "bp_ta_dim_sem": bp_pca_ta_dim_sem,
        "sp_ta_dim": sp_pca_ta_dim,
        "sp_ta_dim_sem": sp_pca_ta_dim_sem,
        "pca_ta_large.components_": pca_big.components_,
        "pca_ta_small.components_": pca_small.components_,
        "ta_dim": pca_ta_dim,
        "ta_dim_sem": pca_ta_dim_sem,
        "full_dim": pca_dim,
        "full_dim_sem": pca_dim_sem
    }    
}

def save(d, path):
    with open(path+f'/{modelname}.pickle', 'wb') as handle:
        pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return None

path = "/auto/users/hellerc/results/nat_pupil_ms/factor_analysis/"
if os.path.isdir(os.path.join(path, str(batch), site)):
   pass
elif os.path.isdir(os.path.join(path, str(batch))):
    os.mkdir(os.path.join(path, str(batch), site))
else:
    os.mkdir(os.path.join(path, str(batch)))
    os.mkdir(os.path.join(path, str(batch), site))

save(results, os.path.join(path, str(batch), site))

if queueid:
    nd.update_job_complete(queueid)