"""
Compute noise correlations for the following conditions:
    all raw data
    all data pupil regressed
    all data lv regressed (not implemented)

    all of the above for big/small pupil conditions
"""

import nems.db as nd
import nems.xforms as xforms
import nems.xform_helper as xhelp
import nems_lbhb.baphy as nb
from nems_lbhb.preprocessing import mask_high_repetion_stims, create_pupil_mask, fix_cpn_epochs, movement_mask
import nems_lbhb.baphy_io as io
from nems_lbhb.baphy_experiment import BAPHYExperiment
from nems.recording import Recording
import sys
sys.path.append('/auto/users/hellerc/code/projects/nat_pupil_ms/')
sys.path.append('/auto/users/hellerc/code/projects/nat_pupil_ms/dprime/')
import load_results as ld
sys.path.append('/auto/users/hellerc/code/charlieTools/')
import charlieTools.simulate_data as sim
import charlieTools.preprocessing as preproc
import charlieTools.noise_correlations as nc
from charlieTools.nat_sounds_ms.preprocessing import get_pupil_range
from charlieTools.nat_sounds_ms.decoding import load_site
import pandas as pd
import numpy as np
import scipy.stats as ss
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from itertools import combinations
import os
import nems
import copy
import logging

log = logging.getLogger(__name__)

if 'QUEUEID' in os.environ:
    queueid = os.environ['QUEUEID']
    nems.utils.progress_fun = nd.update_job_tick

else:
    queueid = 0

if queueid:
    log.info("Starting QUEUEID={}".format(queueid))
    nd.update_job_start(queueid)

# first system argument is the cellid
site = sys.argv[1]  # very first (index 0) is the script to be run
# second systems argument is the batch
batch = sys.argv[2]
# third system argument in the modelname
modelname = sys.argv[3]

perstim = 'perstim' in modelname
balanced = 'bal' in modelname
filt = 'fft' in modelname
if balanced:
    raise DeprecationWarning("Removed pupil balancing option (temporarily?) CRH 05/08/2020")

keys = modelname.split('_')
boxcar = False
evoked = False
fs4 = False
subtract_xforms_pred = False
for k in keys:
    if 'fft' in k:
        low_c = np.float(k.split('-')[0][3:])
        high_c = np.float(k.split('-')[1])
    if 'boxcar' in k:
        boxcar = True
    if 'ev' in k:
        evoked = True
    if 'fs4' in k:
        fs4 = True
    if k.startswith("sub"):
        subtract_xforms_pred = True
        xforms_modelname = k[3:].replace("*", "_")

path = '/auto/users/hellerc/results/nat_pupil_ms/noise_correlations_final/'

log.info('Computing noise correlations for site: {0} with options: \n \
            balanced pupil epochs: {1}'.format(site, balanced))

log.info("Saving results to: {}".format(path))


batch = int(batch)
if filt & (fs4 == False):
    fs = 100
    xforms_modelname = 'ns.fs100.pup-ld-st.pup-hrc-psthfr_sdexp.SxR.bound_jk.nf10-basic'
    #xforms_modelname = 'ns.fs100.pup-ld-st.pup-hrc-psthfr-ev_slogsig.SxR-lv.1xR-lvlogsig.2xR_jk.nf5.p-pupLVbasic.a0:0'
else:
    fs = 4

if batch == 294:
    xforms_modelname = xforms_modelname.replace('pup-ld', 'pup.voc-ld')

cellid, _ = io.parse_cellid({'batch': batch, 'cellid': site})

# just load raw data
if not subtract_xforms_pred:
    # only load model fit if using for regression
    if 0: #batch != 331:
        options = {'cellid': site, 'rasterfs': fs, 'batch': batch, 'pupil': True, 'stim': False}
        if batch == 294:
            options['runclass'] = 'VOC'
        rec = nb.baphy_load_recording_file(**options)

    else:
        manager = BAPHYExperiment(cellid=site, batch=batch)
        options = {'rasterfs': 4, 'resp': True, 'stim': False, 'pupil': True, 'pupil_variable_name': 'area'}
        rec = manager.get_recording(**options)
        rec['resp'] = rec['resp'].rasterize()
        if batch == 294:
            stims = [s for s in rec['resp'].epochs.name.unique() if ('STIM_' in s) & ('Pips' not in s)]
            rec = rec.and_mask(stims)
            rec = rec.apply_mask(reset_epochs=True)
            
    rec['resp'] = rec['resp'].rasterize()
    if 'cells_to_extract' in rec.meta.keys():
        if rec.meta['cells_to_extract'] is not None:
            log.info("Extracting cellids: {0}".format(rec.meta['cells_to_extract']))
            rec['resp'] = rec['resp'].extract_channels(rec.meta['cells_to_extract'])
    if batch == 331:
        rec = fix_cpn_epochs(rec)
    if (batch == 294) | (batch==331):
        epochs = [epoch for epoch in rec.epochs.name.unique() if 'STIM_' in epoch]
    else:
        epochs = [epoch for epoch in rec.epochs.name.unique() if 'STIM_00' in epoch]
    
    rec = rec.and_mask(epochs)

    rec = rec.apply_mask(reset_epochs=True)

else:
    # load xforms model, subtract prediction
    log.info("Load recording from xforms model {}".format(xforms_modelname))
    # 'psth' signal returned in this function is 'pred'
    rec = preproc.generate_state_corrected_psth(batch=batch, modelname=xforms_modelname, cellids=cellid, 
                                        siteid=site, recache=False)

# filtering / pupil regression must always go first!
if subtract_xforms_pred:
    log.info('Regress state by subtracting xforms model pred')
    mod_data = rec['resp']._data - rec['psth']._data + rec['psth_sp']._data
    rec['resp'] = rec['resp']._modified_copy(mod_data)

if filt:
    log.info("Band-pass filter spike counts between {0} and {1} Hz".format(low_c, high_c))
    # before filtering, pad nans with 0 (these are non-validation stimuli. eg. stims with only one rep)
    # for filtering to be approx in real time, need these periods too
    resp_data = rec['resp']._data.copy()
    nan_idx = np.isnan(resp_data[0, :])
    resp_data[:, nan_idx] = 0
    rec['resp'] = rec['resp']._modified_copy(resp_data)

    rec = preproc.bandpass_filter_resp(rec, low_c, high_c, boxcar=boxcar)

# also mask evoked periods only ?
if evoked:
    # double check that mask is cast to bool
    if 'mask' not in rec.signals:
        rec = rec.create_mask(True)
    rec['mask'] = rec['mask']._modified_copy(rec['mask']._data.astype(bool))
    rec = rec.and_mask(['PreStimSilence', 'PostStimSilence'], invert=True)
    if batch != 331:
        nspont_bins = rec['resp'].extract_epoch('PreStimSilence').shape[-1]
    else:
        nspont_bins = 0
else:
    nspont_bins = 0
rec = rec.apply_mask(reset_epochs=True)

# get the pupil range for each stimulus pair, then save the mean of this for the site
X, sp_bins, X_pup, pup_mask = load_site(site=site, batch=batch)
ncells = X.shape[0]
nreps_raw = X.shape[1]
nstim = X.shape[2]
nbins = X.shape[3]
sp_bins = sp_bins.reshape(1, sp_bins.shape[1], nstim * nbins)
nstim = nstim * nbins
X_pup = X_pup.reshape(1, nreps_raw, nstim)
pup_mask = pup_mask.reshape(1, nreps_raw, nstim)
pupil_range = get_pupil_range(X_pup, pup_mask)
mean_pupil_range = pupil_range['range'].iloc[:-1].mean()

log.info("Mask large and small pupil")
# handle using xforms pupil mask
rec_bp = rec.copy()
ops = {'state': 'big', 'method': 'median', 'epoch': ['REFERENCE'], 'collapse': True}
rec_bp = create_pupil_mask(rec_bp, **ops)
ops = {'state': 'small', 'method': 'median', 'epoch': ['REFERENCE'], 'collapse': True}
rec_sp = rec.copy()
rec_sp = create_pupil_mask(rec_sp, **ops)

eps = np.unique([s for s in rec.epochs.name if 'STIM' in s]).tolist()

log.info("Extracting spike count dictionaries for big pupil, \n \
    small pupil, and all pupil trials.")
real_dict_all = rec['resp'].extract_epochs(eps)
real_dict_small = rec_sp['resp'].extract_epochs(eps, mask=rec_sp['mask'])
real_dict_big = rec_bp['resp'].extract_epochs(eps, mask=rec_bp['mask'])

if perstim:
    log.info("Computing noise correlations separately for each stimulus bin")
    for idx, stim in enumerate(eps):
        log.info(f"Epoch {idx+1} / {len(eps)}")
        for b in range(real_dict_all[stim].shape[-1]):
            _df_all = nc.compute_rsc({stim: real_dict_all[stim][:, :, [b]]}, chans=rec['resp'].chans)
            _df_big = nc.compute_rsc({stim: real_dict_big[stim][:, :, [b]]}, chans=rec['resp'].chans)
            _df_small = nc.compute_rsc({stim: real_dict_small[stim][:, :, [b]]}, chans=rec['resp'].chans)      
            _df_all['stim'] = stim+'_'+str(b+nspont_bins)
            _df_big['stim'] = stim+'_'+str(b+nspont_bins)
            _df_small['stim'] = stim+'_'+str(b+nspont_bins)
            if (idx==0) & (b==0):
                df_all = _df_all
                df_big = _df_big
                df_small = _df_small      
            else:
                df_all = pd.concat([df_all, _df_all])
                df_big = pd.concat([df_big, _df_big])
                df_small = pd.concat([df_small, _df_small])

else:
    df_all = nc.compute_rsc(real_dict_all, chans=rec['resp'].chans)
    df_big = nc.compute_rsc(real_dict_big, chans=rec['resp'].chans)
    df_small = nc.compute_rsc(real_dict_small, chans=rec['resp'].chans)
    df_all['stim'] = 'all'
    df_big['stim'] = 'all'
    df_small['stim'] = 'all'

cols = ['all', 'p_all', 'gm_all', 'bp', 'p_bp', 'gm_bp', 'sp', 'p_sp', 'gm_sp', 'site', 'stim']
df = pd.DataFrame(columns=cols, index=df_all.index)
df['all'] = df_all['rsc']
df['p_all'] = df_all['pval']
df['gm_all'] = df_all['gmean']
df['bp'] = df_big['rsc']
df['p_bp'] = df_big['pval']
df['gm_bp'] = df_big['gmean']
df['sp'] = df_small['rsc']
df['p_sp'] = df_small['pval']
df['gm_sp'] = df_small['gmean']
df['site'] = site
df['stim'] = df_all['stim']
df['mean_pupil_range'] = mean_pupil_range

dtypes = {
    'all': 'float64',
    'bp': 'float64',
    'sp': 'float64',
    'p_all': 'float64',
    'p_bp': 'float64',
    'p_sp': 'float64',
    'site': 'object',
    'stim': 'object',
    'mean_pupil_range': 'float64'
}
df = df.astype(dtypes)

log.info("save noise corr results for model: {0}, batch: {1}, site: {2}".format(modelname, batch, site))

if os.path.isdir(os.path.join(path, str(batch), site)):
    df.to_csv(os.path.join(path, str(batch), site, modelname+'.csv'))
elif os.path.isdir(os.path.join(path, str(batch))):
    os.mkdir(os.path.join(path, str(batch), site))
    df.to_csv(os.path.join(path, str(batch), site, modelname+'.csv'))
else:
    os.mkdir(os.path.join(path, str(batch)))
    os.mkdir(os.path.join(path, str(batch), site))
    df.to_csv(os.path.join(path, str(batch), site, modelname+'.csv'))

if queueid:
    nd.update_job_complete(queueid)
