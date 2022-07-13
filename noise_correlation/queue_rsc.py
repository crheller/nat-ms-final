import nems.db as nd
from global_settings import HIGHR_SITES, CPN_SITES

batches = [294, 322, 331] #[289, 294, 331]
batches = [322, 331]
force_rerun = True
boxcar = True
evoked = True
slow = False
perstim = False

xforms_modelname = "psth.fs4.pup.fpca3-ld-st.pup0.fpc-epcpn-hrc-psthfr_stategain.S_jk.nf10-basic"
xforms_modelname = "psth.fs4.pup.fpca3-ld-st.pup.fpc-epcpn-hrc-psthfr-aev_sdexp2.S_basic"
xforms_modelname = xforms_modelname.replace("_", "*")


modelnames = ['rsc', 'rsc_pr_rm2', 
            'rsc_fft0-0.05', 'rsc_pr_rm2_fft0-0.05',
            'rsc_fft0.05-0.1', 'rsc_pr_rm2_fft0.05-0.1',
            'rsc_fft0.1-0.25', 'rsc_pr_rm2_fft0.1-0.25',
            'rsc_fft0.25-0.5', 'rsc_pr_rm2_fft0.25-0.5',
            'rsc_fft0-0.1', 'rsc_pr_rm2_fft0-0.1',
            'rsc_fft0-0.25', 'rsc_pr_rm2_fft0-0.25',
            'rsc_fft0-0.5', 'rsc_pr_rm2_fft0-0.5',
            'rsc_fft0.5-2', 'rsc_pr_rm2_fft0.5-2',
            'rsc_fft2-4', 'rsc_pr_rm2_fft2-4',
            'rsc_fft0.1-4', 'rsc_pr_rm2_fft0.1-4',
            'rsc_fft0.25-4', 'rsc_pr_rm2_fft0.25-4',
            'rsc_fft0.5-4', 'rsc_pr_rm2_fft0.5-4',
            'rsc_fft4-10', 'rsc_pr_rm2_fft4-10',
            'rsc_fft10-25', 'rsc_pr_rm2_fft10-25',
            'rsc_fft25-50', 'rsc_pr_rm2_fft25-50']

if slow:
    modelnames = ['rsc', 'rsc_pr_rm2', 
            'rsc_fft0-0.05', 'rsc_pr_rm2_fft0-0.05',
            'rsc_fft0.05-0.1', 'rsc_pr_rm2_fft0.05-0.1',
            'rsc_fft0.1-0.25', 'rsc_pr_rm2_fft0.1-0.25',
            'rsc_fft0.25-0.5', 'rsc_pr_rm2_fft0.25-0.5',
            'rsc_fft0-0.1', 'rsc_pr_rm2_fft0-0.1',
            'rsc_fft0-0.25', 'rsc_pr_rm2_fft0-0.25',
            'rsc_fft0-0.5', 'rsc_pr_rm2_fft0-0.5',
            'rsc_fft0.5-2', 'rsc_pr_rm2_fft0.5-2']
            
if boxcar:
    modelnames = [m.replace('fft', 'boxcar_fft') for m in modelnames]

if evoked:
    modelnames = [m.replace('rsc', 'rsc_ev') for m in modelnames]

if slow:
    modelnames = [m.replace('fft', 'fs4_fft') for m in modelnames]

if perstim:
    modelnames = [m+'_perstim' for m in modelnames]

# remove old pupil regression / fft stuff
modelnames = [m for m in modelnames if ('fft' not in m) & ('pr' not in m)]

for batch in batches:

    if batch == 322:
        sites = [s for s in HIGHR_SITES if s not in ['BOL005c', 'BOL006b']]  
        xforms_model = xforms_modelname.replace("-epcpn-hrc", "-hrc")
    if batch == 294:
        sites = ['BOL005c', 'BOL006b']
        raise NotImplementedError("Haven't fit the xforms models yet, I don't think")
    if batch == 331:
        sites = CPN_SITES
        xforms_model = xforms_modelname

    xf_models = [m+f"_sub{xforms_model}" for m in modelnames]
    models = modelnames + xf_models

    script = '/auto/users/hellerc/code/projects/nat-ms-final/noise_correlation/cache_rsc.py'
    python_path = '/auto/users/hellerc/miniconda3/envs/lbhb/bin/python'
    nd.enqueue_models(celllist=sites,
                    batch=batch,
                    modellist=models,
                    executable_path=python_path,
                    script_path=script,
                    user='hellerc',
                    force_rerun=force_rerun,
                    reserve_gb=2)
