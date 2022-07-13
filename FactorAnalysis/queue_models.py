import nems.db as nd
from global_settings import HIGHR_SITES, CPN_SITES

batches = [322] #[294, 322, 331]
modelnames = ["factor_analysis_pca_evoked_shuff"]
force_rerun = False

for batch in batches:

    if batch == 322:
        sites = [s for s in HIGHR_SITES if s not in ['BOL005c', 'BOL006b']]  
    if batch == 294:
        sites = ['BOL005c', 'BOL006b']
    if batch == 331:
        # just fit them all to have it done
        sites, _ = nd.get_batch_sites(331) #CPN_SITES

    script = '/auto/users/hellerc/code/projects/nat-ms-final/FactorAnalysis/cache_population_stats.py'
    python_path = '/auto/users/hellerc/miniconda3/envs/lbhb/bin/python'
    nd.enqueue_models(celllist=sites,
                    batch=batch,
                    modellist=modelnames,
                    executable_path=python_path,
                    script_path=script,
                    user='hellerc',
                    force_rerun=force_rerun,
                    reserve_gb=2,
                    priority=3)