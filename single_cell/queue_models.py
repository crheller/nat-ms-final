import nems.db as nd
from global_settings import CPN_SITES

batch = 294  # 294 / 322
force_rerun = True

script = '/auto/users/hellerc/code/projects/nat-ms-final/single_cell/fit_script.py'
python_path = '/auto/users/hellerc/anaconda3/envs/lbhb/bin/python'

modelnames = ['ns.fs4.pup-ld-st.pup-hrc-psthfr_stategain.SxR_jk.nf10-basic',
              'ns.fs4.pup-ld-st.pup0-hrc-psthfr_stategain.SxR_jk.nf10-basic']
modelnames = ['ns.fs4.pup-ld-st.pup-hrc-psthfr_sdexp.SxR.bound_jk.nf10-basic',
              'ns.fs4.pup-ld-st.pup0-hrc-psthfr_sdexp.SxR.bound_jk.nf10-basic']

if batch == 294:
    modelnames = [m.replace('pup-ld', 'pup.voc-ld') for m in modelnames]

if batch == 331:
    modelnames = [m.replace('-hrc', '-epcpn-hrc') for m in modelnames]

cellids = nd.get_batch_cells(batch).cellid.tolist()

if batch == 331:
    cellids = [c for c in cellids if c[:7] in CPN_SITES]

if batch == 294:
    cellids = [c for c in cellids if c.split('-')[0] in ['BOL005c', 'BOL006b']]

nd.enqueue_models(celllist=cellids,
                  batch=batch,
                  modellist=modelnames,
                  executable_path=python_path,
                  script_path=script,
                  user='hellerc',
                  force_rerun=force_rerun,
                  reserve_gb=1)