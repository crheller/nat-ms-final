# save dataframe of results to /auto/users/hellerc
import sys
sys.path.append("/auto/users/hellerc/code/projects/nat-ms-final")
from single_cell.mod_per_state import get_model_results_per_state_model
import os

batch = 331
d = get_model_results_per_state_model(batch=batch, state_list=['st.pup', 'st.pup0'], 
                                                  loader='ns.fs4.pup-ld-', 
                                                  fitter='_jk.nf10-basic',
                                                  basemodel='-epcpn-hrc-psthfr_sdexp.SxR.bound')
path = '/auto/users/hellerc/results/nat_pupil_ms/first_order_model_results/'
d.to_csv(os.path.join(path, "d_{0}_pup_sdexp.csv".format(batch)))