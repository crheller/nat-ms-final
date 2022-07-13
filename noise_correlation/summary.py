import load_results as ld
import matplotlib.pyplot as plt
import pandas as pd

rsc_path = '/auto/users/hellerc/results/nat_pupil_ms/noise_correlations_final/'

rsc = "rsc_ev"
modelname331 = "psth.fs4.pup.fpca3-ld-st.pup0.fpc-epcpn-hrc-psthfr*stategain.S*jk.nf10-basic"
modelname322 = "psth.fs4.pup.fpca3-ld-st.pup0.fpc-hrc-psthfr*stategain.S*jk.nf10-basic"

modelname331 = "psth.fs4.pup.fpca3-ld-st.pup.fpc-epcpn-hrc-psthfr-aev*sdexp2.S*basic"
modelname322 = "psth.fs4.pup.fpca3-ld-st.pup.fpc-hrc-psthfr-aev*sdexp2.S*basic"

rsc_df331 = ld.load_noise_correlation(rsc, xforms_model=modelname331, batch=331, path=rsc_path)
rsc_df322 = ld.load_noise_correlation(rsc, xforms_model=modelname322, batch=322, path=rsc_path)

rsc_df = pd.concat([rsc_df331, rsc_df322])
rg = rsc_df.groupby(by=["site", "batch"]).mean()

f, ax = plt.subplots(1, 3, figsize=(12, 4), sharex=True, sharey=True)

ax[0].scatter(rg["bp"], rg["sp"], s=25)

m = rg.index.get_level_values(1) == 331
ax[1].scatter(rg[m]["bp"], rg[m]["sp"], s=25)
ax[1].set_title("batch 331")

m = rg.index.get_level_values(1) == 322
ax[2].scatter(rg[m]["bp"], rg[m]["sp"], s=25)
ax[2].set_title("batch 322")

for a in ax:
    a.plot([-1, 1], [-1, 1], "k--")
    a.set_xlabel("Large pupil")
    a.set_ylabel("Small pupil")


ax[2].set_xlim((-0, 0.15))
ax[2].set_ylim((-0, 0.15))

f.tight_layout()
