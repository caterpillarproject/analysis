import abundance_matching as am  # N_iter is specified
import DwarfIrregularTable as dit

import testidea  # no Niter specified here
import reionization # no Niter specified here

am.plotMstar_v_Mhalo()   # Figure 1
reionization.plot_frac_lum_mdef()  # Figure 2
testidea.ngreater_v_minlum(mw=True)   # Figure 3

   
testidea.plot_ngreater_2panel(min_lum=[10**3, 10**4],mstar_axes=True)
testidea.plot_P_at_least_one_2panel()

dit.plot_sum_distr_3panel(subset=True,re=True)  # Figure Final



