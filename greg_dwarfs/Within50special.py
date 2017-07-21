



def ngreater_specialLMC(self):
    # Component 1: MW background
    halo_mass = 1.4*10**12
    ratio_MW = 0.022158  # ratio of halos in virial radius to halos in current volume
    samples = self.generate_stellar_mass_samples(halo_mass,N=N_iter, factor=ratio_MW)


    # Component 2: SMC background
    


    N_visible, N_std, per10, per90 = Nvisible_From_Ngreater(min_lum, lums)
    return N_visible, N_std, per10, per90



