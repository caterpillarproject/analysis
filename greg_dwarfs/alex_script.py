import sys
sys.path.append("/nfs/blank/h4231/gdooley/analysis-caterpillar")

from testidea import convert_Mvir
import abundance_matching as am

mw = True

N_visibles = []
for model in [am.GK16_grow(hostSHMF=mw), am.GarrisonKimmel(hostSHMF=mw), am.Moster(hostSHMF=mw), am.Brook(hostSHMF=mw),am.GarrisonKimmel(hostSHMF=mw,reionization=True), am.GK16_grow(hostSHMF=mw,reionization=True),am.Moster(hostSHMF=mw,reionization=True), am.Brook(hostSHMF=mw,reionization=True)]:
    halo_mass = 1.4 * 10**12
    halo_mass = convert_Mvir(halo_mass,model)
    N_visible = model.ngreater(halo_mass, np.array([1000., 2.e5]))[0]
    print N_visible, N_visible[0] - N_visible[1]
    N_visibles.append(N_visible)

