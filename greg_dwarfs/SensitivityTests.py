import numpy as np
import abundance_matching as am
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PlotParams import *

# How does the mean number of satellites above 10^4 change for Moster, GK16, and Brook for each of the following changes:


"""
1) reionization model: at 4 different redshifts
2) Moster - shift infall distribution by 0.5 Gyrs
3) GK16 - increase the sigma by a lot
4) change M* of field halos by 10% up and down
5) change slope of SHMF - force it to -1.9, -1.85 and get best fit K
6) magnitude of SHMF (no, why would it be uncertain?)
7) Change from Field SHMF to Host SHMF  (must use rmaxcut=False)
8) take SHMF of half of the halos and the other half. like jackknife?
9) host halo mass +/- 10%

"""

vmax_ach = 9.48535156 
vmax_filt = 23.54248047 


N_iter = 25000  #25000
StellarMasses = 10**6*np.array([0.14,0.54,0.77,0.82,1.3,1.4,1.6,2.1,2.7,3.5,3.5,4.5,6.0,6.4,6.6,7.8,8.3,16,16,17,19,37,44,47,51,52,62,76,100,100,190,270])
# my metric is the sum of the 5 largest systems with Moster
star_masses = StellarMasses[-5:]


def get_sum_fields(model, min_lum, star_masses):
    #dm_masses = model.stellar_to_halo_mass(star_masses,a=1.0)
    dm_masses = model.mstar_to_mhalo(star_masses,a=1.0)
    if model.isPoisson:
        mean = model.get_field_total(dm_masses,min_lum)
    else:
        samples = model.get_field_distr(dm_masses,min_lum,N=N_iter)
        mean = np.mean(samples)
    return mean



def get_baseline_moster(lum=4):  # not updated to mstar_to_mhalo
    if lum == 4:
        return  11.44
    if lum == 3:
        return 15.37
    if lum == 5:
        return 5.47
    """
    model = am.Moster(reionization = True)
    norm = get_sum_fields(model,10**lum,star_masses)
    print norm, 'norm for 10^4'
    norm = get_sum_fields(model,10**lum,star_masses)
    print norm, 'norm for 10^4'
    norm = get_sum_fields(model,10**3,star_masses)
    print norm, 'norm for 10^3'
    norm = get_sum_fields(model,10**3,star_masses)
    print norm, 'norm for 10^3'
   
    model = am.Moster(reionization = True)
    norm = get_sum_fields(model,10**5,star_masses)
    print norm, 'norm for 10^5'
    norm = get_sum_fields(model,10**5,star_masses)
    print norm, 'norm for 10^5'

    return norm
    """


def get_baseline_GK16(lum=4):   # not updated to mstar_to_mhalo. commented values are stellar_to_halo_mass()
    if lum ==5:
        return 6.79 #return 8.7
    if lum == 4:
        return 12.79 #return 16.34
    if lum == 3:
        return 17.39 #return 22.13
    """

    model = am.GK16_grow(reionization = True)
    norm = get_sum_fields(model,10**4,star_masses)
    print norm, 'norm for 10^4'
    norm = get_sum_fields(model,10**4,star_masses)
    print norm, 'norm for 10^4'
    norm = get_sum_fields(model,10**3,star_masses)
    print norm, 'norm for 10^3'
    norm = get_sum_fields(model,10**3,star_masses)
    print norm, 'norm for 10^3'

    norm = get_sum_fields(model,10**5,star_masses)
    print norm, 'norm for 10^5'
    norm = get_sum_fields(model,10**5,star_masses)
    print norm, 'norm for 10^5'


    return norm
    """
def get_baseline_Brook(lum=4):  # not updated to mstar_to_mhalo
    if lum == 5:
        return 2.71
    if lum == 4:
        return 5.07
    if lum == 3:
        return 8.41

    """
    model = am.Brook(reionization = True)
    norm = get_sum_fields(model,10**5,star_masses)
    print norm, 'norm for 10^5'
    norm = get_sum_fields(model,10**5,star_masses)
    print norm, 'norm for 10^5'
    norm = get_sum_fields(model,10**4,star_masses)
    print norm, 'norm for 10^4'
    norm = get_sum_fields(model,10**4,star_masses)
    print norm, 'norm for 10^4'

    norm = get_sum_fields(model,10**3,star_masses)
    print norm, 'norm for 10^3'
    norm = get_sum_fields(model,10**3,star_masses)
    print norm, 'norm for 10^3'

    return norm
    """

def get_baseline_reioncat(lum=4): # not updated to mstar_to_mhalo
    if lum == 4:
        return 11.30
    if lum == 3:
        return 15.90
    if lum == 5:
        return 5.46


    """
    model=am.Moster(reionization=True,catreion=True,  z=13, vmax_ach=vmax_ach,vmax_filt=vmax_filt)
    z13 = get_sum_fields(model,10**lum,star_masses)

    norm = get_sum_fields(model,10**lum,star_masses)
    print norm, 'norm for 10^4'
    norm = get_sum_fields(model,10**lum,star_masses)
    print norm, 'norm for 10^4'
    norm = get_sum_fields(model,10**3,star_masses)
    print norm, 'norm for 10^3'
    norm = get_sum_fields(model,10**3,star_masses)
    print norm, 'norm for 10^3'

    model=am.Moster(reionization=True,catreion=True,  z=13, vmax_ach=vmax_ach,vmax_filt=vmax_filt)
    norm = get_sum_fields(model,10**5,star_masses)
    print norm, 'norm for 10^5'
    norm = get_sum_fields(model,10**5,star_masses)
    print norm, 'norm for 10^5'


    return norm
    """




# 1) reionization model: at 4 different redshifts
def change_reion_z(lum=4):
    z13 = get_baseline_reioncat(lum)

    model=am.Moster(reionization=True,catreion=True,z=14, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
    z14 = get_sum_fields(model,10**lum,star_masses)
    model=am.Moster(reionization=True,catreion=True,z=11, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
    z11 = get_sum_fields(model,10**lum,star_masses)
    model=am.Moster(reionization=True,catreion=True,z=9, vmax_ach=vmax_ach, vmax_filt=vmax_filt)
    z9 = get_sum_fields(model,10**lum,star_masses)

    del14 = z14/z13 -1
    del11 = z11/z13 - 1
    del9 = z9/z13 - 1
    print del14*100, del11*100, del9*100
    return del14*100, del11*100, del9*100
# 10.069396252602369, 18.065579458709213, 31.179736294240112
"""
10.5652 z=13
11.4308 z=12
11.9286 z=11
12.9832 z=9
2.4768   z=13 with 15km/s and 50km/s cut-offs instead of 10 and 25
"""


# 1) reionization model: at 4 different redshifts
def change_reion_vmax(lum=4):
    z13 = get_baseline_reioncat(lum)

    #model=am.Moster(reionization=True,catreion=True,z=13, vmax_ach=15, vmax_filt=30)
    model=am.Moster(reionization=True,catreion=True,z=13, vmax_ach=vmax_ach*1.25, vmax_filt=vmax_filt*1.25)
    vup = get_sum_fields(model,10**lum,star_masses)
    #model=am.Moster(reionization=True,catreion=True,z=13, vmax_ach=5, vmax_filt=20)
    model=am.Moster(reionization=True,catreion=True,z=13, vmax_ach=vmax_ach*.75, vmax_filt=vmax_filt*.75)
    vdown = get_sum_fields(model,10**lum,star_masses)

    delup = vup/z13 -1
    deldown = vdown/z13 - 1
    print delup*100, deldown*100 
    return delup*100, deldown*100


def change_reion_vmax_brook(lum=4):
    z13 = get_baseline_Brook(lum)

    #model=am.Moster(reionization=True,catreion=True,z=13, vmax_ach=15, vmax_filt=30)
    model=am.Brook(reionization=True,catreion=True,z=13, vmax_ach=vmax_ach*1.25, vmax_filt=vmax_filt*1.25)
    vup = get_sum_fields(model,10**lum,star_masses)
    #model=am.Moster(reionization=True,catreion=True,z=13, vmax_ach=5, vmax_filt=20)
    model=am.Brook(reionization=True,catreion=True,z=13, vmax_ach=vmax_ach*.75, vmax_filt=vmax_filt*.75)
    vdown = get_sum_fields(model,10**lum,star_masses)

    delup = vup/z13 -1
    deldown = vdown/z13 - 1
    print delup*100, deldown*100 
    return delup*100, deldown*100




#2) Moster - shift infall distribution by 0.23 Gyrs
# tinfall = time since z=0. how long ago
def change_infall(lum=4):
    norm = get_baseline_moster(lum=lum)   # get_sum_fields(model,10**4,star_masses)

    model=am.Moster(reionization=True, t_offset = 0.23)
    earlier = get_sum_fields(model,10**lum,star_masses) # closer to bigbang.

    model=am.Moster(reionization=True, t_offset = -0.23)
    later = get_sum_fields(model,10**lum,star_masses) # closer to z=0

    #model=am.Moster(reionization=True, t_offset = 1.0)
    #earlier1 = get_sum_fields(model,10**4,star_masses)

    #model=am.Moster(reionization=True, t_offset = -1)
    #later1 = get_sum_fields(model,10**4,star_masses)

    del_earlier= earlier/norm -1
    del_later= later/norm -1
    #del_earlier1= earlier1/norm -1
    #del_later1= later1/norm -1
    print del_earlier*100, del_later*100, 'earlier (larger t), later (smaller t)'
    return del_earlier*100, del_later*100  #, del_earlier1*100, del_later1*100
# values for +.5, -.5, +1, -1    
# 1.9029528578829114, -1.804524261785545, 3.9060611293386227, -3.6746675876359935



#3) GK16 - increase the sigma by a lot    
def change_scatter(lum=4):
    norm = get_baseline_GK16(lum)
    #model = am.GarrisonKimmel16(reionization=True, sigma=0.4)
    #sigdown = get_sum_fields(model,10**4,star_masses)
    
    #model = am.GarrisonKimmel16(reionization=True, sigma=1.2)
    #sigup = get_sum_fields(model,10**4,star_masses)
    
    model = am.GK16_grow(reionization=True, plotgamma=0.0)
    gk162 = get_sum_fields(model,10**lum,star_masses)
    
    model = am.GK16_grow(reionization=True, plotgamma=-.5)
    gk165 = get_sum_fields(model,10**lum,star_masses)

    #del_sigdown = sigdown/norm - 1
    #del_sigup = sigup/norm - 1
    del_gk162 = gk162/norm - 1
    del_gk165 = gk165/norm - 1
    #return del_sigdown*100, del_sigup*100, del_gk162*100, del_gk165*100
    print del_gk162*100, del_gk165*100, 'gamma =0 and gamma = -0.5'
    return del_gk162*100, del_gk165*100
    # 8.3766425593837965, -8.2617591387025726, 5.0373644649167826, -18.566254909898571



# 4) change M* of field halos by +/- 10%
def change_mstar(lum=4):
    base = get_baseline_moster(lum) 
    model=am.Moster(reionization=True)

    down = get_sum_fields(model,10**lum,star_masses*.75)
    up = get_sum_fields(model,10**lum,star_masses*1.25)
    delup = up/base - 1.
    deldown = down/base - 1
    print delup*100, deldown*100
    return delup*100, deldown*100
# delup: 10.312, 10.385 = 10 deldown: -11.93, -11.99 = -12




## std error = 0.02 for alpha = -1.8054104588, K = 0.000673927822796


## first 16
# alpha = -1.91225193546
# K = 0.00481748488821
# stderr = 0.0888785222922


## second 17
# alpha = -1.77700887695
# K = 0.000398018832465
# 0.0287554270974 std err




# 5) use the uncertainty on the fit on alpha and K and get the extremes
# need to get new values for the SHMF
def shmf_uncertainty(lum=4):
    norm = get_baseline_moster(lum)   # get_sum_fields(model,10**4,star_masses)

# slope = 1.88, K =  0.00246631667892
# slope = 1.75, K =  0.000193663974252

    # first16
    model=am.Moster(reionization=True,hostSHMF=True, hostAlpha= 1.88 , hostK= 0.00246631667892) 
    first16 = get_sum_fields(model,10**lum,star_masses)

    # second17
    model=am.Moster(reionization=True,hostSHMF=True, hostAlpha= 1.75, hostK= 0.000193663974252)
    second17 = get_sum_fields(model,10**lum,star_masses)

    del_first16= first16/norm -1
    del_second17= second17/norm -1
    print del_first16*100, del_second17*100
    return del_first16*100, del_second17*100
    # (-25.844947584753818, 10.982159819007631)


#7) Change from Field SHMF to Host SHMF  (must use rmaxcut=False)
def field_v_host():
    model=am.Moster(reionization=True) # self.alpha = 1.812, self.K = .000761
    field = get_sum_fields(model,10**4,star_masses)
    model=am.Moster(reionization=True,hostSHMF=True)
    host = get_sum_fields(model,10**4,star_masses)

    # 7.0 as lower limit as all 33 halos
    model=am.Moster(reionization=True,hostSHMF=True, hostAlpha=1.811243 , hostK= 0.000753) #1.811243 & 0.000753
    field2 = get_sum_fields(model,10**4,star_masses)

    # 7.5 as lower limit   1.805410 & 0.000674
    model=am.Moster(reionization=True,hostSHMF=True, hostAlpha=1.805410, hostK= 0.000674) 
    field3 = get_sum_fields(model,10**4,star_masses)

    del_host= host/field -1
    del_field2= field2/field -1
    del_field3= field3/field -1
    return del_host*100, del_field2*100, del_field3*100
    # -18.3101763625654, 0.97363223727522374, 2.2608364228543731
    
    
    
#9) host halo mass +/- 25%
def change_mhost(lum=4):
    base = get_baseline_moster(lum)   # get_sum_fields(model,10**4,star_masses)
    model=am.Moster(reionization=True)

    dm_masses =  model.stellar_to_halo_mass(star_masses,a=1.0)
    star_masses_up = model.getStellarMass(dm_masses*1.22,a=1)
    star_masses_down = model.getStellarMass(dm_masses*0.82,a=1)

    down = get_sum_fields(model,10**lum,star_masses_down)
    up = get_sum_fields(model,10**lum,star_masses_up)
    delup = up/base - 1.
    deldown = down/base - 1
    print delup*100, deldown*100
    return delup*100, deldown*100
  # (11.196523485517407, -10.275454907460302)  # need higher N_iter. shouldn't it be 10% exactly?




#get_baseline_moster()
#get_baseline_GK16()


# 1) reionization model: at 4 different redshifts
#change_reion_z(lum=4)
# z = 14: - 5.4
# z = 11: + 7.6
# z = 9: + 20.2


#change_reion_z(lum=3)
#change_reion_z(lum=3)
#change_reion_z(lum=3)
# z = 14: - 12.5   (-12)
# z = 11: 25.8    (26)
# z = 9:  64.9   (65)

#change_reion_z(lum=5)
#change_reion_z(lum=5)
#change_reion_z(lum=5)
# z = 14: ~ -.3%
# z = 11: +.5 %
# z = 9:  + 1.5 %



# change both vmax values up or down 25%
#change_reion_vmax(lum=4)
# up:  -38
# down: 27
#change_reion_vmax(lum=3)
# up: -47
# down: +70

#change_reion_vmax(lum=5)
#change_reion_vmax(lum=5)
#change_reion_vmax(lum=5)
# up: - 13 %
# down: + 3 %



#2) Moster - shift infall distribution by 0.23 Gyrs
#change_infall()
# tinfall up .23 (closer to bigbang): .943 for 25000, .8973. conclusion: 0.9
# tinfall down .23 (closer to z=0): -.687 for 25000, -.673. conclusion: -0.7
#change_infall(lum=3)
# up: +.7:   about the same < 1%
# down: -.7:  < 1%



#3) GK16 - increase the sigma by a lot
#change_scatter()
# gamma = 0:  8.80759971371, 9.26722623463, 9.00357863231, 9.11679354545. Conclusion: 9
# gamma = -0.5: -22.4347712929, -22.5367948468, -22.6463660616, -22.8639469061.  Conclusion: -23

# with new barber
# gamma = 0: + 16.    For 10**3: -2
# gamma = -0.5: -18.  For 10**3: -9

# with new barber and with mstar_to_mhalo
#  gamma = 0: +24.47        For 10**3:  + 4.8
#  gamma = -0.5: -32        For 10**3:   -24


# additional 8%,    and 7%
# additional -14%,   and -15%



# 4) change M* of field halos by +/- 25%
# change_mstar()
#change_mstar(lum=3)
# From Roediger:  random errors of 23-30%  .09 - .13 dex. .1 dex is the low end. could be up to 0.3 or 0.6 dex for 10^8-10^11 msun
# e^.25 = 1.28, e^.4 = 1.49
# maybe I should go with 25-30%??
# Error not uncommon of 25% in Mstar
# 10 and -12 for both lum=3 and 4



# 5) use the uncertainty on the fit on alpha and K and get the extremes
#shmf_uncertainty()
#shmf_uncertainty(lum=4)
# slope = 1.88, K =  0.00246631667892.   -12.115978456, -12.5669658887, -12.4879712747. average is  -12
# slope = 1.75, K =  0.000193663974252.   10.0631956912,  9.83303411131, 9.96373429084. average is   10


# lum=4:  -12 and 10
# lum=3:  -9 and 8


#7) Change from Field SHMF to Host SHMF  (must use rmaxcut=False)
#field_v_host()

    
#9) host halo mass +/- 25%
#change_mhost()
# change_mhost(lum=3)
# +/- 26 % is the answer for both lum=3 and 4


"""
#### infall times host halo
host_jackknife = np.array([7.65699519351 ,7.64654593114 ,7.66179279557 ,7.67879718798 ,7.66082868531 ,7.65314905407 ,7.59670797158 ,7.64976459894 ,7.61520337655 ,7.58810134737 ,7.62625843549 ,7.63277126785 ,7.62451654038 ,7.62943857486 ,7.60847144661 ,7.70246419899 ,7.62726831506 ,7.63775148588 ,7.63515592862 ,7.64112107136 ,7.62391950938 ,7.6398573557 ,7.65368177523,7.6751677347 ,7.65961250838,7.64259220479,7.58926407723,7.61825527621,7.65134366953,7.59558271396,7.6391682698 ,7.6734597656 ,7.63869603811])

mean = 7.63869603811
N = len(host_jackknife)
sig2 = (N-1) * np.sum((host_jackknife - mean)**2) / N

# sig2 = 0.021 this is for host halos



#### infall times field halo
field_jackknife = np.array([7.45075566082,7.43644909453,7.4339121995,7.44960669155,7.44128527192,7.43586377931,7.43347250956,7.62177028084,7.44557222835,7.4479588702,7.44593685404,7.47307580853,7.38480617675,7.44448698234,7.44694501506,7.39853429166,7.45705235559,7.4446268784,7.42250326437,7.4609546539,7.44394706416,7.44494287612,7.44828756174,7.44599329834,7.43861544726,7.42092190491,7.41141755149,7.44872444074,7.43939498695,7.56899967318,7.42886069686,7.43838931246,7.44828756174])

mean = 7.44828756174
N = len(field_jackknife)
sig2 = (N-1) * np.sum((field_jackknife - mean)**2) / N

# sig2 = 0.054 this is for field halos
# sig = sqrt(.054) = 0.234


sampling bias:
mm = mean - (N-1)*(np.mean(field_jackknife) - mean)
print mm
print np.mean(field_jackknife) - mean


"""





"""
SHMF jackknife scheme
get list of alpha and k


alphas = np.array([-1.8188900910003025, -1.8238898420685112, -1.8171319545955282, -1.816073787537215, -1.8243071405090254, -1.8206196910538617, -1.8141440129138615, -1.8242862572905008, -1.8149209200804637, -1.8195231846314297, -1.8219590805762769, -1.8238328037024143, -1.8208614336796412, -1.8227162335688878, -1.8212733298721282, -1.8127832541046363, -1.8233171357020499, -1.8207073135395699, -1.8435490084503228, -1.8261974589998453, -1.8199921230271705, -1.812869674221977, -1.820475228417872, -1.8217583927084549, -1.8151435674913168, -1.8239571057137931, -1.8302446662114957, -1.8228915784639301, -1.8207549240806387, -1.8185756967342974, -1.8106611121809761, -1.8106611121809761])

Klist = np.array([0.0007324191372882771, 0.00081205205612041491, 0.00071022917161181703, 0.00070468489832668096, 0.00081486408543557647, 0.00076162915081789658, 0.00067394158174228038, 0.00081204541675318758, 0.00067436339713238837, 0.00073973610903751226, 0.00078457357511230511, 0.0008143596082301962, 0.00076657382668522569, 0.00079466175516889736, 0.00076704390536096049, 0.00066067877969638617, 0.00079653026453200442, 0.00076180259714869226, 0.0011635764316249343, 0.00083894146092186334, 0.00074884349245938275, 0.00065592584538909425, 0.00075866201755552887, 0.00077767346192275796, 0.00068917083930255011, 0.00081422120739210377, 0.00091089908382822589, 0.00079835850200502218, 0.00076225039309826034, 0.00072863057776608017, 0.00063484586246424432, 0.00063484586246424432])



mean_alpha = -1.8106611121809761
mean_k = 0.00063484586246424432
N = len(alphas)
sig2_alpha = (N-1) * np.sum((alphas - mean_alpha)**2) / N

sig2_k = (N-1) * np.sum((Klist - mean_k)**2) / N

Results:
sig2_alpha =  0.0042153878872658226
sig_alpha = 0.064926018569336455

sig2_k = 8.1072022284576664e-07
sig_k = 0.00090040003489880356

# plot alpha vs k:
plt.scatter(alphas, Klist, color='black')
plt.scatter(mean_alpha, mean_k, color='red')

# fit a line to the data. take my mean alpha and add .065 and subtract .065
m,b = np.polyfit(np.append(mean_alpha,alphas), np.append(mean_k, Klist), deg=1)

x_array = np.arange(np.min(alphas), np.max(alphas), .002)
y_array = m*x_array+b
plt.plot(x_array, y_array, lw=1)

plt.xlabel('alpha')
plt.ylabel('k')
plt.ylim((np.min(Klist),np.max(Klist) ))
plt.xlim((np.min(alphas), np.max(alphas)))
plt.savefig('k_v_alpha')
plt.close()




# upper and lower alpha limits:
upper: sig_alpha + (- mean_alpha) = 1.8755871307503127.   or 1.88.  k = 0.0016741701144383364       0.00246631667892 Kfix
lower (-mean_alpha) - sig_alpha = 1.7457350936116396      or 1.75.  k = -0.00031265764620198866     0.000193663974252







"""
