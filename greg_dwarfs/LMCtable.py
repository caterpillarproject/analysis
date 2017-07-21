import numpy as np
import abundance_matching as am
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PlotParams import *
import RadialDependence as rd
import DwarfMethods as dm



def get_ngreater(mstar,min_lum,model):
        dm_masses = model.mstar_to_mhalo(mstar,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        return model.ngreater(dm_masses,min_lum, percentile=True)

def get_P_at_least_one(mstar,min_lum,model):
        dm_masses = model.mstar_to_mhalo(mstar,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        return model.P_at_least_one(dm_masses,min_lum)

def get_rvir(StellarMasses,model):
        dm_masses = model.mstar_to_mhalo(StellarMasses,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        #dummy, rvir = dm.convert_Mhalo_d_d(dm_masses, model.delta_crit, 103.86)
        return model.mass_to_rvir(dm_masses)  # this does properly convert from one radius definition to the virial definition


# MUST BE TESTED. THIS FAILED
def get_in_fov(radii, dists, camera_radius = 0.56):
        Z = 1.5
        fov = camera_radius * np.pi/180 # convert to radians
        fov_radius = fov * dists # convert to kpc at dist of target
        R = fov_radius/radii
        frac = rd.get_normed_los_better(R,Z) # this has the proper reionization driven radial function
        return frac


def generate_latex_table(Names, stellar_masses, distances,  reion=True, min_lum = 10**5):
    with open('/nfs/blank/h4231/gdooley/Dropbox/DwarfsOfDwarfs/LMCPaper/TableLMC44.tex' ,'w') as f:
        f.write("\\begin{table*}\n")
        f.write("\\tablewidth{0.97\\textwidth}\n")
        f.write("\\centering\n")
        f.write("\\caption{Mean number of observable satellites with $M_* > 10^5$\msun around isolated field dwarf galaxies}\n")
        f.write("\\label{table:results}\n")
        #f.write(" \\tiny \n ")
        f.write("\\begin{tabular}{lccccccccccccc} \n")
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\textbf{Name} &  \\boldmath$M_* $ &  \\boldmath$D_{\\sun} $ &  \\multicolumn{5}{c}{\\textbf{Brook}} && \\multicolumn{5}{c}{\\textbf{GK16}}  \\\\\n")
        f.write(" \\cline{4-8}\\cline{10-14}  \\\\ \n ")

        ## kpc is on the same line as Rvir
        f.write(" & $[10^8 \\, \\mathrm{M_\odot}]$  &  $[\\mathrm{Mpc}]$ & $\\bar{N}_{\\rm{lum}}$ & 20/80\\% & $\\bar{N}_{\\rm{fov}}^{1.5}$ & $\\bar{N}_{\\rm{fov}}^{2.2}$ &   $R_{\\rm{vir}}$ &   & $\\bar{N}_{\\rm{lum}}$ & 20/80\\% & $\\bar{N}_{\\rm{fov}}^{1.5}$ & $\\bar{N}_{\\rm{fov}}^{2.2} $ &   $R_{\\rm{vir}}$ \\\\\n")
        #   & & $\\bar{N}_{\\rm{lum}}$  & 20/80\\% & $R_{\\rm{vir}}$ &  & $\\bar{N}_{\\rm{lum}}$ & 20/80\\% & $R_{\\rm{vir}}$  \\\\\n" ) 

        print 'def on the newer version.'
        Ngreater = []; Pgt1=[]; Rvir=[]; Nfov1=[]; Nfov2=[]
        #for model in [am.Moster(reionization=reion), am.GarrisonKimmel(reionization=reion), am.GK16_grow(reionization=reion),am.Brook(reionization=reion)]: 
        for model in [am.Brook(lmcSHMF=True, reionization=reion), am.GK16_grow(lmcSHMF=True,reionization=reion)]:
            print model.label, 'on this model'
            Nlum = get_ngreater(stellar_masses,min_lum,model)
            print 'finished ngreater'
            Ngreater.append(Nlum)
            #Ngreater.append(get_ngreater(stellar_masses, min_lum, model))
            Pgt1.append(get_P_at_least_one(stellar_masses,min_lum,model))
            radii = get_rvir(stellar_masses,model)
            Rvir.append(radii)
            Nfov1.append( Nlum[0] * get_in_fov(radii, distances*1000, camera_radius = 1.5/2))
            Nfov2.append( Nlum[0] * get_in_fov(radii, distances*1000, camera_radius = 2.2/2))

        Nsubs = np.array(Ngreater)[:,0]
        Nstd = np.array(Ngreater)[:,1]
        N10 = np.array(Ngreater)[:,2]
        N90 = np.array(Ngreater)[:,3]

        for name, mstar, nsubs, pgt1, nstd, n10, n90, rvir, dist, fov1, fov2 in zip(Names,stellar_masses, Nsubs.T, np.array(Pgt1).T, Nstd.T, N10.T, N90.T, np.array(Rvir).T, distances, np.array(Nfov1).T, np.array(Nfov2).T):
            f.write("\\hline\n")
            if name[0:3]=='SMC' or name[0:3]=='LMC':
                print 'in LMC SMC'
                f.write("%s & %.3g & %s     & $%.2f$ & $%.2g / %.2g$ & %s & %s & %d &    & $%.2f$ & $%.2g/ %.2g$ & %s & %s & %d  \\\\\n" %(name,mstar*1e-8, '\\nodata', nsubs[0],n10[0],n90[0], '\\nodata', '\\nodata',  rvir[0],    nsubs[1],n10[1],n90[1], '\\nodata', '\\nodata', rvir[1]))
            else:
                f.write("%s & %.3g & %.3g     & $%.2f$ & $%.2g / %.2g$ & $%.2f$ & $%.2f$ & %d &    & $%.2f$ & $%.2g/ %.2g$ & $%.2f$ & $%.2f$ & %d  \\\\\n" %(name,mstar*1e-8, dist, nsubs[0],n10[0],n90[0], fov1[0], fov2[0],  rvir[0],    nsubs[1],n10[1],n90[1],fov1[1],fov2[1], rvir[1]))
            #f.write("%s & %.3g      & $%.2f$ & $%.2g / %.2g$ & %.2f &    & $%.2f$ & $%.2g/ %.2g$ & %.2f &    & $%.2f$ & $%.2g/%.2g$ & %.2f  &     & $%.2f$ & $%.2g/%.2g$ & %.2f    \\\\\n" %(name,mstar*1e-6, nsubs[0],n10[0],n90[0],rvir[0],    nsubs[1],n10[1],n90[1],rvir[1],    nsubs[2],n10[2],n90[2],rvir[2],  nsubs[3],n10[3],n90[3],rvir[3]))

        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\vspace{-5mm}\n")
        f.write("\\tablecomments{Mean number of satellites with $M^* > 10^5$\msun expected to exist within the virial volume of known isolated field galaxies as predicted with the Brook and GK16 models. The $20^{\\rm{th}}$ and $80^{\\rm{th}}$ percentile of the satellite abundance distributions are listed, indicating a $>80\%$ chance of at least one satellite in even the smallest galaxies considered. $\\bar{N}_{\\rm{fov}}^{1.5}$ and $\\bar{N}_{\\rm{fov}}^{2.2}$ indicate the mean number of satellites within a field of view of diameter $1.5^{\\circ}$ and $2.2^{\\circ}$ respectively. These are the apertures of Hyper SuprimeCam and DECam. Also shown is the inferred virial radius of each galaxy.} \n")
        f.write("\\end{table*}\n")

#*The LMC listed assumes an isolated galaxy with the same stellar mass as the LMC.

data = np.loadtxt('lmc_analogs.txt', dtype='str', usecols = (0,6,8))
#data=data[0:-1]  # would use to exclude the SMC
DwarfNames = data[1:,0]
Distances =  np.array(data[1:,1], dtype='float')  # in Mpc
StellarMasses = np.array(data[1:,2], dtype='float') # in solar masses

argsorted = np.argsort(StellarMasses)
DwarfNames = DwarfNames[argsorted]
Distances = Distances[argsorted]
StellarMasses = StellarMasses[argsorted]


#generate_latex_table(DwarfNames, StellarMasses, Distances,  reion=True, min_lum = 10**5)
generate_latex_table(DwarfNames, StellarMasses, Distances,  reion=True, min_lum = 10**4)
