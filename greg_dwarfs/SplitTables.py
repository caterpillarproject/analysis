import numpy as np
import abundance_matching as am
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PlotParams import *

min_lum = 10**4
DwarfNames = ['Leo T', 'KKR 3', 'Phoenix', 'KKH 86','Antlia','KKR 25','Aquarius','DDO 113','ESO 294- G 010','Sagittarius dIrr',
'ESO 410- G 005','KKH 98','Leo A','GR 8','Pegasus dIrr','UGC 9128','UGC 4879','UGCA 86','DDO 99','UKS 2323-326',
'UGC 8508','NGC 4163','Sextans A','DDO 125','DDO 190','Sextans B','IC 3104','NGC 3109','NGC 6822','IC 1613',
'IC 4662','IC 5152']

StellarMasses = 10**6*np.array([0.14,0.54,0.77,0.82,1.3,1.4,1.6,2.1,2.7,3.5,3.5,4.5,6.0,6.4,6.6,7.8,8.3,16,16,17,19,37,44,47,51,52,62,76,100,100,190,270])

def get_ngreater(mstar,min_lum,model):
        dm_masses = model.mstar_to_mhalo(mstar,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        return model.ngreater(dm_masses,min_lum)
        
def get_P_at_least_one(mstar,min_lum,model):
        dm_masses = model.mstar_to_mhalo(mstar,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        return model.P_at_least_one(dm_masses,min_lum)

def get_rvir(StellarMasses,model):
        dm_masses = model.mstar_to_mhalo(StellarMasses,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        return model.mass_to_rvir(dm_masses)



"""
%run DwarfIrregularTable.py
StellarMasses = 10**6*np.array([0.14,0.54,0.77,0.82,1.3,1.4,1.6,2.1,2.7,3.5,3.5,4.5,6.0,6.4,6.6,7.8,8.3,16,16,17,19,37,44,47,51,52,62,76,100,100,190,270])
mstar = StellarMasses[-5:]
min_lum = 10**5
model = am.Brook(reionization=True)
get_P_at_least_one(mstar,min_lum,model)

"""


def generate_latex_table1(Names,Mstar,reion=True):
    with open('/nfs/blank/h4231/gdooley/Dropbox/DwarfsOfDwarfs/Paper/Table1half.tex' ,'w') as f:
        f.write("\\begin{table*}\n")
        f.write("\\tablewidth{0.97\\textwidth}\n")
        f.write("\\centering\n")
        f.write("\\caption{Local Group dIrrs}\n")
        f.write("\\label{table:results1}\n")
        f.write("\\begin{tabular}{lcccccccc} \n")
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\textbf{Name} &  \\boldmath$M_* \\, [10^6 \\, \\mathrm{M_\odot}]$ & \\multicolumn{3}{c}{\\textbf{Moster}}  & & \\multicolumn{3}{c}{\\textbf{GK14}}  \\\\\n" )
        f.write(" \\cline{3-5} \\cline{7-9} \\\\ \n " )
        f.write(" &  & $\\bar{N}_{\\rm{lum}}$ &  $P( \\geq 1)$ & $R_{\\rm{vir}} \, \\rm{[kpc]}$ & &  $\\bar{N}_{\\rm{lum}}$ & $P( \\geq 1)$ & $R_{\\rm{vir}} \, \\rm{[kpc]}$  \\\\\n")
        
        print 'def on the newer version'
        Ngreater = []; Pgt1=[]; Rvir=[]
        for model in [am.Moster(reionization=reion), am.GarrisonKimmel(reionization=reion)]:
            print model.label, 'on this model'
            Ngreater.append(get_ngreater(StellarMasses,min_lum,model))
            Pgt1.append(get_P_at_least_one(StellarMasses,min_lum,model))
            Rvir.append(np.round(get_rvir(StellarMasses,model)))

        Nsubs = np.array(Ngreater)[:,0]
        Nstd = np.array(Ngreater)[:,1]

        for name, mstar, nsubs, pgt1, nstd, rvir in zip(Names,Mstar, Nsubs.T, np.array(Pgt1).T, Nstd.T, np.array(Rvir).T):
            f.write("\\hline\n")
            f.write("%s & %.2f      & %.2f $\pm$ %.2f & %.2f & %d &     & %.2f $\pm$ %.2f & %.2f & %d  \\\\\n" %(name,mstar*1e-6, nsubs[0],nstd[0],pgt1[0],rvir[0],    nsubs[1],nstd[1],pgt1[1],rvir[1]))
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\tablecomments{Mean number of satellites with $M^* > 10^4$\msun expected to exist within the virial volume of known local group Dwarf Irregular and Dwarf Spheroidal galaxies for the Moster and GK14 AM models. One standard deviation values given as a $\pm$. Also shown is the probability of finding at least one satellite around each galaxy, and the inferred virial radius of each galaxy.} \n")
        f.write("\\end{table*}\n")




def generate_latex_table2(Names,Mstar,reion=True):
    with open('/nfs/blank/h4231/gdooley/Dropbox/DwarfsOfDwarfs/Paper/Table2half.tex' ,'w') as f:
        f.write("\\begin{table*}\n")
        f.write("\\tablewidth{0.97\\textwidth}\n")
        f.write("\\centering\n")
        f.write("\\caption{Local Group dIrrs}\n")
        f.write("\\label{table:results2}\n")
        f.write("\\begin{tabular}{lcccccccc} \n")
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\textbf{Name} &  \\boldmath$M_* \\, [10^6 \\, \\mathrm{M_\odot}]$ & \\multicolumn{3}{c}{\\textbf{GK16}} &  & \\multicolumn{3}{c}{\\textbf{Brook}}  \\\\\n" )
        f.write(" \\cline{3-5} \\cline{7-9} \\\\ \n " )

        f.write(" &  & $\\bar{N}_{\\rm{lum}}$ &  $P( \\geq 1)$ & $R_{\\rm{vir}} \, \\rm{[kpc]}$ & & $\\bar{N}_{\\rm{lum}}$ & $P( \\geq 1)$ & $R_{\\rm{vir}} \, \\rm{[kpc]}$  \\\\\n")
        
        print 'def on the newer version'
        Ngreater = []; Pgt1=[]; Rvir=[]
        for model in [ am.GK16_grow(reionization=reion),am.Brook(reionization=reion) ]:
            print model.label, 'on this model'
            Ngreater.append(get_ngreater(StellarMasses,min_lum,model))
            Pgt1.append(get_P_at_least_one(StellarMasses,min_lum,model))
            Rvir.append(np.round(get_rvir(StellarMasses,model)))

        Nsubs = np.array(Ngreater)[:,0]
        Nstd = np.array(Ngreater)[:,1]

        for name, mstar, nsubs, pgt1, nstd, rvir in zip(Names,Mstar, Nsubs.T, np.array(Pgt1).T, Nstd.T, np.array(Rvir).T):
            f.write("\\hline\n")
            f.write("%s & %.2f      & %.2f $\pm$ %.2f & %.2f & %d &     & %.2f $\pm$ %.2f & %.2f & %d  \\\\\n" %(name,mstar*1e-6, nsubs[0],nstd[0],pgt1[0],rvir[0],    nsubs[1],nstd[1],pgt1[1],rvir[1]))
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\tablecomments{Same as Table~\\ref{table:results1} except for the GK16 and Brook AM models.}\n")
        f.write("\\end{table*}\n")



# TABLE 2
#generate_latex_table1(DwarfNames, StellarMasses)
generate_latex_table2(DwarfNames, StellarMasses)

