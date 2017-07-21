import numpy as np
import abundance_matching as am
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PlotParams import *
import RadialDependence as rd


#And  XVIII  - distance 1355,   stellar mass 0.63
#And XXVIII -  distance 661    stellar mass 0.21

min_lum = 10**4

#Sag DIG  for Sagittarius dIrr
DwarfNames = ['Leo T', 'And XXVIII', 'KKR 3', 'Tucana', 'Leo P', 'And XVIII', 'Phoenix', 'KKH 86','Antlia','KKR 25','Aquarius','DDO 113','Cetus','ESO 294-G010','Sagittarius dIrr','ESO 410-G005','KKH 98','Leo A','GR 8','Pegasus dIrr','UGC 9128','UGC 4879', 'KK 258',  'UGCA 86','DDO 99','UKS 2323-326','UGC 8508', 'KKs 3', 'NGC 4163','WLM','Sextans A','DDO 125','DDO 190','Sextans B','IC 3104','NGC 3109','NGC 6822','IC 1613','IC 4662','IC 5152']

StellarMasses = 10**6*np.array([0.14,0.21,  0.54,0.56, 0.56, 0.63, 0.77,0.82,1.3,1.4,1.6,2.1,2.6,2.7,3.5,3.5,4.5,6.0,6.4,6.6,7.8,8.3,14, 16,16,17,19,23,37,43,44,47,51,52,62,76,100,100,190,270])

Distances = np.array([ 417, 661, 2188, 887, 1620, 1355, 415, 2582, 1349, 1905, 1072, 2951, 755, 2032, 1067, 1923, 2523, 798, 2178, 920, 2291, 1361, 2230, 2965, 2594, 2208, 2582, 2120,  2858, 933, 1432, 2582, 2793, 1426, 2270, 1300, 459, 755, 2443, 1950 ])  # in kpc

N_iter = 20000  # make 20000 for the real thing. 5000 min acceptable

def get_ngreater(mstar,min_lum,model):
        dm_masses = model.mstar_to_mhalo(mstar,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        return model.ngreater(dm_masses,min_lum, percentile=True)

def get_P_at_least_one(mstar,min_lum,model):
        dm_masses = model.mstar_to_mhalo(mstar,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        return model.P_at_least_one(dm_masses,min_lum)

def get_rvir(StellarMasses,model):
        dm_masses = model.mstar_to_mhalo(StellarMasses,a=1)  # stellar_to_halo_mass(mstar,a=1.0)
        #dummy, rvir = dm.convert_Mhalo_d_d(dm_masses, model.delta_crit, 103.86)
        #return rvir
        return model.mass_to_rvir(dm_masses)  # this does the correct conversion

# just want the multiplicative fraction
def get_in_fov(radii, dists, camera_radius = 0.56):
        Z = 1.5
        fov = camera_radius * np.pi/180 # convert to radians
        fov_radius = fov * dists # convert to kpc at dist of target
        R = fov_radius/radii
        frac = rd.get_normed_los_better(R,Z) # rd.get_normed_los_protected(R,Z) 
        return frac
        

# \multicolumn{2}{r}{} $\%$

"""
%run DwarfIrregularTable.py
StellarMasses = 10**6*np.array([0.14,0.54,0.77,0.82,1.3,1.4,1.6,2.1,2.7,3.5,3.5,4.5,6.0,6.4,6.6,7.8,8.3,16,16,17,19,37,44,47,51,52,62,76,100,100,190,270])
mstar = StellarMasses[-5:]
min_lum = 10**5
model = am.Brook(reionization=True)
get_P_at_least_one(mstar,min_lum,model)

"""


def generate_latex_table(Names,reion=True):
    with open('/nfs/blank/h4231/gdooley/Dropbox/DwarfsOfDwarfs/Paper/Table1.tex' ,'w') as f:
        f.write("\\begin{table}\n")
        f.write("\\tablewidth{0.97\\textwidth}\n")
        f.write("\\centering\n")
        f.write("\\caption{Mean number of satellites with $M_* > 10^4$\msun around Local Group dwarf field galaxies}\n")
        f.write("\\label{table:results}\n")
        #f.write(" \\tiny \n ")
        f.write("\\begin{tabular}{lcccccccccccccccc} \n")
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\textbf{Name} &  \\boldmath$M_* $ & \\multicolumn{3}{c}{\\textbf{Moster}} && \\multicolumn{3}{c}{\\textbf{GK14}} && \\multicolumn{3}{c}{\\textbf{GK16}} && \\multicolumn{3}{c}{\\textbf{Brook}} \\\\\n" )  # N Behroozi & N Sawala
        f.write(" \\cline{3-5}\\cline{7-9}\\cline{11-13} \\cline{15-17} \\\\ \n "  )

        ## kpc is on the same line as Rvir
        f.write(" & $[10^6 \\, \\mathrm{M_\odot}]$  & $\\bar{N}_{\\rm{lum}}$ & 20/80\\% &  $P( \\geq 1)$ &   & $\\bar{N}_{\\rm{lum}}$ & 20/80\\% & $P( \\geq 1)$ & & $\\bar{N}_{\\rm{lum}}$  & 20/80\\% & $P( \\geq 1)$ &  & $\\bar{N}_{\\rm{lum}}$ & 20/80\\% & $P( \\geq 1)$  \\\\\n" )  # N Behroozi & N Sawala

        print 'def on the newer version.'
        Ngreater = []; Pgt1=[]; Rvir=[]
        for model in [am.Moster(reionization=reion), am.GarrisonKimmel(reionization=reion), am.GK16_grow(reionization=reion),am.Brook(reionization=reion)]:  #am.Behroozi(reionization=reion), am.Sawala()
            print model.label, 'on this model'
            Ngreater.append(get_ngreater(StellarMasses,min_lum,model))
            Pgt1.append(get_P_at_least_one(StellarMasses,min_lum,model))
            Rvir.append(get_rvir(StellarMasses,model))
            

        Nsubs = np.array(Ngreater)[:,0]
        Nstd = np.array(Ngreater)[:,1]
        N10 = np.array(Ngreater)[:,2]
        N90 = np.array(Ngreater)[:,3]

        for name, mstar, nsubs, pgt1, nstd, n10, n90, rvir in zip(Names,StellarMasses, Nsubs.T, np.array(Pgt1).T, Nstd.T, N10.T, N90.T, np.array(Rvir).T)[0:24]:
            f.write("\\hline\n")

            f.write("%s & %.3g      & $%.2f$ & $%.1g / %.1g$ & %.2f &    & $%.2f$ & $%.1g/ %.1g$ & %.2f &    & $%.2f$ & $%.1g/%.1g$ & %.2f  &     & $%.2f$ & $%.1g/%.1g$ & %.2f    \\\\\n" %(name,mstar*1e-6, nsubs[0],n10[0],n90[0],pgt1[0],    nsubs[1],n10[1],n90[1],pgt1[1],    nsubs[2],n10[2],n90[2],pgt1[2],  nsubs[3],n10[3],n90[3],pgt1[3]))    #, nsubs[4],nstd[4], nsubs[5],nstd[5] ))

        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
        f.write("\\clearpage\n")



        ##### SECOND PAGE OF THE TABLE ######

        f.write("\\begin{table}\n")
        f.write("\\tablewidth{0.97\\textwidth}\n")
        f.write("\\centering\n")
        f.write("\\contcaption{Mean number of satellites with $M_* > 10^4$\msun around Local Group dwarf field galaxies}\n")
        f.write("\\begin{tabular}{lcccccccccccccccc} \n")
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\textbf{Name} &  \\boldmath$M_* $ & \\multicolumn{3}{c}{\\textbf{Moster}} && \\multicolumn{3}{c}{\\textbf{GK14}} && \\multicolumn{3}{c}{\\textbf{GK16}} && \\multicolumn{3}{c}{\\textbf{Brook}} \\\\\n" )  # N Behroozi & N Sawala
        f.write(" \\cline{3-5}\\cline{7-9}\\cline{11-13} \\cline{15-17} \\\\ \n "  )
        ## kpc is on the same line as Rvir
        f.write(" & $[10^6 \\, \\mathrm{M_\odot}]$  & $\\bar{N}_{\\rm{lum}}$ & 20/80\\% &  $P( \\geq 1)$ &   & $\\bar{N}_{\\rm{lum}}$ & 20/80\\% & $P( \\geq 1)$ & & $\\bar{N}_{\\rm{lum}}$  & 20/80\\% & $P( \\geq 1)$ &  & $\\bar{N}_{\\rm{lum}}$ & 20/80\\% & $P( \\geq 1)$  \\\\\n" ) 


        for name, mstar, nsubs, pgt1, nstd, n10, n90, rvir in zip(Names,StellarMasses, Nsubs.T, np.array(Pgt1).T, Nstd.T, N10.T, N90.T, np.array(Rvir).T)[24:]:
            f.write("\\hline\n")
            f.write("%s & %.3g     & $%.2f$ & $%.1g / %.1g$ & %.2f &    & $%.2f$ & $%.1g/ %.1g$ & %.2f &    & $%.2f$ & $%.1g/%.1g$ & %.2f  &     & $%.2f$ & $%.1g/%.1g$ & %.2f    \\\\\n" %(name,mstar*1e-6, nsubs[0],n10[0],n90[0],pgt1[0],    nsubs[1],n10[1],n90[1],pgt1[1],    nsubs[2],n10[2],n90[2],pgt1[2],  nsubs[3],n10[3],n90[3],pgt1[3])) 


        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\tablecomments{Mean number of satellites with $M^* > 10^4$\msun expected to exist within the virial volume of known Local Group dwarf irregular and dwarf spheroidal galaxies as predicted with various AM models. The $20^{\\rm{th}}$ and $80^{\\rm{th}}$ percentile of the satellite abundance distributions are included in the second column. Also shown is the probability of finding at least one satellite around each galaxy, $P( \\geq 1)$.}\n")
        f.write("\\end{table}\n")





def generate_latex_table2(Names, reion=True):
    with open('/nfs/blank/h4231/gdooley/Dropbox/DwarfsOfDwarfs/Paper/Table2.tex' ,'w') as f:
        f.write("\\begin{table}\n")
        f.write("\\tablewidth{0.97\\textwidth}\n")
        f.write("\\centering\n")
        f.write("\\caption{Mean number of satellites with $M_* > 10^4$\msun within a $0.56^\\circ$ radius field of view around Local Group dwarf field galaxies}\n")
        f.write("\\label{table:results2}\n")
        #f.write(" \\tiny \n ")
        f.write("\\begin{tabular}{lcccccccccccccccc} \n")
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\textbf{Name} &  \\boldmath$D_{\\sun} \\, [\\mathrm{kpc}]$ & \\multicolumn{3}{c}{\\textbf{Moster}} && \\multicolumn{3}{c}{\\textbf{GK14}} && \\multicolumn{3}{c}{\\textbf{GK16}} && \\multicolumn{3}{c}{\\textbf{Brook}} \\\\\n" )  # N Behroozi & N Sawala
        f.write(" \\cline{3-5}\\cline{7-9}\\cline{11-13} \\cline{15-17} \\\\ \n "  )

        ## kpc is on the same line as Rvir
        f.write(" &  & $\\bar{N}_{\\rm{lum}}$ &  $\\bar{N}_{\\rm{fov}}$ & $R_{\\rm{vir}} \\, \\rm{[kpc]}$  &  & $\\bar{N}_{\\rm{lum}}$ & $\\bar{N}_{\\rm{fov}}$ & $R_{\\rm{vir}} \\, \\rm{[kpc]}$ & & $\\bar{N}_{\\rm{lum}}$  & $\\bar{N}_{\\rm{fov}}$ & $R_{\\rm{vir}} \\, \\rm{[kpc]}$ &  & $\\bar{N}_{\\rm{lum}}$ & $\\bar{N}_{\\rm{fov}}$ & $R_{\\rm{vir}} \\, \\rm{[kpc]}$   \\\\\n" ) 

        print 'def on the newer version.'
        Ngreater = []; Pgt1=[]; Rvir=[]; Nfov=[]
        for model in [am.Moster(reionization=reion), am.GarrisonKimmel(reionization=reion), am.GK16_grow(reionization=reion),am.Brook(reionization=reion)]:  #am.Behroozi(reionization=reion), am.Sawala()
            print model.label, 'on this model'
            Nlum = get_ngreater(StellarMasses,min_lum,model)
            Ngreater.append(Nlum)
            Pgt1.append(get_P_at_least_one(StellarMasses,min_lum,model))
            radii = get_rvir(StellarMasses,model)
            Rvir.append(radii)
            Nfov.append( Nlum[0] * get_in_fov(radii, Distances))
            
        #print Rvir, 'virial radius of all'
        #print Nfov, 'N fov all'
        Nsubs = np.array(Ngreater)[:,0]
        Nstd = np.array(Ngreater)[:,1]
        N10 = np.array(Ngreater)[:,2]
        N90 = np.array(Ngreater)[:,3]

        for name, dist, nsubs, pgt1, nstd, n10, n90, rvir, fov in zip(Names,Distances, Nsubs.T, np.array(Pgt1).T, Nstd.T, N10.T, N90.T, np.array(Rvir).T, np.array(Nfov).T)[0:24]:
            f.write("\\hline\n")

            f.write("%s & %d    & $%.2f$ & %.2f & %d &    & $%.2f$ & %.2f & %d &    & $%.2f$ & %.2f & %d &     & $%.2f$ & %.2f & %d    \\\\\n" %(name, dist, nsubs[0],fov[0],rvir[0],    nsubs[1],fov[1],rvir[1],    nsubs[2],fov[2],rvir[2],  nsubs[3],fov[3],rvir[3]))    #, nsubs[4],nstd[4], nsubs[5],nstd[5] ))

        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")
        f.write("\\clearpage\n")


        ##### SECOND PAGE OF TABLE #####

        f.write("\\begin{table}\n")
        f.write("\\tablewidth{0.97\\textwidth}\n")
        f.write("\\centering\n")
        f.write("\\contcaption{Mean number of satellites with $M_* > 10^4$\msun within a $0.56^\\circ$ radius field of view around Local Group dwarf field galaxies}\n")
        f.write("\\begin{tabular}{lcccccccccccccccc} \n")
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\textbf{Name} &  \\boldmath$D_{\\sun} \\, [\\mathrm{kpc}]$ & \\multicolumn{3}{c}{\\textbf{Moster}} && \\multicolumn{3}{c}{\\textbf{GK14}} && \\multicolumn{3}{c}{\\textbf{GK16}} && \\multicolumn{3}{c}{\\textbf{Brook}} \\\\\n" )  # N Behroozi & N Sawala
        f.write(" \\cline{3-5}\\cline{7-9}\\cline{11-13} \\cline{15-17} \\\\ \n "  )
        ## kpc is on the same line as Rvir
        f.write(" &  & $\\bar{N}_{\\rm{lum}}$ &  $\\bar{N}_{\\rm{fov}}$ & $R_{\\rm{vir}} \\, \\rm{[kpc]}$  &  & $\\bar{N}_{\\rm{lum}}$ & $\\bar{N}_{\\rm{fov}}$ & $R_{\\rm{vir}} \\, \\rm{[kpc]}$ & & $\\bar{N}_{\\rm{lum}}$  & $\\bar{N}_{\\rm{fov}}$ & $R_{\\rm{vir}} \\, \\rm{[kpc]}$ &  & $\\bar{N}_{\\rm{lum}}$ & $\\bar{N}_{\\rm{fov}}$ & $R_{\\rm{vir}} \\, \\rm{[kpc]}$   \\\\\n" )  # N Behroozi & N Sawala

        for name, dist, nsubs, pgt1, nstd, n10, n90, rvir, fov in zip(Names,Distances, Nsubs.T, np.array(Pgt1).T, Nstd.T, N10.T, N90.T, np.array(Rvir).T, np.array(Nfov).T)[24:]:
            f.write("\\hline\n")

            f.write("%s & %d    & $%.2f$ & %.2f & %d &    & $%.2f$ & %.2f & %d &    & $%.2f$ & %.2f & %d &     & $%.2f$ & %.2f & %d    \\\\\n" %(name, dist, nsubs[0],fov[0],rvir[0],    nsubs[1],fov[1],rvir[1],    nsubs[2],fov[2],rvir[2],  nsubs[3],fov[3],rvir[3]))    #, nsubs[4],nstd[4], nsubs[5],nstd[5] ))


        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\tablecomments{$\\bar{N}_{\\rm{fov}}$ indicates the mean number of luminous satellites with $M_* > 10^4$\msun within a field of view of radius $0.56^\\circ$ (corresponding to a footprint of equal area as the Solo Dwarfs Project) centered on target Local Group dwarf galaxies, as predicted with various AM models. Galaxies with a larger heliocentric distance and smaller AM model inferred virial radius, $R_{\\rm{vir}}$, will have a larger fraction of their volume surveyed in the field of view. The total mean number of satellites within each galaxy's virial volume, $\\bar{N}_{\\rm{lum}}$, is listed for comparison.}\n")
        f.write("\\end{table}\n")




def add_sum_distr_ax(ax,min_lum,subset=True,fixedDM=False, re=False):
    if subset:
        StellarMassesSubset = StellarMasses[-5:]
    else:
        StellarMassesSubset = StellarMasses

    import scipy.misc

    for model in [am.GarrisonKimmel(reionization=re), am.GK16_grow(reionization=re),am.Moster(reionization=re),am.Brook(reionization=re)]:  #am.Sawala(),am.GarrisonKimmel16(reionization=re)
        if fixedDM:  # fix the dm mass to just moster
            tmpmodel = am.Moster()
            dm_masses = tmpmodel.mstar_to_mhalo(StellarMassesSubset,a=1)  #stellar_to_halo_mass(StellarMassesSubset,a=1.0)
            if subset:
                print np.log10(dm_masses), 'moster masses'
                print model.mstar_to_mhalo(StellarMassesSubset,a=1.0), model.label, 'masses'
        else:
            dm_masses = model.mstar_to_mhalo(StellarMassesSubset,a=1) #model.mstar_to_mhalo(StellarMassesSubset,a=1)  #stellar_to_halo_mass(StellarMassesSubset,a=1.0)
        print model.label
        if model.isPoisson:
            mean = model.get_field_total(dm_masses,min_lum)
            nocc = np.arange(0,mean*3+1)
            prob = (mean**nocc * np.e**-mean /scipy.misc.factorial(nocc))
            tmp=np.repeat(nocc,2)
        else:
        # for monte carlo methods, need to do 10,000 instances, looping over each halo and collecting all 10,000 numbers
        # then loop over all halos, and keep adding to the running total list.
        # then make histogram of the final 10,000 long distribution, normalized to 1
            samples = model.get_field_distr(dm_masses,min_lum,N=N_iter)
            distr,nocc = np.histogram(samples,bins=np.arange(min(samples),max(samples)+2))
            prob = distr/float(len(samples))
            tmp=np.repeat(nocc[0:-1],2)  # last value of bins is not inclusive
        ax.plot(np.append(tmp[1:],tmp[-1]+1), np.repeat(prob,2),linewidth=linewidth,color=model.color,label=model.label)

        #if isinstance(model, am.Moster):
        #    samples = model.get_field_distr(dm_massesFull,min_lum,N=N_iter)
        #    distr,nocc = np.histogram(samples,bins=np.arange(min(samples),max(samples)+2))
        #    prob = distr/float(len(samples))
        #    tmp=np.repeat(nocc[0:-1],2)  # last value of bins is not inclusive
        #    dall, = ax.plot(np.append(tmp[1:],tmp[-1]+1), np.repeat(prob,2),linewidth=linewidth,color='grey',linestyle='--')
    return


# want to generate mean and variance of the sum of all distributions
#for non monte carlo ones, just get the sum total, and plot poisson distribution
def plot_sum_distr(subset=True):
    plt.figure()
    ax = plt.subplot(111)
    add_sum_distr_ax(ax,min_lum=10**3,subset=True)

    plt.ylabel('Probability',fontsize=label_font)
    plt.xlabel('N satellites total',fontsize=label_font)
    plt.legend(loc='upper right',frameon=False,fontsize=legend_size)
    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.gcf().subplots_adjust(bottom=0.15)
    
    if subset:
        plt.savefig('largest10_field_sats')
    else:
        plt.savefig('total_field_sats')



# determine dark matter mass instead of stellar mass from just Moster.
# what are those masses? 
# refine field halos used for SHMF to those field halos that most closely mirror these 10 largest
# make a version of the plot below with the mass fixed to just that of Moster

def plot_sum_distr_3panel(subset=True,fixedDM=False,re=False):
    #w,h = plt.figaspect(1.5)
    f=plt.figure()
    h=f.get_figheight()
    w=f.get_figwidth()
    plt.close()

    fig = plt.figure(figsize=(w,h*1.5))
    ax1 =fig.add_subplot(3,1,1)
    ax2 =fig.add_subplot(3,1,2)
    ax3 =fig.add_subplot(3,1,3)
    #f,(ax1,ax2,ax3) = plt.subplots(nrows=3)

    add_sum_distr_ax(ax1,min_lum=10**3,subset=subset,fixedDM=fixedDM,re=re)
    add_sum_distr_ax(ax2,min_lum=10**4,subset=subset,fixedDM=fixedDM,re=re)
    add_sum_distr_ax(ax3,min_lum=10**5,subset=subset,fixedDM=fixedDM,re=re)
    
    ax1.set_ylabel('Probability',fontsize=label_font)
    ax2.set_ylabel('Probability',fontsize=label_font)
    ax3.set_ylabel('Probability',fontsize=label_font)
    ax3.set_xlabel('N satellites total',fontsize=label_font)


    #ax2.legend([dall], ['Moster + Reion \nAll Field Galaxies'],loc='upper right',frameon=False,fontsize=legend_size-3)

    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    ax3.tick_params(axis='both', which='major', labelsize=tick_size)
    #plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(bottom=0.10)
    #plt.gcf().subplots_adjust(top=0.05)

    ax1.text(.35,.82,'$M_* > 10^3 \, \mathrm{M_\odot}$',transform=ax1.transAxes,fontsize=legend_size)
    ax2.text(.35,.8,'$M_* > 10^4 \, \mathrm{M_\odot}$',transform=ax2.transAxes,fontsize=legend_size)
    ax3.text(.35,.8,'$M_* > 10^5 \, \mathrm{M_\odot}$',transform=ax3.transAxes,fontsize=legend_size)
    
    if subset:
        ax1.set_xlim((0,45))
        ax2.set_xlim((0,35))
        ax3.set_xlim((0,20))
        ax1.legend(loc='upper right',frameon=False,fontsize=legend_size-3)
    else:
        ax3.legend(loc='upper right',frameon=False,fontsize=legend_size-3) 

    extra=''
    if fixedDM:
        extra = 'fixed_DM_mass'
    if re:
        extra+='reionization'
    if subset:
        plt.savefig('largest10_field_sats_3panel'+extra+'.pdf')
    else:
        plt.savefig('total_field_sats_3panel'+extra+'.pdf')





def plot_sum_distr_2panel(subset=True,fixedDM=False):
    f,(ax1,ax2) = plt.subplots(nrows=2)
    add_sum_distr_ax(ax1,min_lum=10**3,subset=subset,fixedDM=fixedDM)
    add_sum_distr_ax(ax2,min_lum=10**4,subset=subset,fixedDM=fixedDM)
    
    ax1.set_ylabel('Probability',fontsize=label_font)
    ax2.set_ylabel('Probability',fontsize=label_font)
    ax2.set_xlabel('N satellites total',fontsize=label_font)
    ax1.legend(loc='upper right',frameon=False,fontsize=legend_size)
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.gcf().subplots_adjust(bottom=0.15)

    ax1.text(.45,.85,'$M_* > 10^3 \, \mathrm{L_\odot}$',transform=ax1.transAxes,fontsize=legend_size)
    ax2.text(.45,.85,'$M_^* > 10^4 \, \mathrm{L_\odot}$',transform=ax2.transAxes,fontsize=legend_size)
    
    extra=''
    if fixedDM:
        extra = 'fixed_DM_mass'
    if subset:
        plt.savefig('largest10_field_sats_2panel'+extra+'.pdf')
    else:
        plt.savefig('total_field_sats_2panel'+extra+'.pdf')



# FIGURE 6
#plot_sum_distr_3panel(subset=True,re=True)
#plot_sum_distr_3panel(subset=False,re=True)

# TABLE 2
#generate_latex_table(DwarfNames)
generate_latex_table2(DwarfNames)




#plot_sum_distr()
#plot_sum_distr_2panel(subset=True)
#plot_sum_distr_2panel(subset=False)
#plot_sum_distr_2panel(subset=True,fixedDM=True)  # fixed DM doesn't really work, it must be changed to 350 or 200
#plot_sum_distr_2panel(subset=False,fixedDM=True)



        #f.write("\\textbf{Name} &  \\boldmath$M_* \\, [10^6 \\, \\mathrm{M_\odot}]$ & \\multicolumn{2}{c}{\\textbf{Moster}}  & \\multicolumn{2}{c}{\\textbf{GK14}} &  \\multicolumn{2}{c}{\\textbf{GK16}} & \\multicolumn{2}{c}{\\textbf{Brook}} \\\\\n" )  # N Behroozi & N Sawala


        #f.write(" &  & $\\bar{N}_{\\rm{lum}}$ &  $P( \\geq 1)$  & $\\bar{N}_{\\rm{lum}}$ & $P( \\geq 1)$ & $\\bar{N}_{\\rm{lum}}$  & $P( \\geq 1)$ & $\\bar{N}_{\\rm{lum}}$ & $P( \\geq 1)$ \\\\\n" )  # N Behroozi & N Sawala

        ## kpc is on a different line than Rvir
        #f.write(" &  & $\\bar{N}_{\\rm{lum}}$ &  $P( \\geq 1)$ & $R_{\\rm{vir}}$  &  & $\\bar{N}_{\\rm{lum}}$ & $P( \\geq 1)$ & $R_{\\rm{vir}}$ & & $\\bar{N}_{\\rm{lum}}$  & $P( \\geq 1)$ & $R_{\\rm{vir}}$ &  & $\\bar{N}_{\\rm{lum}}$ & $P( \\geq 1)$ & $R_{\\rm{vir}}$   \\\\\n" )  # N Behroozi & N Sawala
        #f.write("  & & & &  $\\rm{[kpc]}$ & & & & $\\rm{[kpc]}$ & & & & $\\rm{[kpc]}$ & & & &   \\\\\n  ")



#f.write("%s & %.2f      & %.2f $\pm$ %.2f & %.2f & %d &    & %.2f $\pm$ %.2f & %.2f & %d &    & %.2f $\pm$ %.2f & %.2f & %d &     & %.2f $\pm$ %.2f & %.2f & %d    \\\\\n" %(name,mstar*1e-6, nsubs[0],nstd[0],pgt1[0],rvir[0],    nsubs[1],nstd[1],pgt1[1],rvir[1],    nsubs[2],nstd[2],pgt1[2],rvir[2],  nsubs[3],nstd[3],pgt1[3],rvir[3]))    #, nsubs[4],nstd[4], nsubs[5],nstd[5] ))

#f.write("%s & %.2f & %.2f $\pm$ %.2f & %.2f   & %.2f $\pm$ %.2f & %.2f    & %.2f $\pm$ %.2f & %.2f   & %.2f $\pm$ %.2f & %.2f \\\\\n" %(name,mstar*1e-6, nsubs[0],nstd[0],pgt1[0],    nsubs[1],nstd[1],pgt1[1],    nsubs[2],nstd[2],pgt1[2],  nsubs[3],nstd[3],pgt1[3]))    #, nsubs[4],nstd[4], nsubs[5],nstd[5] ))

"""
NGC 6822 & 100.00 & 3.88 & 60.94 & 3.24 & 2.27 & 11.41 \\
\hline
IC 1613 & 100.00 & 3.88 & 60.94 & 3.24 & 2.27 & 11.43 \\
\hline
IC 4662 & 190.00 & 5.12 & 97.71 & 4.43 & 3.01 & 15.40 \\
\hline
IC 5152 & 270.00 & 5.96 & 125.09 & 5.18 & 3.51 & 17.86 \\
"""







# Leo P: 1.62 Mpc or 1620 kpc. Stellar mass is 5.6*10^5
# KKs3  - stellar mass is 2.3e7. or 23e6.  distances is 2.12 Mpc or 2120 kpc
# KK258  - 2230 kpc, stellar mass = ??     . from discovery paper, 2230 kpc. stellar mass is 10**7.16 or 14e6


# total stelar mass formed during life is 2.2e7
# Perseus - distances is 785, stellar mass = ??   Also known as Perseus I. likely a distant satellite of M31. No stellar mass estimates yet.
# tucana is 10**6.62 in Karach but .56e6
# HIZSS3 A(B) - 1675 kpc. binary galaxy system, only HI masses measured. unable to include.

