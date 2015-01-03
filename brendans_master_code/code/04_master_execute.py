import numpy as np
import subprocess as sub
import sys,os,glob
import mergertrees.MTCatalogue as MT
import readhalos.RSDataReader as rstar
import readsnapshots.readsnapHDF5_greg as rhdf5
import brendanlib.grifflib as glib
import pickle, platform

""" 
Contact: Brendan Griffen <brendan.f.griffen@gmail.com>

This code is made for running a simulation suite on LX11 caterpillar halos.
It creates initial conditions from MUSIC with different geometries and submits
Gadget jobs to the cluster. These geometries are in each folder as follows:

halos/H*/H*{geometry}*LX11*NVIR{nrvir}/contamination_suite/
geometry: set by one of the following options:
HX_BA_XX: Original MUSIC bounding box (e.g. the exact bounding box of lagr volume).
HX_BB_XX: 1.2 bounding box extent
HX_BC_XX: 1.4 bounding box extent
HX_BD_XX: 1.6 bounding box extent
HX_CA_XX: Convex Hull Volume
HX_EA_XX: Original MUSIC ellipsoid (e.g. the exact bounding box of lagr volume).
HX_EB_XX: 1.1 padding.
HX_EC_XX: 1.2 padding.
nrvir: set by input below.

Ideally, don't run without Brendan's permission as it carries out serious file alterations.
Once these halos have completed, runanyrockstar.py will take care of the Rockstar halo/MT generation.
See modules/brendanlib/grifflib.py for key functions.

Please comment out which piece of code you would like to run. 
Be warned: all mods using 'glib' will do something if you comment them out.

"""

if "antares" not in platform.node():
	print "YOU ARE NOT RUNNING FROM ANTARES - GOODBYE!"
	sys.exit()

nrvir = 5
lx_list = ["1"]  	    # Set which levels, you would like to iterate over (for contamination AND/OR full runs)
SUBMIT_GADGET = True 	# Set to false if you don't want to submit Gadget runs (will only copy files if False)
mass_bin = sys.argv[1]

base_path = "/bigbang/data/AnnaGroup/caterpillar/halos/"+mass_bin+"_mass_halos/"
music_path = "/bigbang/data/bgriffen/lib/music_pad/MUSIC"
lagr_path = "/bigbang/data/AnnaGroup/caterpillar/ics/lagr/"
gadget_file_path = "/home/bgriffen/exec/"

if not os.path.isdir(base_path):
    print "CAN NOT FIND BASE PATH:",base_path
    sys.exit()

if nrvir == 4:
    suite_names = ["EX"]

    #["BA","BB","EA","CA","EB","EC"]
if nrvir == 5:
    suite_names = ["EX"]

    #["EA","CA","BA"]

print "+EXECUTING +"

# MAKES ALL CONTAMINATION TEST FOLDERS
glib.make_destination_folders(base_path,suite_names,lx=11,nrvir=nrvir)

# CREATE SUITE LIST FOR CONTAMINATION STUDY ONLY (e.g. halos/.../H{haloid}/contamination_suite/)
contamination_paths = glob.glob(base_path + "H*/contamination_suite/*NV"+str(nrvir)+"*")
# Reminder: be sure to set lx_list correctly. "1" for LX11 etc.
#glib.run_music(contamination_paths,music_path,lagr_path,lx_list)

glib.run_gadget(contamination_paths,gadget_file_path,lx_list,submit=SUBMIT_GADGET)
# Subfind: has not been tested yet.
#glib.run_subfind(contamination_paths,gadget_file_path)

# CREATE HIGHER LEVELS (e.g. in halos/.../H{haloid}/), loads correct geometry from pickled dictionary.
halo_geometries =  pickle.load( open( base_path+"geometries.p", "rb" ) )
# Generates initial conditions for higher level once geomtry is set in geometries.p
#glib.run_music_higher_levels(halo_geometries,base_path,music_path,lagr_path,lx_list=lx_list)

# Runs through available halos in dictionary and runs P-Gadget3 on each halo.
#for halo_name,ic_info in halo_geometries.iteritems():
#    suite_paths = glob.glob(base_path+halo_name+"/H*")
#    glib.run_gadget(suite_paths,gadget_file_path,lx_list,submit=SUBMIT_GADGET)

print
print "------ CONTAMINATION RUNS ---------"
print "------   "+mass_bin+"_mass_halos/   -------"
glib.get_completed_list(contamination_paths)
print 
print "-------- PRODUTION RUNS -----------"
suite_paths = glob.glob(base_path + "H*/H*")
glib.get_completed_list(suite_paths)

print "+ DONE + "
