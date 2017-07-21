import DwarfMethods as dm
import os,sys,subprocess,time

hpaths = dm.get_hpaths(field=True,lx=14)

for hpath in hpaths:   
    subprocess.call("stat -c %y "+hpath+"/analysis/FieldHaloSubstructure.dat", shell=True)
    subprocess.call("stat -c %y "+hpath+"/analysis/AllExtantFieldData.dat", shell=True)



"""
stat -c %y /bigbang/data/AnnaGroup/caterpillar/halos/H1631506/H1631506_EA_Z127_P7_LN7_LX14_O4_NV5/analysis/FieldHaloSubstructure.dat
 '/bigbang/data/AnnaGroup/caterpillar/halos/H264569/H264569_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1725139/H1725139_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H447649/H447649_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H5320/H5320_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H581141/H581141_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H94687/H94687_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1130025/H1130025_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1387186/H1387186_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H581180/H581180_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1354437/H1354437_EA_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1725272/H1725272_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1195448/H1195448_EC_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1292085/H1292085_EX_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H796175/H796175_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H388476/H388476_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1079897/H1079897_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H94638/H94638_EX_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H95289/H95289_BB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1232164/H1232164_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1422331/H1422331_EX_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H196589/H196589_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1268839/H1268839_EB_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1599988/H1599988_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1195075/H1195075_EC_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1631582/H1631582_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1422429/H1422429_EX_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H65777/H65777_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1232423/H1232423_EA_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H196078/H196078_EA_Z127_P7_LN7_LX14_O4_NV5',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1599902/H1599902_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H795802/H795802_EX_Z127_P7_LN7_LX14_O4_NV4',
 '/bigbang/data/AnnaGroup/caterpillar/halos/H1104787/H1104787_EB_Z127_P7_LN7_LX14_O4_NV4']
"""
