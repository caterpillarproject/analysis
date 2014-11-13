from brendanlib.grifflib import getlagrxyz,getcentext
from brendanlib.ellipsoid import EllipsoidTool as ET
import pylab as plt
import os
import numpy as np
import matplotlib.patches as patches
import glob,sys

"""
Contact: Brendan Griffen <brendan.f.griffen@gmail.com>

This script will plot xy,xz,yz projections of the lagrangian volume
including the bounding box from the header files.

"""
base_path = '/bigbang/data/AnnaGroup/caterpillar/ics/lagr/'
candidatelist = glob.glob(base_path + 'H*.head')

def plotboundingbox(xbox,ybox,zbox,extx,exty,extz,color,style,width):
    ax1.add_patch(patches.Rectangle((xbox,ybox),extx,exty,facecolor='none',
                                                         edgecolor=color,
                                                         linestyle=style,linewidth=width))
    ax2.add_patch(patches.Rectangle((xbox,zbox),extx,extz,facecolor='none',
                                                         edgecolor=color,
                                                         linestyle=style,linewidth=width))
    ax3.add_patch(patches.Rectangle((ybox,zbox),exty,extz,facecolor='none',
                                                         edgecolor=color,
                                                         linestyle=style,linewidth=width))

for pointfile in candidatelist:
    haloid = pointfile.split("NRVIR")[0].split("lagr/")[1]
    nrvir = pointfile.split("NRVIR")[1].replace(".head","")

    filename = './images/'+str(haloid)+'_NRVIR'+str(nrvir)+'_LAGR_PROJ.png'
    if os.path.isfile(pointfile) and not os.path.isfile(filename):
        print "reading...",pointfile
        x,y,z = getlagrxyz(pointfile.replace(".head",""))
        centx,centy,centz,extx,exty,extz = getcentext(pointfile)
        if not os.path.isfile(filename):
            fig = plt.figure(figsize=(17,5))
            ax1 = fig.add_subplot(131)
            ax2 = fig.add_subplot(132)
            ax3 = fig.add_subplot(133)
        
            ax1.plot(x,y,'b.',markeredgecolor='b')
            ax2.plot(x,z,'b.',markeredgecolor='b')
            ax3.plot(y,z,'b.',markeredgecolor='b')
        
            ax1.set_xlabel('x-pos')
            ax1.set_ylabel('y-pos')
            ax2.set_xlabel('x-pos')
            ax2.set_ylabel('z-pos')
            ax3.set_xlabel('y-pos')
            ax3.set_ylabel('z-pos')

            #extx = extx*2.0
            xbox = centx - 0.5*extx
            ybox = centy - 0.5*exty
            zbox = centz - 0.5*extz

            plotboundingbox(xbox,ybox,zbox,extx,exty,extz,'r','dotted',3)

            xbox = centx - 0.5*extx*1.2
            ybox = centy - 0.5*exty*1.2
            zbox = centz - 0.5*extz*1.2

            plotboundingbox(xbox,ybox,zbox,extx*1.2,exty*1.2,extz*1.2,'m','dashdot',3)

            xbox = centx - 0.5*extx*1.4
            ybox = centy - 0.5*exty*1.4
            zbox = centz - 0.5*extz*1.4

            plotboundingbox(xbox,ybox,zbox,extx*1.4,exty*1.4,extz*1.4,'g','dashed',3)

            xbox = centx - 0.5*extx*1.6
            ybox = centy - 0.5*exty*1.6
            zbox = centz - 0.5*extz*1.6

            plotboundingbox(xbox,ybox,zbox,extx*1.6,exty*1.6,extz*1.6,'k','dashed',3)

            ax1.set_xlim(x.min()*0.8,x.max()*1.2)
            ax1.set_ylim(y.min()*0.8,y.max()*1.2)
        
            ax2.set_xlim(x.min()*0.8,x.max()*1.2)
            ax2.set_ylim(z.min()*0.8,z.max()*1.2)
        
            ax3.set_xlim(y.min()*0.8,y.max()*1.2)
            ax3.set_ylim(z.min()*0.8,z.max()*1.2)
        
            fig.savefig(filename,bbox_inches='tight')
            plt.close(fig)


