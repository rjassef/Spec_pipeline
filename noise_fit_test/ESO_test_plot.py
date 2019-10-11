#!/usr/bin/env python 

import numpy as np
import astropy.units as u
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import sys

from Fake_ESO_Spec import Fake_ESO_Spec

if len(sys.argv)!=2:
    print("Correct use: python",sys.argv[0],"gmag")
    sys.exit()
gmag = float(sys.argv[1])

z_values = np.arange(0.5,3.0+0.1,0.5)
moons = [0,90,180]
lam_targs = [4000.,6000.,8000.]

line_color=["blue","red","darkgreen"]
fill_color=["cyan","magenta","green"]


fig,ax = plt.subplots(len(moons),len(lam_targs),
                      sharex=True,sharey=True,
                      figsize=(18,18))
fig.subplots_adjust(hspace=0,wspace=0)

fig.suptitle("vanden Berk Composite, g={0:.1f} mag".format(gmag),
             fontsize=36)
fig.text(0.5,0.04,
         "Redshift"
         ,ha='center',fontsize=36)
fig.text(0.04,0.5,
         r'$\delta f_{\lambda}(\rm{Fit})/\delta f_{\lambda}(\rm{Input})$',
         va='center',rotation='vertical',fontsize=36)

for j,lam_targ in enumerate(lam_targs):
    for i,moon_template in enumerate(moons):

        ax[i,j].set_xlim([0,len(z_values)+1])
        ax[i,j].set_ylim([0.60,1.59])
        ax[i,j].tick_params(axis='both', which='major', labelsize=30)
        ax[i,j].tick_params(axis='both', which='minor', labelsize=20)

        for m,moon in enumerate(moons):

            frac_flam_err = []
            for z in (z_values):

                #First, create the fake object.
                fake_spec = Fake_ESO_Spec(gmag=gmag,z=z,moon=moon)

                #Now, parse the MC file. 
                MC_fname = "MCfiles/err_ESO.z{0:.1f}.g{1:.1f}.m{2:d}.mt{3:d}.dat".format(
                    z,gmag,moon,moon_template)
                try:
                    cat = open(MC_fname)
                    nrep, nlam = [int(float(ix)) for ix in cat.readline().split()]
                    cat.close()
                except IOError:
                    print("Could not open file",MC_fname)
                    sys.exit()

                data = np.loadtxt(MC_fname,skiprows=1)
                lam = data[:nlam]
                k = np.argmin(np.abs(lam-lam_targ))
                flam_err = data[nlam+k::nlam]

                frac_flam_err.append(flam_err/fake_spec.flam_err[k].value)

            #We'll make a boxplot. The box extends from the lower
            #to the upper quartile. We'll make the whiskers show
            #the entire range.
            bp = ax[i,j].boxplot(frac_flam_err,whis='range',labels=z_values,
                                 notch=True,patch_artist=True)
            for element in ['boxes', 'whiskers', 'fliers', 'means', 
                            'medians', 'caps']:
                plt.setp(bp[element], color=line_color[m])
            for patch in bp['boxes']:
                patch.set_facecolor(fill_color[m])
                patch.set_alpha(0.3)
                
            if i==0:
                ax[i,j].set_title(
                    r'$\lambda={0:.0f}$\AA'.format(lam_targ),fontsize=36)
                
            if j==0:
                if moon_template==0:
                    label = "New Moon fit"
                elif moon_template==90:
                    label = "Quarter Moon fit"
                else:
                    label = "Full Moon fit"
                ax[i,j].text(0.40,1.50,label,fontsize=20)
                
            if i==0 and j==0 and m==0:
                ax[i,j].text(0.4,1.42,"New Moon input",
                             color=line_color[0],fontsize=20)
                ax[i,j].text(0.4,1.34,"Quarter Moon input",
                             color=line_color[1],fontsize=20)
                ax[i,j].text(0.4,1.26,"Full Moon input",
                             color=line_color[2],fontsize=20)

#plt.show()
plt.savefig("ESO_test.g{0:.1f}.png".format(gmag),dpi=300)

