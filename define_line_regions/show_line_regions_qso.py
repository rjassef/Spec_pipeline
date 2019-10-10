#!/usr/bin/env python 

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u

from Spec_pipeline.Line_Fitter.multi_line import Multi_Line_fit

###

def plot_vline(x,fmt):
    ymin, ymax = plt.ylim()
    plt.plot([x,x],[ymin,ymax],fmt)

###

#Read the vanden Berk template. File is in Flam. 
qso = Table.read("vandenberk_composite.txt",format='ascii.cds')
lam_qso = qso['Wave'].to(u.AA)
flam_qso = qso['FluxD'] * u.erg/u.s/u.cm**2/u.AA

#Emission lines
lines = [
    'Lyb_OVI',
    'LyA_NVred',
    'CIV',
    'HeII',
    'CIII',
    'CII',
    'NeIV',
    'MgII',
    'OII',
    'NeIII',
    'Hb',
    'OIII'
]


k = 0
while k<len(lines):

    print(lines[k])

    #Create the line object.
    line = Multi_Line_fit(lines[k],lines_file="multi_lines.txt")

    #Set x limits to +/-500 A
    xmin = np.min(line.line_center) - 500.*u.AA
    xmax = np.max(line.line_center) + 500.*u.AA
    plt.xlim([xmin.value,xmax.value])

    #Plot spectrum.
    kuse = np.where((lam_qso>=xmin) & (lam_qso<=xmax))
    ymin = 0.8*np.min(flam_qso[kuse])
    ymax = 1.2*np.max(flam_qso[kuse])
    plt.ylim([ymin.value,ymax.value])
    plt.plot(lam_qso,flam_qso,color='xkcd:grey',linewidth=1.0,
             linestyle='solid')
    
    for lam_center in line.line_center.value:
        plot_vline(lam_center,'b--')

    for lam_cont in line.continuum_regions.value.flatten():
        plot_vline(lam_cont,'g--')

    plt.show(block=False)

    answer = ""
    while answer!="yes" and answer!="no":
        print("Replot? (yes/no)")
        answer = input()
            
    if answer=='no':
        k+=1

    plt.close()
