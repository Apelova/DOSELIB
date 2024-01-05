"""
@author: apel04
email: marvin.apel.ext@ptb.de

path = "/home/marvin/Desktop/master"     # Notebook
path = "C:/Users/apel04/Desktop/master/" # Office Pc
"""
import matplotlib.pyplot as plt
from scipy import interpolate
from egs_doselib import *
import pandas as pd
import numpy as np
import matplotlib

path = "/home/marvin/Desktop/master"    # Notebook
path = "C:/Users/apel04/Desktop/master"  # Office Pc

dose_dosy = dose_3d(path+"/Simulationen/DATA/TEST_PTB/DOSXYZ_NRC/ROUGH_GRID/PTB_6MV_old_dosxyz.3ddose")
dose_dosy.set_x_profile(Z=5.25, Y=0)
dose_dosy.set_y_profile(Z=5.25, X=0)

#--- chamber-3ddose-File!
dose_chamber = dose_3d(path+"/Simulationen/DATA/TEST_PTB/EGS_CHAMBER/ROUGH_GRID/PTB_6MeVp_10x10_old_dose_pdd.3ddose")
dose_chamber.add_profile(path+"/Simulationen/DATA/TEST_PTB/EGS_CHAMBER/ROUGH_GRID/PTB_6MeVp_10x10_old_dose_x_profile.3ddose", AXIS="X")
dose_chamber.add_profile(path+"/Simulationen/DATA/TEST_PTB/EGS_CHAMBER/ROUGH_GRID/PTB_6MeVp_10x10_old_dose_y_profile.3ddose", AXIS="Y")

#--- experimentell-Data!
#dose_exp = dose_mcc(path+"/Messungen/6MeV_10x10_Dose_Profiles/profiles_@_5cm_ff_and_fff_SSD_100/FF/combined_ff.mcc")
dose_exp = dose_mcc(path+"/Messungen/6MeV_10x10_Dose_Profiles/profiles_@_5_25cm_ff_SSD_98/all_profiles_scan.mcc")

# %%    
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 
fig, axs = compare_dose( [dose_exp, dose_chamber, dose_dosy], labels=["Experimental","DOSXYZ_NRC", "EGS_CHAMBER"], axes=["Z", "X", "Y"], difference=True, interpol="quadratic", diff_dx=0.1)
axs[0,0].legend(loc="lower left", fontsize=30)

axs[0,0].set_xlim(0,30)
axs[0,1].set_xlim(-15,15)
axs[0,2].set_xlim(-15,15)
#first diff
axs[1,0].set_ylim(-10,10)
axs[1,0].set_yticks(np.arange(-10,15,5))
#second diff
axs[1,1].set_ylim(-10,10)
axs[1,1].set_yticks(np.arange(-10,15,5))
#third diff
axs[1,2].set_ylim(-20,20)
axs[1,2].set_yticks(np.arange(-20,25,10))

# Testing -> random dose, with and without difference, more axes !


# %%
fig, axs = compare_dose( [dose_exp, dose_dosy, dose_chamber], labels=["Experimental","DOSXYZ_NRC", "EGS_CHAMBER"], axes=["Y", "X"], difference=True, interpol="quadratic", diff_dx=0.01)




































