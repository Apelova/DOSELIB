"""
@author: apel04
email: marvin.apel.ext@ptb.de

path = "/home/marvin/Desktop/master"     # Notebook
path = "C:/Users/apel04/Desktop/master/" # Office Pc
"""
import matplotlib.pyplot as plt
from egs_doselib import *
import pandas as pd
import numpy as np
import matplotlib

path = "/home/marvin/Desktop/master"    # Notebook
#path = "C:/Users/apel04/Desktop/master"  # Office Pc

dose_dosy = dose_3d(path+"/Simulationen/DATA/TEST_PTB/DOSXYZ_NRC/PTB_6MV_old_dosxyz.3ddose")
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
def compare_dose( dose_objects, labels, axes=["Z", "X", "Y"], difference=True, figsize=None, metrics=None, colors = ["black", "red", "blue"]):
    """
        A function that return plot of selected axes.
        The limits can be set by calling the output !
        
        If difference is True and (dose_objects) len(2) it assumes the difference to the first element of dose_objects should be calculated for the others
        
        metrics = "average", "profiles"
        
        TODO-> elute input
    """    
    #--- limited input testing
    if len([ each for each in axes  if each not in ["X","Y","Z"] ]) != 0:
        raise TypeError(f"Invalid input for Axis! THe following Axes are invalid {[ each for each in axes  if each not in ['X','Y','Z'] ]} Valid Inputs are [x1X, Y, Z] ")

    #--- set variables for layout
    cols, rows = (len(axes)), (int(difference)+1)

    #--- set figsize depending on requested graphs
    if not figsize:
        figsize = (8*(cols), (4*rows))
        figsize = (10*(cols), (5*rows))
    
    #--- add requested plots
    fig_, axs_ = plt.subplots(rows, cols, figsize=figsize, sharex="col", constrained_layout=True)

    #--- define helping functions !
    def plot_difference(fig_, axs_, position_dose_pair_arrays, labels, metrics=None, colors = ["black", "red", "blue"], difference=True):
        #TODO !Hier mit scipy alle auf die gleichen werte bringen -> linear interpolieren
        
        for i, pair in enumerate(position_dose_pair_arrays):
            # the first one is the reference -> make it appear different!
            if i==0:
                pass

        
    #--- define helping functions !
    def plot_dose_distribution(fig_, axs_, position_dose_pair_arrays, labels, metrics=None, colors = ["black", "red", "blue"], difference=True):
        #--- for each plot of the pairs 
        for i, pair in enumerate(position_dose_pair_arrays):
            # the first one is the reference -> make it appear different!
            if i==0:
                if difference:
                    axs_[0, col].scatter(pair[0], pair[1], c=colors[i], label=labels[i], marker="o", s=45, alpha=0.75)
                else:
                    axs_[col].scatter(pair[0], pair[1], c=colors[i], label=labels[i], marker="o", s=45, alpha=0.75)

            else:                                       
                if difference:
                    axs_[0, col].plot(pair[0], pair[1], c=colors[i], lw=2, label=labels[i])
                else:
                    axs_[col].plot(pair[0], pair[1], c=colors[i], lw=2, label=labels[i])

        if difference:
            plot_difference(fig_, axs_, position_dose_pair_arrays, labels, metrics, colors)

       
    # Das war nur ein test wie es klappt hier jetzt die angepassten dosiswerte nehmen     
    for col in range(cols):
        #--- print original plots
            match axes[col]: 
                case "X":
                    plot_dose_distribution(fig_, axs_, [[each.position.x, each.x_profile.norm] for each in dose_objects], labels, difference=difference)
                
                case "Y":
                    plot_dose_distribution(fig_, axs_, [[each.position.y, each.y_profile.norm] for each in dose_objects], labels, difference=difference)
        
                case "Z":
                    plot_dose_distribution(fig_, axs_, [[each.position.z, each.pdd.norm] for each in dose_objects], labels, difference=difference)


    # do some configuring/eye candy !
    
       
    return fig_, axs_
    
    
    
    
fig, axs = compare_dose( [dose_exp, dose_dosy, dose_chamber], labels=["Experimental","DOSXYZ_NRC", "EGS_CHAMBER"], axes=["Z", "X", "Y"], difference=True)

try:
    axs[0].legend(loc="upper right", fontsize=20)
except:
    axs[0,0].legend(loc="upper right", fontsize=20)

# %%

















