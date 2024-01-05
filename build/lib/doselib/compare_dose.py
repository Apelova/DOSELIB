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

######################
#
# Kann ich sicherstellen das die profile bei konstanter tiefe verglichen werden ?
#
######################


def compare_dose( dose_objects, labels, axes=["Z", "X", "Y"], difference=True, figsize=None, metrics=None, colors = ["black", "red", "blue"], interpol="linear", diff_dx=0.1):
    """
        A function that return plot of selected axes.
        The limits can be set by calling the output !
        
        If difference is True and (dose_objects) len(2) it assumes the difference to the first element of dose_objects should be calculated for the others
        
        metrics = "average", "profiles"
        
        TODO-> elute input
    """    
    def interpolate_arrays(position_array, dose_array, dx=diff_dx, interpol=interpol):
        #--- get position limits
        min_shared_position = max([min(dose) for dose in position_array])  #all intervalls include this point !
        max_shared_position = min([max(dose) for dose in position_array])  #all intervalls include this point !

        #--- set the grid on which interpolation will take place
        position_interpolated = np.arange(min_shared_position, max_shared_position+dx, dx)

        if position_interpolated[-1] > max_shared_position:#removes last element if its outside of data !
            position_interpolated = position_interpolated[:-1]
            
        #--- calculate the interpolated dose for each dose_array
        dose_arrays_interpolated = []
        for position, dose in zip(position_array, dose_array):
            #--- calculate interpol function
            interpol_functon = interpolate.interp1d(position, dose, kind=interpol)
            #--- apply function and store values
            dose_arrays_interpolated.append(interpol_functon(position_interpolated))

        #--- dose_arrays_interpolated contains the interpolated arrays in an array in the same order as it was given by dose_array and position_array !                
        return position_interpolated, dose_arrays_interpolated

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
    if difference:
        fig_, axs_ = plt.subplots(rows, cols, figsize=figsize, sharex="col", constrained_layout=True, gridspec_kw={'height_ratios': [4,1]})
    else:
        fig_, axs_ = plt.subplots(rows, cols, figsize=figsize, sharex="col", constrained_layout=True)


    #--- define helping functions !
    def plot_difference(fig_, axs_, position_dose_pair_arrays, labels, metrics=None, colors = ["black", "red", "blue"], difference=True):
        #--- group positions and pairs together and interpolate dose on a common grid !
        common_position, interpolated_dose = interpolate_arrays( [pair[0] for pair in position_dose_pair_arrays], [pair[1] for pair in position_dose_pair_arrays])
        #--- iterate through all the AXIS to plot !
        for i in range(1, len(interpolated_dose)):
            axs_[1,col].plot(common_position, interpolated_dose[i]-interpolated_dose[0], c=colors[i], lw=2)
            
        
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

    for col in range(cols):
        #--- print original plots
        if axes[col] == "X":
            plot_dose_distribution(fig_, axs_, [[each.position.x, each.x_profile.norm] for each in dose_objects], labels, difference=difference)
                
        if axes[col] == "Y":
            plot_dose_distribution(fig_, axs_, [[each.position.y, each.y_profile.norm] for each in dose_objects], labels, difference=difference)
    
        if axes[col] == "Z":
            plot_dose_distribution(fig_, axs_, [[each.position.z, each.pdd.norm] for each in dose_objects], labels, difference=difference)

    #--- add grid
    for i in range(cols):
        if difference:
            for j in range(rows):
                axs_[j,i].grid()
        else:
            axs_[i].grid()

    #--- add x-labels
    axes_labels = {"X": "Position along X-Axis [cm]", "Y": "Position along Y-Axis [cm]", "Z": "Position along Z-Axis [cm]"}

    for i, each in enumerate(axes):
        if difference:
            axs_[1,i].set_xlabel(axes_labels[each], fontsize=18)
        else:
            axs_[i].set_xlabel(axes_labels[each], fontsize=18)

    #--- add y-labels
    if difference:
        axs_[0,0].set_ylabel("relative Dose [%]", fontsize=18)
        axs_[1,0].set_ylabel("Δ [%]", fontsize=18)
        axs_[1,0].set
    else:
        axs_[0].set_ylabel("relative Dose [%]", fontsize=18)
        
    #--- return fig and axs
    return fig_, axs_

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




































