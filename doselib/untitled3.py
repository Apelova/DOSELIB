"""
@author: apel04
email: marvin.apel.ext@ptb.de

path = "/home/marvin/Desktop/master"     # Notebook
path = "C:/Users/apel04/Desktop/master/" # Office Pc
"""
import matplotlib.pyplot as plt
from doselib import *
import pandas as pd 
import numpy as np

def norm(arr):
    return np.array(arr/max(arr)*100)

path = "/home/marvin/Desktop/master"     # Notebook
path = "C:/Users/apel04/Desktop/master/" # Office Pc


def gamma(axes_reference, dose_reference, axes_evaluation, dose_evaluation,
             dose_percent_threshold=3, distance_mm_threshold=0.3, dx=0.1, dD=0.1,
             fill_step=0.01, CUTOFF=20, show_plt=False, derivatives=False):
    """ Assumptions:
            - both profiles are on a 0.1mm scaled grid ! -> could interpolate if necessary not implemented yet
            - only tested for 1d Profiles !
    """
    # Required Assisting functions !
    #--------------------------------------------------------------------------
    def interpol_position(axis, step=0.05):
        x_grid = np.arange(left_lim, right_lim, step=step)
        f_interpol_position = interpolate.interp1d(axis, axis, fill_value='extrapolate', kind="linear")
        return f_interpol_position(x_grid)

    def interpol_dose(axis, dose, step=0.05):
        x_grid = np.arange(left_lim, right_lim, step=step)
        f_interpol_position = interpolate.interp1d(axis, dose, fill_value='extrapolate', kind="linear")
        return f_interpol_position(x_grid)

    def euclid_distance_gamma(x_r, D_r, x_e , D_e, dx=0.3, dD=3):
        return np.sqrt(((D_r-D_e)/dD)**2 + ((x_r-x_e)/dx)**2)
    
    # Step 0a get limits for which Dose > 20 in both profiles ! 
    #--------------------------------------------------------------------------
    a, b = axes_reference[dose_reference>=CUTOFF][[0,-1]]
    c, d = axes_evaluation[dose_evaluation>=CUTOFF][[0,-1]]
    a, b, c, d = round(a, 2), round(b, 2), round(c, 2), round(d, 2) #to avoid floating point errors
    left_lim, right_lim = max([a,c]), min([b,d])
    gamma_array = np.zeros(len(axes_reference))
    
    try: #if they are not on the same grid switch them 
        ix = axes_reference.tolist().index(left_lim) #used to fill array
    except:
        raise ValueError("Try to Switch the Input doseobjects! If this doesnt the error most likely occurs to the limted resolution of 1 digit after '.' in the axes.")
        
    del a, b, c, d
    
    if show_plt:    
        axes_for_plot = axes_reference

    # Step 0b select observed Data (!!!!! Maintain order in selecting the indices ! first dose then axes !!!!!)
    #--------------------------------------------------------------------------
    dose_reference  = dose_reference[(axes_reference>left_lim)&(axes_reference<right_lim)]
    axes_reference  = axes_reference[(axes_reference>left_lim)&(axes_reference<right_lim)]
    dose_evaluation = dose_evaluation[(axes_evaluation>left_lim)&(axes_evaluation<right_lim)]
    axes_evaluation = axes_evaluation[(axes_evaluation>left_lim)&(axes_evaluation<right_lim)]

    #Plot
    if show_plt: #shows that a fine interpolation is required !  
        if derivatives:
            fig, axs = plt.subplots(1,4, figsize=(36, 10))    
        else:
            fig, axs = plt.subplots(1,3, figsize=(30, 10))   
        axs[0].scatter(axes_reference, dose_reference, ec="black", fc="white", s=50)
        axs[0].scatter(axes_evaluation, dose_evaluation, marker="x", c="red", s=50) 
        
    # Step 1a - interpol dose on fine grid !
    #--------------------------------------------------------------------------
    dose_evaluation_ip= interpol_dose(axes_evaluation, dose_evaluation, step=fill_step)
    
    # Step 1b - interpol position on fine grid !
    #--------------------------------------------------------------------------
    axes_evaluation_ip = interpol_position(axes_evaluation, step=fill_step)
    #Plot
    if show_plt:
        axs[1].scatter(axes_evaluation_ip, dose_evaluation_ip, marker="x", c="red", s=50) 
        axs[1].scatter(axes_reference, dose_reference, ec="black", fc="white", s=50) 

    # Step 2 Implement Brute-Force Gamma-Index Calculation
    #TODO! additionally use numba !    
    count = 0
    #--------------------------------------------------------------------------        
    for x_R, D_R in zip(axes_reference, dose_reference):
        min_gam = np.inf
        for i in range(len(axes_evaluation_ip)):
            temp = euclid_distance_gamma(
                                         x_r = x_R , 
                                         D_r = D_R ,
                                         x_e = axes_evaluation_ip[i] , 
                                         D_e = dose_evaluation_ip[i],
                                         dx=dx, dD=dD
                                        )
            min_gam = temp if temp < min_gam else min_gam
            
    # Step 3 Fill the zero-array such that gamma is only available at x where D(x) > 20
    #--------------------------------------------------------------------------    
        gamma_array[ix+count] = min_gam
        count +=1
    del count


            
    return gamma_array  


def plot(reference, to_be_evaluated):
    gamma_x = gamma(axes_reference = reference.position.x, dose_reference  = reference.x_profile.norm, 
                     axes_evaluation= to_be_evaluated.position.x, dose_evaluation = to_be_evaluated.x_profile.norm, 
                     dD=dD, dx=dx) #3% und 3mm
     #Y-Axis
     gamma_y = gamma(axes_reference = reference.position.y, dose_reference  = reference.y_profile.norm, 
                     axes_evaluation= to_be_evaluated.position.y, dose_evaluation = to_be_evaluated.y_profile.norm, 
                     dD=dD, dx=dx) #3% und 3mm
    
    
     Gamma_PR_x = round(len(gamma_x[gamma_x<1])/len(gamma_x)*100,2)
     Gamma_PR_y = round(len(gamma_y[gamma_y<1])/len(gamma_y)*100,2)




plot(dose_exp, dose_PP)




# %%
dose_03 = dose_3d("C:/Users/apel04/Desktop/Data_Backup/MASTER_OLD/6MV/TEST_Voxel_Volume/final_analysis/correct_ssd/0_3x0_3x0_3/x_profile_voxel_03_03_03.3ddose")
dose_03.add_profile("C:/Users/apel04/Desktop/Data_Backup/MASTER_OLD/6MV/TEST_Voxel_Volume/final_analysis/correct_ssd/0_3x0_3x0_3/y_profile_voxel_03_03_03.3ddose", "y")#---
dose_PinPoint = dose_chamber_log([
    "C:/Users/apel04/Desktop/Data_Backup/MASTER_OLD/6MV/TEST_Voxel_Volume/final_analysis/correct_ssd/PinPoint/PinPoint31022_simulation_x_complete.egslog",
    "C:/Users/apel04/Desktop/Data_Backup/MASTER_OLD/6MV/TEST_Voxel_Volume/final_analysis/correct_ssd/PinPoint/PinPoint31022_simulation_y_complete.egslog"
    ], delete_deviants=True, error_limit=2)
#---
dose_exp_pp = dose_mcc("N:/STRTHP/Messungen_LINAC_151605/6MeV_FF_complete/PP_PTW31021_MESSUNG/PROFILES_FS_10x10_SSD_100_dx_0_1_dt_02.mcc", average_profiles=True)
























