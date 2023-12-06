"""
@author: Apelova 

A few examples that show how to use the egs_doselib!
In Order to analyse them comment change the if Statements of each block.
"""
#--- Import all classees and methods from egs_doselib!
from egs_doselib import *

#--- Import data
import matplotlib.pyplot as plt #to visualize data
import pandas as pd #to store data
import os #used for setting the paths

# TODO! 
#--- Change this to the path were you downloaded the repository too
origin = os.getcwd() + "/DATA_EXAMPLES/"

dose = dose_3d(origin+"/3ddose/Complete_dosxyz.3ddose", INFO=False) #INFO=True shows information on reading the data
###############################################################################
#--- Example 1 Read a 3ddose - File output from DOSXYZNRC
###############################################################################
if False:
    #alternatively one can print the dose_3d instance ! uncomment the line below to test this!
    #print(dose)
    #access the Default Percentage Depth Dose at (X,Y) = (0,0) and X/Y-Profiles !
    fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(8,6), sharex="col")
    #--- absolute PDD (top left)
    axs[0,0].plot(dose.position.z, dose.pdd, c="red")
    axs[0,0].set_ylabel("absolute Dose in [Gy]")
    axs[0,0].grid(True)
    #--- normalized PDD (bottom left)
    axs[1,0].plot(dose.position.z, norm(dose.pdd), c="red")
    axs[1,0].set_ylabel("relative Dose [%]")
    axs[1,0].set_xlabel("Position along Z-Axis [cm]")
    axs[1,0].grid(True)
    #--- absolute Profiles (top right)
    axs[0,1].plot(dose.position.x, dose.x_profile, label="X-Profile", c="blue")
    axs[0,1].plot(dose.position.y, dose.y_profile, label="Y-Profile", c="lightgreen")
    axs[0,1].set_ylabel("absolute Dose [Gy]")
    axs[0,1].legend(loc="lower center")
    axs[0,1].grid(True)
    #--- normalized Profiles (bottom rigt)
    axs[1,1].plot(dose.position.x, norm(dose.x_profile), label="X-Profile", c="blue")
    axs[1,1].plot(dose.position.y, norm(dose.y_profile), label="Y-Profile", c="lightgreen")
    axs[1,1].set_ylabel("relative Dose [%]")
    axs[1,1].set_xlabel("Position Perpendicular to Z-Axis [cm]")
    axs[1,1].legend(loc="lower center")
    axs[1,1].grid(True)
    fig.tight_layout()

###############################################################################
#--- Example 2  Change the class attributes .pdd, .x_profile, .y_profile
###############################################################################
if False:
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(18,5))
    if False: #Change Position of PDD
        # add unchanged PDD to graph at X,Y=0,0
        axs[0].plot(dose.position.z, dose.pdd, label="PDD through (X,Y)=(0,0)", c="blue")    
        # change the pdd to a different locationin the XY Plane and add to the plane!
        dose.set_pdd(X=5, Y=5, MUTE=True)
        axs[0].plot(dose.position.z, dose.pdd, label="PDD through (X,Y)=(5,5)", c="red")    
        #--- cometics for the plot    
        axs[0].legend(loc="upper right")
        
    if False: #Set profiles together
        #do the same with an x and y profile !
        axs[1].plot(dose.position.x, dose.x_profile, label=f"Profile at {dose.profile_depth_x}", c="blue", ls="dotted", lw=2)
        axs[1].plot(dose.position.y, dose.y_profile, label=f"Profile at {dose.profile_depth_y}", c="lightgreen", ls="dotted", lw=2)
        #--- set the new profiles at the same time
        dose.set_profiles(10.25, MUTE=True)
        axs[1].plot(dose.position.x, dose.x_profile, label=f"Profile at {dose.profile_depth_x}", c="blue", lw=2)
        axs[1].plot(dose.position.y, dose.y_profile, label=f"Profile at {dose.profile_depth_y}", c="lightgreen", lw=2)
        #--- cometics for the plot
        axs[1].legend(loc="lower center")        
    
    if False: #Set profiles independenntly
        axs[2].plot(dose.position.x, dose.x_profile, label=f"Profile at {dose.profile_depth_x}", c="blue", ls="dotted", lw=2)
        axs[2].plot(dose.position.y, dose.y_profile, label=f"Profile at {dose.profile_depth_y}", c="lightgreen", ls="dotted", lw=2)
        dose.set_x_profile(Z= 15, Y=0)
        dose.set_y_profile(Z= 7, X=0)
        axs[2].plot(dose.position.x, dose.x_profile, label=f"Profile at {dose.profile_depth_x}", c="blue", lw=2)
        axs[2].plot(dose.position.y, dose.y_profile, label=f"Profile at {dose.profile_depth_y}", c="lightgreen", lw=2)
        axs[2].legend(loc="lower center")        
        
    if False: # gather Inforamtion about available Profiles
        print("Available Depths for X-direction:\n", dose.available_depths_x,"\n")
        print("Available Depths for Y-direction:\n", dose.available_depths_y)
###############################################################################    
#--- Example 3 Extracting Slices of a Plane !
###############################################################################
if False:    
    fig, axs = plt.subplots(ncols=2,nrows=3, figsize=(4,6))
    #--- first and last XY-PLANE
    axs[0,0].imshow(dose.get_plane(AXIS="Z", POSITION_ON_AXIS=dose.position.z[0] ))
    axs[0,1].imshow(dose.get_plane(AXIS="Z", POSITION_ON_AXIS=dose.position.z[-1]))
    
    #--- first and last YZ-PLANE
    axs[1,0].imshow(dose.get_plane(AXIS="X", POSITION_ON_AXIS=dose.position.x[0] ))
    axs[1,1].imshow(dose.get_plane(AXIS="X", POSITION_ON_AXIS=dose.position.x[-1]))
    #--- first and last XZ-PLANE
    axs[2,0].imshow(dose.get_plane(AXIS="Y", POSITION_ON_AXIS=dose.position.y[0] ))
    axs[2,1].imshow(dose.get_plane(AXIS="Y", POSITION_ON_AXIS=dose.position.y[-1]))
    fig.tight_layout()

    # display axis as help
    # recall that all positions go from minmus to plus !
    #the arrows point from minus to plus !
    axs[0,0].set_xlabel("X \u2192", fontsize=20)
    axs[0,0].set_ylabel("Y \u2192", fontsize=20)
    axs[0,1].set_xlabel("X \u2192", fontsize=20)
    #
    axs[1,0].set_xlabel("\u2190 Y", fontsize=20)
    axs[1,0].set_ylabel("Z \u2192", fontsize=20)
    axs[1,1].set_xlabel("\u2190 Y", fontsize=20)
    #
    axs[2,0].set_xlabel("X \u2192", fontsize=20)
    axs[2,0].set_ylabel("Z \u2192", fontsize=20)
    axs[2,1].set_xlabel("X \u2192", fontsize=20)
    
    #get riff of all ticks
    for i in range(3):
        for j in range(2):
            axs[i,j].set_xticks([])
            axs[i,j].set_yticks([])

###############################################################################
#--- Example 4 Getting values at a specific Voxel
###############################################################################
if False:
    if False:
        matrix = dose.dose_matrix
        print(f"Shape of the Matrix: {matrix.shape} = ({len(dose.position.z)}, {len(dose.position.y)}, {len(dose.position.x)})")
        # Slice matrix by using the np.array indexing style Dose(x,y,z) == matrix[z_index, y_index, x_index]
        # to get the Dose at the Point (X,Y,Z) = ? one can either use the indices or the 
        print(f"Using indices:\
       D(x= 2.75, Y=-7.5, Z=-13.5) ==> dose_matrix[5, 15, 3] = {matrix[5, 15, 3]}")
        #
        z_index=dose.find_closest_index(dose.position.z, 2.75)
        x_index=dose.find_closest_index(dose.position.z, -7.5)
        y_index=dose.find_closest_index(dose.position.z, -13.5)
        print(f"Using Values:\
        D(x= 2.75, Y=-7.5, Z=-13.5) ==> dose_matrix[z_index, z_index, z_index] = {matrix[5, 15, 3]}")
        #####
        # This works one error_matrix as well
        #####
        
    # one can easily obtain a pdd by using slicing techniques
    if False:
       plt.plot(matrix[:, 30, 30], c = "lightgreen")
       plt.plot(dose.pdd, ls="dashed", c="black")

###############################################################################
#--- Example 6 Get PDD or Profile with Error in Percent
###############################################################################
if True:
    dose_values, error_in_percent = dose.get_pdd()
    dose_values, error_in_percent = dose.get_x_profile()
    dose_values, error_in_percent = dose.get_pdd()
    

###############################################################################
#
#                                  END 
#       all of the above features are the fundamentals of egs_doselib.
#     Everything below are just usefull things I needed to use for my work.           
#
###############################################################################


###############################################################################
#--- Example 6 Combine Segmented 3ddose files from EGS_CHAMBER_sepperate_profiles_template.egsinp (from /3ddose/)
###############################################################################
if False:
    dose = dose_3d(origin+"/3ddose/Sepperate_dose_pdd.3ddose", INFO=True)
    # only a pdd is available in this data-set
    print(f"\
It is easy to see that the dose-file only contains a PDD by regarding the voxel bounadries.\
or by testing the length of the profiles !\n\
 - The PDD has {len(dose.pdd)} elements \n\
 - The X-profile has {len(dose.x_profile)} element \n\
 - The Y-profile has {len(dose.y_profile)} element \n\
")
    # add the profiles
    dose.add_profile(origin+"/3ddose/Sepperate_dose_x_profile.3ddose", AXIS="X")
    dose.add_profile(origin+"/3ddose/Sepperate_dose_y_profile.3ddose", AXIS="Y")
    #--- show the availale profile depths
    print(f"X-Profile available at {dose.available_depths_x.tolist()}")
    print(f"Y-Profile available at {dose.available_depths_y.tolist()}")
    # --- 
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(6,3))
    axs[0].plot(dose.position.z, dose.pdd, c="red")
    axs[1].plot(dose.position.x, dose.x_profile, c="blue")
    axs[1].plot(dose.position.y, dose.y_profile, c="lightgreen")
    # Technically its possible to add multiple x/y profiles (different dephts) and load them with
    
###############################################################################
#--- Example 7 Read a .mcc-File from Meassured Data
###############################################################################

if False:
    dose = dose_mcc(origin+"/mcc/pdd_and_profiles_at_2_depths.mcc", INFO=True)
    fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10,5))
    axs[0].plot(dose.position.z, dose.pdd, c="red")
    axs[1].plot(dose.position.x, dose.x_profile, c="blue")
    axs[1].plot(dose.position.y, dose.y_profile, c="lightgreen")
    dose.set_profiles(DEPTH=14.15, MUTE=False)
    axs[1].plot(dose.position.x, dose.x_profile, c="blue")
    axs[1].plot(dose.position.y, dose.y_profile, c="lightgreen")
    #--- the pth pdd are identical because it is an artificially created .mcc file !

###############################################################################
#                                 TODO
#        I encourage you to play around with the example data given.
#              For more Information Please view the Sourcecode.
#    If there are any Bugs please report the issue at the GitHub repository
#               https://github.com/Apelova/EGS_DOSE_TOOLS
###############################################################################