"""
@author: apel04
"""
import matplotlib.pyplot as plt
import pandas as pd 

def get_data_from_mcc(mcc_path, norm_dose=True):
    """
        This Routine is specifically designed to read the dose meassured using the PTW-BEAMSCAN system in the meassuring protocol
        PDD and Profiles. The Position and meassured Dose of each Meassurement is stored in form of a Pandas DataFrame and all
        of the Dataframes are gathered in a dictionary. 
        In Case two meassurements are part of the .mcc the code assumes those are the profiles and assigns the keys 'x' and 'y'. 
        If only one exists its supposed to be a z-scan (PDD).
        If three exists the order is equal to z -> x? -> y?
        Id more than three exist the keys are increasing integers.         
        
        Input:   mcc_path - <string> To the .mcc-File
        Output:  A dictionary where each entry is corresponding to one scan. 
    """
    #--- test how many messurements are in the file
    with open(mcc_path) as file:
        count = file.readlines().count("\t\tBEGIN_DATA\n")
        if count in [3, 1]: # if count ==3 or count == 1
            axis="z"
        elif count==2:
            axis="y"
        else:
            axis = 0
    
    #--- iterate line by line and extract data
    with open(mcc_path) as file:
        dose_df, temp_dose_df, = {}, {"position":[], "dose":[]}
        line = file.readline()
        while line:
            if line == "\t\tBEGIN_DATA\n":
                line = file.readline() #skip this line to start with actual data
                while line:
                    #--- if the meassurement ends store the data and refresh variables & set axis
                    if line == "\t\tEND_DATA\n":
                        dose_df[axis] = pd.DataFrame(temp_dose_df) 
                        #--- normalize dose if input says so, default is True
                        if norm_dose:
                            dose_df[axis].dose = dose_df[axis].dose/max(dose_df[axis].dose)
                        #--- reset temporary storage
                        temp_dose_df = {"position":[], "dose":[]}
                        if axis=="z":
                            axis="y"
                        elif axis=="y":
                                axis="x"
                        elif axis=="x":
                            pass
                        else:
                            axis += 1 
                        break
                    #--- if not extract float values from one line of the mcc file (order: position, dose, ?)
                    else:
                        output, temp = [], ""
                        for char in line:
                            if char not in  ["\t", "\n"]:
                                temp += char
                            else:
                                #--- only append if temp is not empty
                                if temp:
                                    output.append(float(temp))
                                    temp = ""
                        temp_dose_df["position"].append(output[0])
                        temp_dose_df["dose"].append(output[1])
                        line = file.readline()
            #--- continue scanning through the file for more meassurements        
            line = file.readline()
        return dose_df
        
# %%
test1="C:/Users/apel04/Desktop/master/Data for Python/mcc/6MeV_10x10_PDD_hori_alignement.mcc"
test2="C:/Users/apel04/Desktop/master/Data for Python/mcc/6MeV_10x10_Profile_vert_alignement.mcc"
test3="C:/Users/apel04/Desktop/master/Data for Python/mcc/deleteme.mcc"

#--- read data
data = get_data_from_mcc(test3, norm_dose=False)
def plot_profiles_and_difference(arr_A, arr_B, xvals=range(-80,81,1), norm_dose=True):
    if norm_dose:
        diff = (arr_A-arr_B)/max(max(arr_A), max(arr_B))*100
        ylab0 = "relative absorbed dose to water [%]"
        ylab1 = "Dosedifference between both orientations [%]"
        
    else:
        diff = arr_A - arr_B
        ylab0 = "absolute absorbed dose to water [cGy]"
        ylab1 = "Dosedifference between both orientations [cGy]"
        
    #--- plot data
    fig, axs = plt.subplots(2,1, figsize=(14,7), sharex=True)
    axs[0].plot(xvals, arr_A, color="blue", ls="dashed", lw=2, label="X-Profile")
    axs[0].plot(xvals, arr_B, color="red" , ls="dashed", lw=2, label="Y-Profile")
    axs[1].plot(xvals, diff, color="black", lw=3)
    #--- set axis labels
    axs[1].set_xlabel("lateral position perpendicular to central axis[cm]", size=12)
    axs[0].set_ylabel(ylab0, size=10)
    axs[1].set_ylabel(ylab1, size=10)
    #---set limits
    axs[0].set_xlim(-80,80)
    #---set custom ticks
    axs[0].set_xticks(range(-80,90,10))
    #axs[0].set_yticks(range(-80,90,10))
    #axs[1].set_yticks(range(-80,90,10))
    #--- add grid
    axs[0].grid(True)
    axs[1].grid(True)
    #--- add label
    axs[0].legend(loc="lower center", fontsize=20)
    #--- squeeze plots together
    fig.tight_layout()
    #--- allign ylabels
    fig.align_ylabels(axs)
    

plot_profiles_and_difference(data["y"].dose, data["x"].dose, norm_dose=False)
plot_profiles_and_difference(data["y"].dose, data["x"].dose, norm_dose=True)
     




























































