"""
@author: Apelova 

Module Task: get Data from .mcc-File as output from Mephysto.

Important
- needs implementation to figure out xyz axis ! from grantry positon and if in plane or crossplane !

TODO
- is it possible to get information about the z-depth for the profiles ? ja !!!
SCAN_CURVETYPE=INPLANE_PROFILE !
SCAN_DEPTH=52.50 !

The order z-> x -> y is based on the scanning procedure inside the beamscan phantom !
"""
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

class _xyz_array_Object:
    def __init__(self):
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([])

class dose_mcc:
    """ doc ... """
    def __init__(self, PATH, INFO=False):
        self.origin = PATH
        self.position = _xyz_array_Object()
        self.__load_data()
        #--- set voxel position to voxelcenter
        self.__restructure_data()
        #--- change scale from mm to cm
        self.position.x = self.position.x/10
        self.position.y = self.position.y/10
        self.position.z = self.position.z/10
         
    def __load_data(self):
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
        with open(self.origin) as file:
            self.number_of_scans = file.readlines().count("\t\tBEGIN_DATA\n")
            
            if self.number_of_scans in [3, 1]: # if count ==3 or count == 1
                self.scan_axis ="z"
                
            elif self.number_of_scans==2:
                self.scan_axis ="y"
            
            else:
                self.scan_axis = 0

        with open(self.origin) as file:
            self.dose_df, temp_dose_df, = {}, {"position":[], "dose":[]}
            line = file.readline()

        #--- iterate line by line and extract data
            while line:
                if line == "\t\tBEGIN_DATA\n":
                    line = file.readline() #skip this line to start with actual data
                    while line:
                        #--- if the meassurement ends store the data and refresh variables & set axis
                        if line == "\t\tEND_DATA\n":
                            self.dose_df[self.scan_axis] = pd.DataFrame(temp_dose_df) 
                            #--- reset temporary storage
                            temp_dose_df = {"position":[], "dose":[]}
                            if self.scan_axis=="z":
                                self.scan_axis="y"
                            elif self.scan_axis=="y":
                                    self.scan_axis="x"
                            elif self.scan_axis=="x":
                                pass
                            else:
                                self.scan_axis += 1 
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
            

    def __restructure_data(self):
        print(f"reading Data from {self.origin}...")
        if list(self.dose_df.keys()) == ["z", "y", "x"]:
            self.__load_x()
            self.__load_y()
            self.__load_z()
        
        elif list(self.dose_df.keys()) == ["y", "x"]:
            self.__load_x()
            self.__load_y()
        
        elif list(self.dose_df.keys()) == ["z"]:
            self.__load_z()
            
        else:
            print("Unusual amount of meassurements. To access the data use self.dose_df !")
        
        
    def __load_x(self):
        self.position.x = np.array(self.dose_df["x"].position)
        self.x_profile = np.array(self.dose_df["x"].dose)
        print("Successfully loaded x-Profile !")
    
    def __load_y(self):
        self.position.y = np.array(self.dose_df["y"].position)
        self.y_profile = np.array(self.dose_df["y"].dose)
        print("Successfully loaded y-Profile !")
    
    def __load_z(self):
        self.position.z = np.array(self.dose_df["z"].position)
        self.pdd = np.array(self.dose_df["z"].dose)
        print("Successfully loaded percentage depth dose (PDD) !")



# =============================================================================
# test1="C:/Users/apel04/Desktop/master/Data_for_python/mcc/6MeV_10x10_PDD_hori_alignement.mcc"
# test2="C:/Users/apel04/Desktop/master/Data_for_python/mcc/6MeV_10x10_Profile_vert_alignement.mcc"
# test3="C:/Users/apel04/Desktop/master/Data_for_python/mcc/delete_me_all_axis.mcc"
# dose = dose_mcc(test3)
# 
# # x profile
# plt.plot(dose.position.x, dose.x_profile)
# plt.show()
# 
# # y profile
# plt.plot(dose.position.y, dose.y_profile)
# plt.show()
# 
# # pdd
# plt.plot(dose.position.z, dose.pdd)
# plt.show()
# =============================================================================


















































