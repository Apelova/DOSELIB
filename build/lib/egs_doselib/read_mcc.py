"""
@author: Apelova 
https://github.com/Apelova/EGS_DOSE_TOOLS

Definition of dose_mcc-class to get Data from .mcc-File as output from Mephysto.

So far only Supports one PDD pro dose_mcc-instance! 
"""
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

class _xyz_array_Object:
    def __init__(self):
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([])
        
class _mcc_scan:
    def __init__(self, GANTRY_angle, SCAN_curvetype, SCAN_angle, SCAN_depth, DATA):
        self.GANTRY_angle=GANTRY_angle 
        self.SCAN_curvetype=SCAN_curvetype 
        self.SCAN_angle=SCAN_angle
        self.SCAN_DEPTH = SCAN_depth
        self.DATA=DATA
        self.AXIS=None
        self.POSITION=None
        self.DOSE=None
        self.Error=None
        
class dose_mcc:
    """
    """
    def __init__(self, PATH, INFO=False):
        self.origin = PATH
        #--- extract data from scans and store the raw information in a dictionary .all_scans
        self.__load_scans()
        #--- scan to dose object
        self.__set_scan_parameter()
        #--- set default profiles
        #--- load available profiles PDD, X/Y-Profiles
        self.__load_profiles()
    
        if INFO:
            print(self)
        
        
         
    def __load_scans(self):
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
        self.all_scans = {}
        #extract the raw data and store it as _mcc_scan object in a dictionary     
        with open(self.origin) as file:
            current_scan = 1
        #--- iterate line by line and extract data
            line = file.readline()    
            while line:
                if line == f'\tBEGIN_SCAN  {current_scan}\n':
                    in_scan = True
                    while in_scan:
                        #--- extract data data
                        if line == f"\tEND_SCAN  {current_scan}\n":
                            #next scan if no begin is found
                            current_scan += 1 
                            break
                        elif "GANTRY=" in line:
                            GANTRY_ANGLE=line.split("=")[1][:-1]

                        elif "SCAN_CURVETYPE" in line:
                            SCAN_CURVETYPE=line.split("=")[1][:-1]

                        elif "SCAN_ANGLE" in line:
                            SCAN_ANGLE=line.split("=")[1][:-1]
                        
                        elif "SCAN_DEPTH" in line:
                            SCAN_DEPTH=line.split("=")[1][:-1]
                            
                        #--- extract meassured values
                        elif line == "\t\tBEGIN_DATA\n":
                            line = file.readline()#skip BEGIN_DATA
                            SCAN_DATA = []
                            while True:
                                if line == "\t\tEND_DATA\n":
                                    try:
                                        self.all_scans[current_scan] = _mcc_scan(GANTRY_ANGLE, SCAN_CURVETYPE, SCAN_ANGLE, SCAN_DEPTH, SCAN_DATA)
                                    except:
                                        self.all_scans[current_scan] = _mcc_scan(GANTRY_ANGLE, SCAN_CURVETYPE, SCAN_ANGLE, "unavailable", SCAN_DATA)
                                    current_scan += 1 # next scan after getting data
                                    in_scan = False
                                    break
                               
                                SCAN_DATA.append([float(i) for i in  line.split("\t\t") if i != ""])        
                                #--- for while true    
                                line = file.readline()
                        #for while in_scan
                        line = file.readline()
                #for while line
                line = file.readline()
    
    
    def __set_scan_parameter(self):
        """set axis, position on axis, dose and error
           convert position from mm to cm
        """
        for i in self.all_scans:
            #--- assign xyz label
            if self.all_scans[i].SCAN_curvetype == "PDD":
                self.all_scans[i].AXIS="Z"
            else:    
                self.all_scans[i].AXIS = self.__get_axis_label(i)
            #---convert
            self.all_scans[i].POSITION = np.array(self.all_scans[i].DATA)[:,0]/10 #/10 mm to cm 
            self.all_scans[i].DOSE = np.array(self.all_scans[i].DATA)[:,1]
            try:
                self.all_scans[i].ERROR= np.array(self.all_scans[i].DATA)[:,2]
            except: #some meassurements dont store the error !
                self.all_scans[i].ERROR= np.array([0 for i in self.all_scans[i].DOSE])
            #--- invert data such that coordinate systems match !
            if "-" in self.all_scans[i].AXIS:
                self.all_scans[i].POSITION = self.all_scans[i].POSITION[::-1]*-1 #necessary if scan is not symmetric !
                self.all_scans[i].DOSE = self.all_scans[i].DOSE[::-1]
                self.all_scans[i].ERROR = self.all_scans[i].ERROR[::-1]
                #after fixing coordinate set axis label -X/-Y/-Z -> X/Y/Z
                self.all_scans[i].AXIS = self.all_scans[i].AXIS[1::]
                
    def __get_axis_label(self, i):
        """ Implementation only viabe for GANTRY_ANGLE==0, the axis definition is made such that x from read_mcc== x from read_3ddose"""        
        if self.all_scans[i].GANTRY_angle=="0.00":
            if self.all_scans[i].SCAN_curvetype == "CROSSPLANE_PROFILE":
                return "Y"
            
            elif self.all_scans[i].SCAN_curvetype == "INPLANE_PROFILE":
                return "-X"
    
            else:
                return "undefined"
            
        else:
            return "undefined"
        
    def __load_profiles(self):
        self.pdd, self.x_profile, self.y_profile= None, None, None
        self.position = _xyz_array_Object()
        self.available_depths_x, self.available_depths_y = [], []
        
        for i in self.all_scans:
            action = "loaded"                

            if self.all_scans[i].AXIS == "Z" and self.pdd==None:
                self.position.z = self.all_scans[i].POSITION
                self.pdd = self.all_scans[i].DOSE
                self.pdd_error = self.all_scans[i].ERROR
                print("Successfully set Percentage Depth Dose (PDD)!")
            
            elif self.all_scans[i].AXIS == "X":
                if type(self.x_profile)==type(None):
                    self.position.x = self.all_scans[i].POSITION
                    self.x_profile = self.all_scans[i].DOSE
                    self.x_profile_error = self.all_scans[i].ERROR
                    action = "set"
                self.available_depths_x.append((f"Number of Scan {i}",float(self.all_scans[i].SCAN_DEPTH)/10))
                print(f"Successfully {action} X-Profile at depth Z={float(self.all_scans[i].SCAN_DEPTH)/10}cm!")
            
            elif self.all_scans[i].AXIS == "Y":
                if type(self.y_profile)==type(None):
                    self.position.y = self.all_scans[i].POSITION
                    self.y_profile = self.all_scans[i].DOSE
                    self.y_profile_error = self.all_scans[i].ERROR
                    action = "set"
                self.available_depths_y.append((f"Number of Scan {i}",float(self.all_scans[i].SCAN_DEPTH)/10))
                print(f"Successfully {action} Y-Profile at depth Z={float(self.all_scans[i].SCAN_DEPTH)/10}cm!")
            
            else:
                print(f"Warning---Axis of Scan {i} is unknown!")
        
    def set_x_profile(self, NUMSCAN, MUTE=True):
        try:
            int(NUMSCAN)
        except:
            raise ValueError(f"Input vor Variable NUMSCAN '{NUMSCAN}' must be an integer value but is not.")
        
        if self.all_scans[NUMSCAN].AXIS != "X":
            raise ValueError(f"Scan Number {NUMSCAN} does not represent an X-Scan !!! but a {self.all_scans[NUMSCAN].AXIS}-Scan!")
        self.position.x = self.all_scans[NUMSCAN].POSITION
        self.x_profile = self.all_scans[NUMSCAN].DOSE
        self.x_profile_error = self.all_scans[NUMSCAN].ERROR
        if not MUTE:
            print(f"Successfully set X-Profile at depth Z={float(self.all_scans[NUMSCAN].SCAN_DEPTH)/10}cm!")

    
    def set_y_profile(self, NUMSCAN, MUTE=True):
        try:
            int(NUMSCAN)
        except:
            raise ValueError(f"Input vor Variable NUMSCAN '{NUMSCAN}' must be an integer value but is not.")
        
        if self.all_scans[NUMSCAN].AXIS != "Y":
            raise ValueError(f"Scan Number {NUMSCAN} does not represent an Y-Scan !!! but a {self.all_scans[NUMSCAN].AXIS}-Scan!")
        self.position.y = self.all_scans[NUMSCAN].POSITION
        self.y_profile = self.all_scans[NUMSCAN].DOSE
        self.y_profile_error = self.all_scans[NUMSCAN].ERROR
        if not MUTE:
            print(f"Successfully set Y-Profile at depth Z={float(self.all_scans[NUMSCAN].SCAN_DEPTH)/10}cm!")

    def set_profiles(self, DEPTH, MUTE=True):
        try:
            int(DEPTH)
            float(DEPTH)
        except:
            raise ValueError(f"Input vor Variable DEPTH '{DEPTH}' must be an integer or float value but is not.")
        
        #get index if DEPTH is available at axisd and then the corresponding NUMSCAN than turn it into int
        X_NUMSCAN = int(self.available_depths_x[[i[1] for i in self.available_depths_x].index(DEPTH)][0].split("Scan")[-1])
        Y_NUMSCAN = int(self.available_depths_y[[i[1] for i in self.available_depths_y].index(DEPTH)][0].split("Scan")[-1])
        self.set_x_profile(X_NUMSCAN, MUTE)
        self.set_y_profile(Y_NUMSCAN, MUTE)
        
    #--- show information about the project
    def __str__(self):
        return f""" ##############################################################################
 Successfully read the .mcc-File and initialized the default PDD and profiles. \n
   - Percentage Depth Dose loaded for values of z in [{self.position.z[0]}, {self.position.z[-1]}]
   - Available x-Profile depths are {[ i [1] for i in self.available_depths_x]} at {[ i[0] for i in self.available_depths_x]} 
   - Available y-Profile depths are {[ i [1] for i in self.available_depths_y]} at {[ i[0] for i in self.available_depths_y]}\n
 By default the first available Profiles are set as .x_profile or .y_profile!
 To Change the currently set x/y_profile use .set_x/y_profile(Number of Scan).   
 
 The Data of each SCAN is stored together with some Information about the Scan in a dictionary .all_scans. 
 The Keys if this dictionary are the number of Scans as integer values.
 
 To save a profile one can just use the Pandas Datafram class e.g pd.DataFrame(dose_3d_object.y_profile).to_csv(destination_path)      
##############################################################################"""  