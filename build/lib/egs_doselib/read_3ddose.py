"""
@author: Apelova 
https://github.com/Apelova/EGS_DOSE_TOOLS

The function add_profile is a unique feature I do not recommend you to use.
This is a routine I added for myself which allows me to use multiple 3ddose files from egs++.
Each of them only scores the profile at one z and therfore saves storage/computation time.
"""
import pandas as pd
import numpy as np
import os


class xyz_array:
    def __init__(self, array):
        self.x = np.array(array["x"])
        self.y = np.array(array["y"])
        self.z = np.array(array["z"])


def norm(arr):
    return np.array(arr/max(arr)*100)

class dose_3d:
    """
        A class that enables one to read .3ddose-Files with python and read certain profiles or arbitrary voxels.
        For more Info see Github 'https://github.com/Apelova/EGS_DOSE_TOOLS'.
    """
    def __init__(self, PATH, INFO=False):
        self.origin = PATH
        if self.__test_path(INFO):        
            #--- boolean value that decides how dose_3d instance is treated/printed
            self.has_multiple_inputs = False
            #--- load data
            self.__read_3ddose_file() 
            #--- set voxel position to voxelcenter
            self.__replace_boundary_with_voxel_center() 
            #--- change datastructure
            self.__dicts_to_xyz_arrays()
            #--- on read initialize pdd @(0,0) and profiles along x and y at z=5cm
            self.__set_profiles()
            #--- gather information for read
            self.__get_profile_depths()
            if INFO:
            #--- output information on read
                print(self)
        else:
            raise TypeError("The Input PATH is invalid!" )

    def __str__(self):
        """Print Information about the 3ddose object """
        if self.has_multiple_inputs:
            return f"""
##############################################################################
  Successfully read the 3ddose-Files and initialized the default PDD and profiles. \n
  X-Profile is available within the Boundaries below for Z {self.available_depths_x.tolist()}:
   {self.boundaries.x.tolist()}\n

  Y-Profile is available within the Boundaries below for Z {self.available_depths_y.tolist()}:
   {self.boundaries.y.tolist()}\n

  Z-Profile (PDD) is available within the Boundaries:
   {self.boundaries.z.tolist()}\n
   
  NOTE THAT IT IS ONLY POSSIBLE TO SELECT VALUES ON THOSE PROFILES SINCE THEY ARE READ FROM MULTIPLE 3ddose-files         

  You can access the Dose-values quickly through OBJECT_NAME.pdd, OBJECT_NAME.x_profile, OBJECT_NAME.y_profile.
  Alternatnatively you can use set_profiles or set_pdd to change the position of each Doseprofile or use get_plane to get a 2D-np.array representing a plane.
  If you want to access individual Elements use the class attribute  'dose_matrix' an access it with dose_matrix[z_index, y_index, x_index].
  
  To save a profile one can just call the pandas Datafram class method pd.DataFrame(profile).to_csv()  e.g pd.DataFrame(dose_3d_object.y_profile).to_csv(destination_path)      
##############################################################################
""" 
        else:
            return f"""
##############################################################################
  Successfully read the 3ddose-File and initialized the default PDD and profiles. \n
  Voxel Boundaries in X-direction:
  {self.boundaries.x.tolist()}\n

  Voxel Boundaries in Y-direction:
   {self.boundaries.y.tolist()}\n

  Voxel Boundaries in Z-direction:
   {self.boundaries.z.tolist()}\n
            
  You can access the Dose-values quickly through OBJECT_NAME.pdd, OBJECT_NAME.x_profile, OBJECT_NAME.y_profile.
  Alternatnatively you can use set_profiles or set_pdd to change the position of each Doseprofile or use get_plane to get a 2D-np.array representing a plane.
  If you want to access individual Elements use the class attribute  'dose_matrix' an access it with dose_matrix[z_index, y_index, x_index].
  
  
  To save a profile one can just call the pandas Datafram class method pd.DataFrame(profile).to_csv()  e.g pd.DataFrame(dose_3d_object.y_profile).to_csv(destination_path)      
##############################################################################
"""     
    def __test_path(self, INFO_INPUT):
        if type(self.origin) != str:
            raise TypeError("The Input for Parameter PATH is not at string !")
        
        if type(INFO_INPUT) != bool:
            raise TypeError("The Input for Parameter INFO is not at boolean value!")
            
        return os.path.isfile(self.origin) and self.origin.endswith(".3ddose")

    #__ makes this function inaccessible to the user -> sets it as private
    def __read_3ddose_file(self):        
        #--- initialize data structures
        self.voxel_in_axis = {"x":[], "y":[], "z":[]}
        self.boundaries = {"x":[], "y":[], "z":[]}
        
        #---iterate through 3ddose file 
        with open(self.origin) as file:
            #--- set voxel count
            self.voxel_in_axis["x"], self.voxel_in_axis["y"], self.voxel_in_axis["z"] = [int(x) for x in file.readline().split(" ") if x!=""]
            #--- read voxel boundaries
            for axis in self.boundaries:
                self.boundaries[axis] =[float(x) for x in file.readline().split(" ") if x not in ["", "\n"] ]
            #--- read dose and error 
            self.dose_array = np.array([float(x) for x in file.readline().split()])
            self.dose_matrix = self.__matrix_from_array(self.dose_array)
            self.error_matrix =  self.__matrix_from_array(np.array([float(x) for x in file.readline().split()]))

    def __replace_boundary_with_voxel_center(self):
        self.position = {axis: [sum(self.boundaries[axis][i:i+2])/2 for i in range(len(self.boundaries[axis])-1)] for axis in self.boundaries}
        
    def __matrix_from_array(self,data_array_3d):
        nx, ny, nz = list(self.voxel_in_axis.values())
        matrix_array = []

        for iz in range(nz):
            #--- cut data arrays into arrays that include the valuues for a plane
            plane_array = data_array_3d[iz*nx*ny:(iz+1)*nx*ny]
            #--- add the data row by row into the plane
            for iy in range(ny):
                if iy==0:
                    plane = np.array([plane_array[0:nx]])
                else:
                    plane = np.vstack([plane, plane_array[iy*nx: (iy+1)*nx]])
            #--- add plane to an array that stores the planes as elements of the array
            matrix_array.append(plane)
        #--- stack the planes to a matrix
        matrix = np.dstack([matrix_array])
        del matrix_array, plane_array, plane #in case garbage collector doesnt 
        return matrix
  
    def __dicts_to_xyz_arrays(self):
        self.voxel_in_axis = xyz_array(self.voxel_in_axis)
        self.boundaries = xyz_array(self.boundaries)
        self.position = xyz_array(self.position)

    def find_closest_index(self, array, value):
        """ used to find the index in the array that is closest to the value """
        idx = abs(array - value).tolist().index( min(abs(array - value)) )
        if np.sum( abs(array-value) == min(abs(array - value)))> 1:
            print(f"WARNING --- multiple voxels in same distance to input {value}. Selected Voxel with center at {array[idx]}")
        return idx 
        
    #--- setter functions
    def __set_profiles(self):
        """ set default pdd and profiles on read """
        self.set_pdd(MUTE=True)
        self.set_x_profile(MUTE=True)
        self.set_y_profile(MUTE=True)
                
    def set_pdd(self, X=0, Y=0, MUTE=False):
        if type(X) and type(Y) in [float, int]:
            if self.has_multiple_inputs:
                raise NotImplementedError("The Instance has multiple .3ddose Inputs and thus this function is unavailable. To get the profile adde use .pdd!")
        
            x_index = self.find_closest_index(self.position.x, X)
            y_index = self.find_closest_index(self.position.y, Y)
            self.pdd = self.dose_matrix[:, y_index, x_index]
            self.pdd_error = self.error_matrix[:, y_index, x_index]
            if not MUTE:
                print(f"Loaded Z-Profile/PDD perpendicular to XY-plane going through the Voxel at (X,Y)=({self.position.x[x_index]},{self.position.y[y_index]})!\n")
        else:
            raise TypeError("The Parameters X and Y accept only integer and Float values !")
                
            
    def set_x_profile(self, Z=0, Y=0, MUTE=False):
        """ default depth for doseprofile is 5cm"""
        if type(Z) and type(Y) in [float, int]:
            if self.has_multiple_inputs:
                raise NotImplementedError("The Instance has multiple .3ddose Inputs and thus this function is unavailable. To get the profile adde use .x_profile!")
        
            y_index = self.find_closest_index(self.position.y, Y)
            z_index = self.find_closest_index(self.position.z, Z)
            self.x_profile = self.dose_matrix[z_index, y_index, :]
            self.x_profile_error = self.error_matrix[z_index, y_index, :]
            self.profile_depth_x = self.position.z[z_index]
            if not MUTE:
                print(f"Loaded X-Profile at the depth Z={self.position.z[z_index]} going through Y={self.position.y[y_index]}!\n")
        else:
            raise TypeError("The Parameters Z and Y accept only integer and Float values !")
        
        
    def set_y_profile(self, Z=0, X=0, MUTE=False):
        """ default depth for doseprofile is 5cm"""
        if type(Z) and type(X) in [float, int]:
            if self.has_multiple_inputs:
                raise NotImplementedError("The Instance has multiple .3ddose Inputs and thus this function is unavailable. To get the profile adde use .y_profile!")
            
            x_index = self.find_closest_index(self.position.x, X)
            z_index = self.find_closest_index(self.position.z, Z)
            self.y_profile = self.dose_matrix[z_index, :, x_index]
            self.y_profile_error = self.error_matrix[z_index, :, x_index]
            self.profile_depth_y = self.position.z[z_index]
            if not MUTE:
                print(f"Loaded Y-Profile at the depth Z={self.position.z[z_index]} going through X={self.position.x[x_index]}!\n")
        else:
            raise TypeError("The Parameters Z and Y accept only integer and Float values !")
        
    #--- getter functions
    def get_pdd(self):
        """ return the currently set Z-Profile/(PDD) and the corresponding statistical error"""
        return self.pdd, self.pdd_error
    
    def get_x_profile(self):
        """ return the currently set X-Profile and the corresponding statistical error"""
        return self.x_profile, self.x_profile_error
        
    def get_y_profile(self):
        """ return the currently set Y-Profile and the corresponding statistical error"""
        return self.y_profile, self.y_profile_error
             
    def get_plane(self, AXIS, POSITION_ON_AXIS):  
        """Returns a plane perpendicular to AXIS at the POSITION_ON_AXIS as 2d-array that."""
        if AXIS.lower()=="x":
           if not self.__value_is_in_boundaries(AXIS, POSITION_ON_AXIS):
               raise ValueError(f"Invalid input for Parameter POSITION_ON_AXIS: {POSITION_ON_AXIS} is not in the boundaries [{self.boundaries.x[0]}, {self.boundaries.x[-1]}]")
           X_INDEX = self.find_closest_index(self.position.x, POSITION_ON_AXIS)    
           return self.dose_matrix[:, :, X_INDEX]
        
        elif AXIS.lower()=="y":
           if not self.__value_is_in_boundaries(AXIS, POSITION_ON_AXIS):
               raise ValueError(f"Invalid input for Parameter POSITION_ON_AXIS: {POSITION_ON_AXIS} is not in the boundaries [{self.boundaries.y[0]}, {self.boundaries.y[-1]}]")
           
           Y_INDEX = self.find_closest_index(self.position.y, POSITION_ON_AXIS)    
           return self.dose_matrix[:, Y_INDEX, :]
        
        elif AXIS.lower()=="z":
           if not self.__value_is_in_boundaries(AXIS, POSITION_ON_AXIS):
               raise ValueError(f"Invalid input for Parameter POSITION_ON_AXIS: {POSITION_ON_AXIS} is not in the boundaries [{self.boundaries.z[0]}, {self.boundaries.z[-1]}]")
            
           Z_INDEX = self.find_closest_index(self.position.z, POSITION_ON_AXIS)    
           return self.dose_matrix[Z_INDEX, :, :]
        
        else:
            raise TypeError(f"Invalid input for Parameter AXIS: '{AXIS}' is not an option viable options are x or X, y or Y, z or Z")
        
    #--- assisting functions
    def __value_is_in_boundaries(self, AXIS, POSITION_ON_AXIS):
        if AXIS.lower()=="x":
           return self.boundaries.x[0] < POSITION_ON_AXIS < self.boundaries.x[-1]
        
        elif AXIS.lower()=="y":
           return self.boundaries.y[0] < POSITION_ON_AXIS < self.boundaries.y[-1]
        
        elif AXIS.lower()=="z":
           return self.boundaries.z[0] < POSITION_ON_AXIS < self.boundaries.z[-1]
    
    def add_profile(self, PATH, AXIS):
        """If the initial 3ddose-file is just a PDD you can add profiles with add_profile!"""        
        #--- read the additional 3ddose file and assign its profiles accordingly
        profile = dose_3d(PATH)
        if AXIS.lower()=="x":
            if len(profile.position.x) <= 10:
                print(f"Warning---The X-profile you want to add has only {len(profile.position.x)} Entries. I hope you know what you are doing.")
                #--- set the initial positions for x and y such that adding profiles wont allow setting faulty z-profiles (when adding the first profile only)
            self.has_multiple_inputs = True
            self.boundaries.x = profile.boundaries.x
            self.position.x = profile.position.x
            self.x_profile = profile.x_profile
            self.x_profile_error = profile.x_profile_error
            if np.array_equal(self.available_depths_x, self.position.z):
                self.available_depths_x = np.append([], profile.position.z)
            else:
                self.available_depths_x = np.append(self.available_depths_x, profile.position.z)

            print(f"added x-profile at z {profile.position.z.tolist()}!\n")
            
        elif AXIS.lower()=="y":
            if len(profile.position.y) <= 10:
                print(f"Warning---The Y-profile you want to add has only {len(profile.position.y)} Entries. I hope you know what you are doing.")
            #--- set the initial positions for x and y such that adding profiles wont allow setting faulty z-profiles (when adding the first profile only)
            self.has_multiple_inputs = True
            self.boundaries.y = profile.boundaries.y
            self.position.y = profile.position.y
            self.y_profile = profile.y_profile
            self.y_profile_error = profile.y_profile_error
            if np.array_equal(self.available_depths_y, self.position.z): #different definition between single and lonely pdd
                self.available_depths_y = np.append([], profile.position.z)
            else:
                self.available_depths_y = np.append(self.available_depths_y, profile.position.z)
            print(f"added y-profile at z {profile.position.z.tolist()}!\n")
            
        else:
            raise ValueError(f"Invalid input for Parameter AXIS: '{AXIS}' is not an option viable options are x or X, y or Y!\n")

        del profile
        
    def __get_profile_depths(self):
        self.available_depths_x = self.position.z
        self.available_depths_y = self.position.z

    def set_profiles(self, Z, MUTE=False):
        value_in_x = self.available_depths_x[self.find_closest_index(self.available_depths_x, Z)]        
        if value_in_x in self.available_depths_y:
            self.set_x_profile(Z, Y=0, MUTE=MUTE)
            self.set_y_profile(Z, X=0, MUTE=MUTE)
        else:
            raise ValueError("Depth {Z} is not available for both profiles. For more information use .available_depths_x/y !")
            