"""
@author: Apelova 

Module Task: get Data from .agr-File as output from statdose.

Important
- make a dose_3d from ptroject name
- therefore read all 3 profiles from chamber 3ddose output     

TODO
- 
- add Documentation
    - to save an array as csv use pd.DataFrame(array).to_csv("path") !!
- add console output to show where the voxel lays for a certain input value !!!
- type testing in every function
- add example data to the repo

"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

class xyz_array_from_dict:
    def __init__(self, array):
        self.x = np.array(array["x"])
        self.y = np.array(array["y"])
        self.z = np.array(array["z"])

class dose_3d:
    """
    """
    def __init__(self, PATH, INFO=False):
        self.origin = PATH
        #--- load data
        self.__read_3ddose_file() 
        #--- set voxel position to voxelcenter
        self.__replace_boundary_with_voxel_center() 
        #--- change datastructure
        self.__dicts_to_xyz_arrays()
        #--- on read initialize pdd @(0,0) and profiles along x and y at z=5cm
        self.__set_profiles()
        #--- output information on read
        if INFO:
            print(self)
        #--- boolearn to decide what to print
        self.has_multple_inputs = False
        
    
    def __str__(self, ):
        if self.has_multple_inputs:
            return f"""
##############################################################################
#  Successfully read the 3ddose-Files and initialized the default PDD and profiles. \n#
#  X-Profile is available within the Boundaries:
   {self.boundaries.x.tolist()}
#
#  Y-Profile is available within the Boundaries:
   {self.boundaries.y.tolist()}
#
#  Z-Profile is available within the Boundaries:
   {self.boundaries.z.tolist()}
#   
#  NOTE THAT IT IS ONLY POSSIBLE TO SELECT VALUES ON THOSE PROFILES SINCE THEY ARE READ FROM MULTIPLE 3ddose-files         
#  You can access the Dose-values quickly through OBJECT_NAME.pdd, OBJECT_NAME.x_profile, OBJECT_NAME.y_profile.
#  Alternatnatively you can use set_profiles or set_pdd to change the position of each Doseprofile or use get_plane to get a 2D-np.array representing a plane.
##############################################################################
""" 
            
        else:
            return f"""
##############################################################################
#  Successfully read the 3ddose-File and initialized the default PDD and profiles. \n#
#  Voxel Boundaries in X-direction:
   {self.boundaries.x.tolist()}
#
#  Voxel Boundaries in Y-direction:
   {self.boundaries.y.tolist()}
#
#  Voxel Boundaries in Z-direction:
   {self.boundaries.z.tolist()}
#            
#  You can access the Dose-values quickly through OBJECT_NAME.pdd, OBJECT_NAME.x_profile, OBJECT_NAME.y_profile.
#  Alternatnatively you can use set_profiles or set_pdd to change the position of each Doseprofile or use get_plane to get a 2D-np.array representing a plane.
##############################################################################
"""     
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
        self.voxel_in_axis = xyz_array_from_dict(self.voxel_in_axis)
        self.boundaries = xyz_array_from_dict(self.boundaries)
        self.position = xyz_array_from_dict(self.position)

    def find_closest_index(self, array, value):
        idx = abs(array - value).tolist().index( min(abs(array - value)) )
        if np.sum( abs(array-value) == min(abs(array - value)))> 1:
            print(f"WARNING --- multiple voxels in same distance to input {value}. Selected Voxel with center at {array[idx]}")
        return idx 
          
    #--- setter functions
    def __set_profiles(self):
        self.set_pdd()
        self.set_x_profile()
        self.set_y_profile()
        
    def set_pdd(self, X=0, Y=0):
        x_index = self.find_closest_index(self.position.x, X)
        y_index = self.find_closest_index(self.position.y, Y)
        self.pdd = self.dose_matrix[:, y_index, x_index]
        self.pdd_error = self.error_matrix[:, y_index, x_index]
            
    #TODO! --- test
    def set_x_profile(self, Z=0, Y=0):
        """ default depth for doseprofile is 5cm"""
        y_index = self.find_closest_index(self.position.y, Y)
        z_index = self.find_closest_index(self.position.z, Z)
        self.x_profile = self.dose_matrix[z_index, y_index, :]
        self.x_profile_error = self.error_matrix[z_index, y_index, :]
        
    def set_y_profile(self, Z=0, X=0):
        """ default depth for doseprofile is 5cm"""
        x_index = self.find_closest_index(self.position.x, X)
        z_index = self.find_closest_index(self.position.z, Z)
        self.y_profile = self.dose_matrix[z_index, :, x_index]
        self.y_profile_error = self.error_matrix[z_index, :, x_index]
    
    #--- getter functions
    def get_pdd(self):
        return self.pdd, self.pdd_error
    
    def get_x_profile(self):
        return self.x_profile, self.x_profile_error
        
    def get_y_profile(self):
        return self.y_profile, self.y_profile_error
             
    def get_plane(self, AXIS, POSITION_ON_AXIS):  
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
        #TODO! type testing !
        
        profile = dose_3d(PATH)
        if AXIS.lower()=="x":
            self.has_multple_inputs = True
            self.boundaries.x = profile.boundaries.x
            self.position.x = profile.position.x
            self.x_profile = profile.x_profile
            self.x_profile_error = profile.x_profile_error
            
        elif AXIS.lower()=="y":
            self.has_multple_inputs = True
            self.boundaries.y = profile.boundaries.y
            self.position.y = profile.position.y
            self.y_profile = profile.y_profile
            self.y_profile_error = profile.y_profile_error

        else:
            raise ValueError(f"Invalid input for Parameter AXIS: '{AXIS}' is not an option viable options are x or X, y or Y!")
            
        del profile
#path = "/home/marvin/EGSnrc/EGSnrc/egs_home/dosxyznrc/16MVp_h2o_phantom_beamsource_example.3ddose"#
#dose = dose_3d(path)
#plt.imshow(dose.get_plane("", 10))
































