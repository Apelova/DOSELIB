# -*- coding: utf-8 -*-
"""
Spyder Editor
Aufgaben:
    1. lesen der 3ddose files die aus chamber kommen 
    2. lesen der allgemeinen 3ddose files
    3. zwei funktionen schreiben
        - get_1d_dose
        - get_3_ddose
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

x_path = "Desktop/Master/Data_for_Python/3ddose/EX16MeVp_dose_x_profile.3ddose"
y_path = "Desktop/Master/Data_for_Python/3ddose/EX16MeVp_dose_y_profile.3ddose"
z_path = "Desktop/Master/Data_for_Python/3ddose/EX16MeVp_dose_pdd.3ddose" 


    
#dose_1d = read_1d_dose(x_path)

    
# 3ddimensionale dosis classe schreiben

path = "C:/Users/apel04/Desktop/master/comparison_dose_chamber/NRC_dosxyz.3ddose"

class xyz_array_from_dict:
    def __init__(self, array):
        self.x = np.array(array["x"])
        self.y = np.array(array["y"])
        self.z = np.array(array["z"])

class dose3d:
    """
    """
    def __init__(self, path):
        self.origin = path
        self.__read_3ddose_file() #--- load data
        self.__replace_boundary_with_voxel_center() #--- set voxel position to voxelcenter
        self.dicts_to_xyz_arrays()

        
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
            self.dose_matrix = self.__matrix_from_array(np.array([float(x) for x in file.readline().split()]))
            self.error_matrix =  self.__matrix_from_array(np.array([float(x) for x in file.readline().split()]))

    def __replace_boundary_with_voxel_center(self):
        self.dose_array = np.array([])
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
  
    def dicts_to_xyz_arrays(self):
        self.voxel_in_axis = xyz_array_from_dict(self.voxel_in_axis)
        self.boundaries = xyz_array_from_dict(self.boundaries)
        self.position = xyz_array_from_dict(self.position)

    def find_closest_index(self, array, value):
        """assumes the array is ordered !!!"""
        return abs(array - value).tolist().index( min(abs(array - value)))
          
    #TODO!
    def get_pdd(self, X=0.1, Y=0):
        x_index = self.find_closest_index(self.position.x, X)
        #y_index = self.find_closest_index(self.position.y, Y)
        #print(X, "  found at  ", x_index, " because:")
        #for i in range(len(self.position.x)):
        #    print(i, self.position.x[i])
        
        pass
        
    #TODO!
    def get_profile(self):
        pass

    





        

dose = dose3d(path)

dose.get_pdd()














































