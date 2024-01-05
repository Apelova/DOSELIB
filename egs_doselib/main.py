"""
@author: Apelova 
https://github.com/Apelova/EGS_DOSE_TOOLS

Definition of all classes and methods used in to analyse .3ddose or .mcc files.


"""
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd 
import numpy as np
import os

class _xyz_array_Object:
    def __init__(self, array):
        self.x = np.array(array["x"])
        self.y = np.array(array["y"])
        self.z = np.array(array["z"])

class normable_array(np.ndarray):
    def __new__(cls, arr):
        obj = np.asarray(arr).view(cls)
        return obj

    def __init__(self, arr):
        self.norm = self/max(self)*100

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

class dose_object:
    def __init__(self): 
        #--- Input Attributes
        self.origin             = None
        self.has_multiple_inputs= None
        #--- Geometry Attributes
        self.voxel_in_axis      = None #_xyz_array_Object({"x":[], "y":[], "z":[]})
        self.boundaries         = None #_xyz_array_Object({"x":[], "y":[], "z":[]})
        self.position           = None
        #--- Unprocessed Dose
        self.dose_array         = None
        self.dose_matrix        = None
        #--- Unprocessed Error
        self.error_array        = None
        self.error_matrix       = None
        #--- Information about current and available X,Y Profiles
        self.available_depths_x = None
        self.available_depths_y = None
        self.profile_depth_x    = None #current depth
        self.profile_depth_y    = None #current depth
        #--- currently loaded profiles
        self.pdd                = None #dose arrays as normable_array([])
        self.x_profile          = None #dose arrays as normable_array([])
        self.y_profile          = None #dose arrays as normable_array([])
        #--- corresponding error to each profile
        self.pdd_error          = None #dose errors in percent (relative to each index)
        self.x_profile_error    = None #dose errors in percent (relative to each index)
        self.y_profile_error    = None #dose errors in percent (relative to each index)
        #---  enables reusability of calculated indices -> inaccessibel for user
        self.__l_ix = None
        self.__r_ix= None
        #--- initialize all possible metrics NaN 
        self.pdd_max_abs = np.nan
        #--- PDD metrics
        self.pdd_i_max  = np.nan
        self.Rmax= np.nan
        self.D0= np.nan
        self.D10= np.nan
        self.D20= np.nan
        self.Q= np.nan
        #--- X-Profile metrics
        self.D0_x= np.nan
        self.plateau_limits_x= np.nan
        self.penumbra_limits_x_20= np.nan
        self.penumbra_limits_x_80= np.nan
        self.penumbra_yl= np.nan
        self.penumbra_yr= np.nan
        self.H_x= np.nan
        self.HW_x= np.nan
        self.S_x= np.nan
        self.dCAX_x= np.nan
        #--- Y-Profile metrics
        self.D0_y= np.nan
        self.plateau_limits_y= np.nan
        self.penumbra_limits_y_20= np.nan
        self.penumbra_limits_y_80= np.nan
        self.penumbra_yl= np.nan
        self.penumbra_yr= np.nan
        self.H_y= np.nan
        self.HW_y= np.nan
        self.S_y= np.nan
        self.dCAX_y= np.nan
        
    def find_closest_index(self, array, value, show_warning=True):
        """ used to find the index in the array that is closest to the value """
        idx = abs(array - value).tolist().index( min(abs(array - value)) )
        if np.sum( abs(array-value) == min(abs(array - value)))> 1 and show_warning:
            print(f"WARNING --- multiple voxels in same distance to input {value}. Selected Voxel with center at {array[idx]}")
        return idx 
    
    def __get_interpol_indices(self, z):
        if z==10:
            l_ix, r_ix = self.i_10, self.i_10

        if z==20:
            l_ix, r_ix = self.i_20, self.i_20

        while self.position.z[l_ix] > z:
            l_ix -= 1
        
        while self.position.z[r_ix] < z:
            r_ix += 1
        
        return l_ix, r_ix
    
    def set_pdd_metrics(self):
        """ 
            The values D0, D10, D20 are all relative and can be transforme to 
            absolute dose values through multiplication with self.pdd_max_absma
        """
        
        #--- maximum absolute Dose
        self.pdd_max_abs = max(self.pdd)
        #--- corresponting index to the maximum
        self.pdd_i_max   = self.pdd.norm.tolist().index(100)
        #--- corresponting Position on z-Axis to the maximum
        self.Rmax        = self.position.z[self.pdd_i_max]

        #--- relative Surfacedose
        if (0 in self.position.z):
            self.D0 = round(self.pdd.norm[self.find_closest_index(self.position.z, 0, show_warning=False)], 3) # den wert über index nehmen !
        else:
            f = interpolate.interp1d(self.position.z[:2], self.pdd.norm[:2], kind="linear", fill_value='extrapolate')
            self.D0 = round(f([0])[0], 3)
        
        #--- indices at or close to depth 10cm/20cm 
        self.i_10 = self.find_closest_index(self.position.z, 10, show_warning=False)
        self.i_20 = self.find_closest_index(self.position.z, 20, show_warning=False)        

        #--- get the Dose at 10 cm if neccessary interpolate (quadratic)!
        if (10 in self.position.z):
            self.D10 = round(self.pdd.norm[self.i_10], 3)
        else: 
            l_ix, r_ix = self.__get_interpol_indices(10)
            f = interpolate.interp1d(self.position.z[[l_ix, r_ix]], self.pdd.norm[[l_ix, r_ix]], kind="linear", fill_value='extrapolate')
            self.D10 = round(f([10])[0], 3)

        #--- get the Dose at 20 cm from the index and if it is equidistant to two points average !
        if (20 in self.position.z):
            self.D20 = round(self.pdd.norm[self.i_20], 3)
        else: 
            l_ix, r_ix = self.__get_interpol_indices(20)
            f = interpolate.interp1d(self.position.z[[l_ix, r_ix]], self.pdd.norm[[l_ix, r_ix]], kind="linear", fill_value='extrapolate')
            self.D20 = round(f([20])[0], 3)

        #--- Quality Index after DIN 6800-2
        self.Q = round(1.2661*self.D20/self.D10 - 0.0595, 3)
    
    
    def set_metrics_profile(self, AXIS):
        self.__set_Dose_on_axis(AXIS)
        self.__set_penumbra_limits_20(AXIS)
        self.__set_penumbra_limits_80(AXIS)
        self.__set_penumbra_width(AXIS)
        self.__set_halfwidth(AXIS)
        self.__set_dCAX(AXIS)
        self.__set_homogeneity(AXIS)
        #self.__set_symmetry(AXIS)
        
    def __set_Dose_on_axis(self, AXIS): #hier muss man eigentlich auch interpolieren !
        """ it is expected that dose the profiles include 0 !!"""
        if AXIS =="X":
            self.D0_x = self.x_profile.norm[self.find_closest_index(self.position.x, 0)]
    
        if AXIS =="Y":
            self.D0_y = self.y_profile.norm[self.find_closest_index(self.position.y, 0)]
            
    def get_limits_from_interpolation(self, position, dose_array, dose_value, degree="linear"):
        """ Used for FWHM, 20%, 80% -> Penumbra widths"""
        l_ix = self.find_closest_index(dose_array, dose_array[dose_array > dose_value][0])
        r_ix = self.find_closest_index(dose_array, dose_array[dose_array > dose_value][-1])
        #--- get linear functions for interpolation
        if degree == "linear":
            f_l = interpolate.interp1d(position[[l_ix-1, l_ix]], dose_array[[l_ix-1, l_ix]], fill_value='extrapolate', kind="linear")
            f_r = interpolate.interp1d(position[[r_ix, r_ix+1]], dose_array[[r_ix, r_ix+1]], fill_value='extrapolate', kind="linear")
            
        elif degree == "quadratic":
            f_l = interpolate.interp1d(position[[l_ix-1, l_ix, l_ix+1]], dose_array[[l_ix-1, l_ix, l_ix+1]], fill_value='extrapolate', kind="quadratic")
            f_r = interpolate.interp1d(position[[r_ix-1, r_ix, r_ix+1]], dose_array[[r_ix-1, r_ix, r_ix+1]], fill_value='extrapolate', kind="quadratic")

        else:
            raise TypeError("Invalid polynomial degree. Valid Inputs are 'linear' or 'quadratic'")
        #--- search for root of the functions below to find x !
        f_l_ = lambda x: f_l(x)-dose_value
        f_r_ = lambda x: f_r(x)-dose_value
        #--- linear interpolation to point x where D=0.8'D0
        x_l = fsolve(f_l_, position[l_ix])                
        x_r = fsolve(f_r_, position[r_ix])
        return np.array([x_l[0], x_r[0]]) #return left and right value

    def __set_halfwidth(self, AXIS):
       if AXIS =="X":
            x_l, x_r = self.get_limits_from_interpolation(self.position.x, self.x_profile.norm, 50, degree = "quadratic")
            #--- store outer limtis                
            self.HW_x = np.array([x_l, x_r])
     
       if AXIS =="Y":
            x_l, x_r = self.get_limits_from_interpolation(self.position.y, self.y_profile.norm, 50, degree = "quadratic")
            #--- store outer limtis                
            self.HW_y = np.array([x_l, x_r])    

    def __set_penumbra_limits_80(self, AXIS): #hier ist wahrscheinlich eine lineare interpolation besser geeignet !
        # Hier die werte für x==80 prozent her holen
       if AXIS =="X":
            x_l, x_r = self.get_limits_from_interpolation(self.position.x, self.x_profile.norm, 80, degree = "quadratic")
            self.penumbra_limits_x_80 = np.array([x_l, x_r])
                
       if AXIS =="Y":
            x_l, x_r = self.get_limits_from_interpolation(self.position.y, self.y_profile.norm, 80, degree = "quadratic")
            self.penumbra_limits_y_80 = np.array([x_l, x_r])
                        
    def __set_penumbra_limits_20(self, AXIS):
       if AXIS =="X":
            x_l, x_r = self.get_limits_from_interpolation(self.position.x, self.x_profile.norm, 20, degree="quadratic")
            self.penumbra_limits_x_20 = np.array([x_l, x_r])
    
       if AXIS =="Y":
            x_l, x_r = self.get_limits_from_interpolation(self.position.y, self.y_profile.norm, 20, degree="quadratic")
            self.penumbra_limits_y_20 = np.array([x_l, x_r])

    def __set_penumbra_width(self, AXIS):
       if AXIS =="X":
            self.penumbra_xl= self.penumbra_limits_x_80[0]-self.penumbra_limits_x_20[0]
            self.penumbra_xr= self.penumbra_limits_x_20[1]-self.penumbra_limits_x_80[1]
       if AXIS =="Y":
            self.penumbra_yl= self.penumbra_limits_y_80[0]-self.penumbra_limits_y_20[0]
            self.penumbra_yr= self.penumbra_limits_y_20[1]-self.penumbra_limits_y_80[1]

    def __set_homogeneity(self, AXIS): 
       if AXIS =="X":
            #--- find index/limits inside 0.8 of the FWHM
            self.__l_ix  = self.find_closest_index(self.position.x, self.position.x[self.HW_x[0]*0.8 < self.position.x][0])
            self.__r_ix  = self.find_closest_index(self.position.x, self.position.x[self.HW_x[1]*0.8 < self.position.x][0])
            #--- caluclate homgeinity defined as (D_max-Dmin)/D0                 
            self.H_x = (max(self.x_profile.norm[self.__l_ix : self.__r_ix+1])-min(self.x_profile.norm[self.__l_ix : self.__r_ix+1]))/self.D0_x 
            self.plateau_limits_x = np.array([self.__l_ix, self.__r_ix])
              
       if AXIS =="Y":
            #--- find index/limits inside 0.8 of the FWHM
            self.__l_iy  = self.find_closest_index(self.position.y, self.position.y[self.HW_y[0]*0.8 < self.position.y][0])
            self.__r_iy  = self.find_closest_index(self.position.y, self.position.y[self.HW_y[1]*0.8 < self.position.y][0])
            #--- caluclate homgeinity defined as (D_may-Dmin)/D0                 
            self.H_y = (max(self.y_profile.norm[self.__l_iy : self.__r_iy+1])-min(self.y_profile.norm[self.__l_iy : self.__r_iy+1]))/self.D0_y 
            self.plateau_limits_y = np.array([self.__l_iy, self.__r_iy])
          
    def __set_symmetry(self, AXIS): 
        """Assumes voxels have constant width and are symmetric around 0 -> otherwise linear interpolation on a fixed grid would be necessary"""
        if AXIS =="X":
            ix_0 = self.find_closest_index(self.position.x, 0) #get index at x=0 then calculate abs(D(x)-D(-x))
            N = (min([ix_0-self.__l_ix, self.__r_ix-ix_0])) # Number of steps for which to calculate the Difference
            delta_D = [abs(self.x_profile.norm[ix_0-i]-self.x_profile.norm[ix_0+i]) for i in range(1,N+1)]
            self.S_x = max(delta_D)/self.D0_x
                
        if AXIS =="Y":
            ix_0 = self.find_closest_index(self.position.y, 0) # get index at x=0 then calculate abs(D(x)-D(-x))
            N = (min([ix_0-self.__l_iy, self.__r_iy-ix_0]))    # Number of steps for which to calculate the Difference
            delta_D = [abs(self.y_profile.norm[ix_0-i]-self.y_profile.norm[ix_0+i]) for i in range(1,N+1)]
            self.S_y = max(delta_D)/self.D0_y
        
    def __set_dCAX(self, AXIS):
       if AXIS =="X":
           self.dCAX_x = self.HW_x[0] + (self.HW_x[1] - self.HW_x[0])/2
       if AXIS =="Y":
           self.dCAX_y = self.HW_y[0] + (self.HW_y[1] - self.HW_y[0])/2
    
    
    
class dose_3d(dose_object):
    """
        A class that enables one to read .3ddose-Files with python and read certain profiles or arbitrary voxels.
        For more Info see Github 'https://github.com/Apelova/EGS_DOSE_TOOLS'.
    """
    def __init__(self, PATH, INFO=False):
        dose_object.__init__(self) # maintain Attributes of dose_object
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
        self.voxel_in_axis = _xyz_array_Object(self.voxel_in_axis)
        self.boundaries = _xyz_array_Object(self.boundaries)
        self.position = _xyz_array_Object(self.position)
        
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
            self.pdd = normable_array(self.dose_matrix[:, y_index, x_index])
            self.pdd_error = self.error_matrix[:, y_index, x_index]
         
            if len(self.pdd) > 1: #when reading multiple 3ddose files this makes sure profile metrics are set correctly
                self.set_pdd_metrics()
                
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
            self.x_profile = normable_array(self.dose_matrix[z_index, y_index, :])
            self.x_profile_error = self.error_matrix[z_index, y_index, :]
            self.profile_depth_x = self.position.z[z_index]
            
            if len(self.x_profile) > 1: #when reading multiple 3ddose files this makes sure profile metrics are set correctly
                self.set_metrics_profile("X")
                
            if not MUTE:
                print(f"Loaded X-Profile at the depth Z={self.position.z[z_index]} going through Y={self.position.y[y_index]}!\n")
        else:
            raise TypeError("The Parameters Z and Y accept only integer and Float values !")
        
        #--- indices to dose at 10
        
    def set_y_profile(self, Z=0, X=0, MUTE=False):
        """ default depth for doseprofile is 5cm"""
        if type(Z) and type(X) in [float, int]:
            if self.has_multiple_inputs:
                raise NotImplementedError("The Instance has multiple .3ddose Inputs and thus this function is unavailable. To get the profile adde use .y_profile!")
            
            x_index = self.find_closest_index(self.position.x, X)
            z_index = self.find_closest_index(self.position.z, Z)
            self.y_profile = normable_array(self.dose_matrix[z_index, :, x_index])
            self.y_profile_error = self.error_matrix[z_index, :, x_index]
            self.profile_depth_y = self.position.z[z_index]
            
            if len(self.y_profile) > 1: #when reading multiple 3ddose files this makes sure profile metrics are set correctly
                self.set_metrics_profile("Y")
                
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
            self.set_metrics_profile("X")
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
            self.set_metrics_profile("Y")
            
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


class dose_mcc(dose_object):
    def __init__(self, PATH, INFO=False):
        dose_object.__init__(self) # maintain Attributes of dose_object        self.origin = PATH
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
            self.all_scans[i].DOSE = normable_array(np.array(self.all_scans[i].DATA)[:,1])
            try:
                self.all_scans[i].ERROR= np.array(self.all_scans[i].DATA)[:,2]
            except: #some meassurements dont store the error !
                self.all_scans[i].ERROR= np.array([np.nan for i in self.all_scans[i].DOSE])
            #--- invert data such that coordinate systems match !
            if "-" in self.all_scans[i].AXIS:
                self.all_scans[i].POSITION = self.all_scans[i].POSITION[:-1]*-1 #necessary if scan is not symmetric !
                self.all_scans[i].DOSE = normable_array(self.all_scans[i].DOSE[::-1])
                self.all_scans[i].ERROR = self.all_scans[i].ERROR[::-1]
                #after fixing coordinate set axis label -X/-Y/-Z -> X/Y/Z
                self.all_scans[i].AXIS = self.all_scans[i].AXIS[1::]
                
    def __get_axis_label(self, i):
        """ Implementation only viable for GANTRY_ANGLE==0, the axis definition is made such that x from read_mcc== x from read_3ddose"""        
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
        self.position = _xyz_array_Object({"x":[],"y":[],"z":[]})
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
        self.set_pdd_metrics()
        self.set_metrics_profile("X")
        self.set_metrics_profile("Y")
        
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
        self.set_metrics_profile("X")
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
        self.set_metrics_profile("Y")
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
   - Percentage Depth Dose loaded for values of z in [{self.position.z[0] if self.position.z.size else '-'}, {self.position.z[-1] if self.position.z.size else '-'}]
   - Available x-Profile depths are {[ i [1] for i in self.available_depths_x]} at {[ i[0] for i in self.available_depths_x]} 
   - Available y-Profile depths are {[ i [1] for i in self.available_depths_y]} at {[ i[0] for i in self.available_depths_y]}\n
 By default the first available Profiles are set as .x_profile or .y_profile!
 To Change the currently set x/y_profile use .set_x/y_profile(Number of Scan).   
 
 The Data of each SCAN is stored together with some Information about the Scan in a dictionary .all_scans. 
 The Keys if this dictionary are the number of Scans as integer values.
 
 To save a profile one can just use the Pandas Datafram class e.g pd.DataFrame(dose_3d_object.y_profile).to_csv(destination_path)      
##############################################################################"""  


#######################
#
# Testing
#
#######################
# teste sowohl chamber als auch dosxyz 3ddose results
path = "/home/marvin/Desktop/master"
path = "C:/Users/apel04/Desktop/master"

#test_exp = dose_mcc(path+"/Messungen/6MeV_10x10_Dose_Profiles/profiles_@_5cm_ff_and_fff_SSD_100/FF/combined_ff.mcc")

if False:
    test_dosy = dose_3d(path+"/Simulationen/DATA/TEST_PTB/DOSXYZ_NRC/PTB_6MV_old_dosxyz.3ddose")

if False:
    test_chamber = dose_3d(path+"/Simulationen/DATA/TEST_PTB/EGS_CHAMBER/ROUGH_GRID/PTB_6MeVp_10x10_old_dose_pdd.3ddose")
    test_chamber.add_profile(path+"/Simulationen/DATA/TEST_PTB/EGS_CHAMBER/ROUGH_GRID/PTB_6MeVp_10x10_old_dose_x_profile.3ddose", AXIS="X")
    test_chamber.add_profile(path+"/Simulationen/DATA/TEST_PTB/EGS_CHAMBER/ROUGH_GRID/PTB_6MeVp_10x10_old_dose_y_profile.3ddose", AXIS="Y") 


####################
# Z-Profile metrics
####################
if False:
    fig, axs = plt.subplots(1,3, figsize=(30,10))
    axs[0].set_ylim(20, 101)
    axs[1].set_ylim(20, 101)
    axs[2].set_ylim(20, 101)
    #--- PDD-Test Metrics   
    # Experimental
    axs[0].plot(test_exp.position.z, test_exp.pdd.norm, label=f"Experiment Q: {test_exp.Q}")
    axs[0].scatter([0, 10, 20, test_exp.Rmax], [test_exp.D0, test_exp.D10, test_exp.D20, 100])
    axs[0].legend(loc="lower left", fontsize=20)
    # DOSXYZ_NRC
    axs[1].plot(test_dosy.position.z, test_dosy.pdd.norm, label=f"DOSXYZ_NRC Q: {test_dosy.Q}")
    axs[1].scatter([0, 10, 20, test_dosy.Rmax], [test_dosy.D0, test_dosy.D10, test_dosy.D20, 100])
    axs[1].legend(loc="lower left", fontsize=20)
    # EGS_CHAMBER
    axs[2].plot(test_dosy.position.z, test_dosy.pdd.norm, label=f"EGS_CHAMBER Q: {test_dosy.Q}")
    axs[2].scatter([0, 10, 20, test_dosy.Rmax], [test_dosy.D0, test_dosy.D10, test_dosy.D20, 100])
    axs[2].legend(loc="lower left", fontsize=20)

####################
# X-Profile metrics
####################    
names = ["Experiment", "DOSXYZ_NRC", "EGS_CHAMBER"]
if False:
    fig, axs = plt.subplots(1,3, figsize=(30,10))
    for i, dose in enumerate([test_exp, test_dosy, test_chamber]):
        label = f"{names[i]}\n    H = {round(dose.H_x,3)}\n    S = {round(dose.S_x,3)}\nCAX = {round(dose.dCAX_x,3)}"
        print(dose.position.x.shape, dose.x_profile.shape)
        #axs[i].scatter(dose.position.x, dose.x_profile.norm, facecolor="white",edgecolor="black", marker="o", label=label)
        #axs[i].plot(dose.HW_x, [50,50], c="red")
        #axs[i].plot([dose.penumbra_limits_x_20[0], dose.penumbra_limits_x_80[0]], [20,80], c="red", lw=2)
        #axs[i].plot([dose.penumbra_limits_x_20[1], dose.penumbra_limits_x_80[1]], [20,80], c="red", lw=2)
        #axs[i].plot([dose.dCAX_x, dose.dCAX_x], [-1, 103])
        #axs[i].plot(dose.position.x[dose.plateau_limits_x], [100, 100], lw=4, c="blue")
        #axs[i].set_ylim(0, 101)
        #axs[i].legend(loc="lower center", fontsize=25)
                
####################
# Y-Profile metrics
####################
names = ["Experiment", "DOSXYZ_NRC", "EGS_CHAMBER"]
if False:
    fig, axs = plt.subplots(1,3, figsize=(30,10))
    for i, dose in enumerate([test_exp, test_dosy, test_chamber]):
        label = f"{names[i]}\n    H = {round(dose.H_y,3)}\n    S = {round(dose.S_y,3)}\nCAX = {round(dose.dCAX_y,3)}"
        axs[i].scatter(dose.position.y, dose.y_profile.norm, facecolor="white",edgecolor="black", marker="o", label=label)
        axs[i].plot(dose.HW_y, [50,50], c="red")
        axs[i].plot([dose.penumbra_limits_y_20[0], dose.penumbra_limits_y_80[0]], [20,80], c="red", lw=2)
        axs[i].plot([dose.penumbra_limits_y_20[1], dose.penumbra_limits_y_80[1]], [20,80], c="red", lw=2)
        axs[i].plot([dose.dCAX_y, dose.dCAX_y], [-1, 103])
        axs[i].plot(dose.position.y[dose.plateau_limits_y], [100, 100], lw=4, c="blue")
        axs[i].set_ylim(0, 101)
        axs[i].legend(loc="lower center", fontsize=25)
        
# =============================================================================
# #profile Metrics - X
# print("\n\nX-Axis")
# print(test_3d.D0_x, test_mcc.D0_x)
# print(test_3d.plateau_limits_x, test_mcc.plateau_limits_x)
# print(test_3d.penumbra_limits_x, test_mcc.penumbra_limits_x)
# print(test_3d.H_x, test_mcc.H_x)
# print(test_3d.HW_x, test_mcc.HW_x)
# print(test_3d.S_x, test_mcc.S_x)
# print(test_3d.dCAX_x, test_mcc.dCAX_x)
# 
# 
# #profile Metrics - y
# print("\n\ny-Axis")
# print(test_3d.D0_y, test_mcc.D0_y)
# print(test_3d.plateau_limits_y, test_mcc.plateau_limits_y)
# print(test_3d.penumbra_limits_y, test_mcc.penumbra_limits_y)
# print(test_3d.H_y, test_mcc.H_y)
# print(test_3d.HW_y, test_mcc.HW_y)
# print(test_3d.S_y, test_mcc.S_y)
# print(test_3d.dCAX_y, test_mcc.dCAX_y)
# 
# =============================================================================