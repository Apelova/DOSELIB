"""
@author: Apelova 
https://github.com/Apelova/EGS_DOSE_TOOLS

Definition of all classes and methods used in to analyse .3ddose or .mcc files.


"""
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
import matplotlib
import os

class _xyz_array_Object:
    def __init__(self, array):
        self.x = np.array(array["x"])
        self.y = np.array(array["y"])
        self.z = np.array(array["z"])

class normable_array(np.ndarray):
    """ A normable array is the Object-Class used to represent Dose-Profiles
        The normalized profile can easily be accessed through the attribute .norm
    """
    def __new__(cls, arr):
        obj = np.asarray(arr).view(cls)
        return obj

    def __init__(self, arr):
        self.norm = self/max(self)*100    
        self.error = None

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
        self.ReferenceValue=None
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
        self.__l_iy = None
        self.__r_iy= None
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
        self.penumbra_xl= np.nan
        self.penumbra_xr= np.nan
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
        
    def norm_plateau(self):
        if self.x_profile is not None:
            if len(self.x_profile) > 1: 
                self.x_profile.norm *= 100/np.mean(self.x_profile.norm[self.__l_ix: self.__r_ix])

        if self.y_profile is not None:
            if len(self.y_profile) > 1:
                self.y_profile.norm *= 100/np.mean(self.y_profile.norm[self.__l_iy: self.__r_iy])

    def norm_on_center(self):
        if self.x_profile is not None:
            center_x = self.find_closest_index(self.position.x, 0)
            self.x_profile.norm = self.x_profile.norm/self.x_profile.norm[center_x]*100 

        if self.y_profile is not None:
            center_y = self.find_closest_index(self.position.y, 0)
            self.y_profile.norm = self.y_profile.norm/self.y_profile.norm[center_y]*100 

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
        try:
            self.__set_Dose_on_axis(AXIS)
            self.__set_penumbra_limits_20(AXIS)
            self.__set_penumbra_limits_80(AXIS)
            self.__set_penumbra_width(AXIS)
            self.__set_halfwidth(AXIS)
            self.__set_dCAX(AXIS)
            self.__set_homogeneity(AXIS)
            self.__set_symmetry(AXIS)
        except:
            print(f"Warning --- setting metrics for {AXIS}-axes failed !")

        
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
    def __init__(self, PATH, INFO=False, FIRST_PDD_VOXEL_IS_HALF_VOLUME=False, normalize_x_y_on_max=False, normalize_profiles_on_center=False):
        dose_object.__init__(self) # maintain Attributes of dose_object
        self.origin = PATH
        self.normed_on_max    = normalize_x_y_on_max
        self.normed_on_center = normalize_profiles_on_center
    
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
            #--- normalize X/Y Profiles on start
            if not self.normed_on_max and not self.normed_on_center:
                self.norm_plateau()
            elif self.normed_on_center:
                self.norm_on_center()
            #--- correct dose as Eabs/m when m[0] is 0.5 m[i] -> i index for all voxels without 0
            if FIRST_PDD_VOXEL_IS_HALF_VOLUME:
                self.pdd[0] = 2*self.pdd[0]
                self.pdd.norm[0] = 2*self.pdd.norm[0]
            #--- output information on read
            if INFO:
                print(self)
        else:
            raise TypeError("The Input PATH is invalid!" )
       
        print(" ")

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
        if self.boundaries["z"][0] < 0:
            self.boundaries["z"][0] = np.mean(self.boundaries["z"][0:2])
            print(f"Warning --- First z-Boundary is negative ! Therfore it is replaced with {self.boundaries['z'][0]}. Make sure to test if you have to set FIRST_PDD_VOXEL_IS_HALF_VOLUME==True")
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
            self.pdd.error = self.error_matrix[:, y_index, x_index]
      
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
            self.x_profile.error = self.error_matrix[z_index, y_index, :]
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
            self.y_profile.error = self.error_matrix[z_index, :, x_index]
            self.profile_depth_y = self.position.z[z_index]
            
            if len(self.y_profile) > 1: #when reading multiple 3ddose files this makes sure profile metrics are set correctly
                self.set_metrics_profile("Y")
                
            if not MUTE:
                print(f"Loaded Y-Profile at the depth Z={self.position.z[z_index]} going through X={self.position.x[x_index]}!\n")
        else:
            raise TypeError("The Parameters Z and Y accept only integer and Float values !")
        
    #--- getter functions
    def get_pdd(self):
        """ return the currently set Z-Profile/(PDD) and the codoserresponding statistical error"""
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
        profile = dose_3d(PATH, normalize_x_y_on_max=self.normed_on_max, normalize_profiles_on_center=self.normed_on_center)
        if AXIS.lower()=="x":
            if len(profile.position.x) <= 10:
                print(f"Warning---The X-profile you want to add has only {len(profile.position.x)} Entries. I hope you know what you are doing.")
                #--- set the initial positions for x and y such that adding profiles wont allow setting faulty z-profiles (when adding the first profile only)
            self.has_multiple_inputs = True
            self.boundaries.x = profile.boundaries.x
            self.position.x = profile.position.x
            self.x_profile = profile.x_profile
            self.x_profile_error = profile.x_profile_error
            self.x_profile.error = profile.x_profile_error
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
            self.y_profile.error = profile.y_profile_error
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
    def __init__(self, PATH, INFO=False,  normalize_x_y_on_max=False, average_profiles=False, consider_reference=False, normalize_profiles_on_center=False):
        dose_object.__init__(self) # maintain Attributes of dose_object        self.origin = PATH
        self.origin = PATH
        self.average_profiles_at_init = average_profiles
        self.consider_reference = consider_reference
        self.has_reference = False
        #--- extract data from scans and store the raw information in a dictionary .all_scans
        self.__load_scans()
        #--- scan to dose object
        self.__set_scan_parameter()
        #--- set default profiles
        #--- load available profiles PDD, X/Y-Profiles
        self.__load_profiles()
        #--- norm plateau 
        if not normalize_x_y_on_max and not normalize_profiles_on_center:
            self.norm_plateau()            
        elif normalize_profiles_on_center:
            self.norm_on_center()
        if INFO:
            print(self)  
        print(" ")
            
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
                self.all_scans[i].ReferenceValue = normable_array(np.array(self.all_scans[i].DATA)[:,2])
                self.has_reference = True
                print("MCC-File contains either reference values or experimental errors. Call consider_reference to fix normalized profiles.")
            except:
                try:
                    self.all_scans[i].ERROR= np.array(self.all_scans[i].DATA)[:,2]
                except: #some meassurements dont store the error !
                    self.all_scans[i].ERROR= np.array([np.nan for i in self.all_scans[i].DOSE])
            
            if self.has_reference:
                try:
                    self.all_scans[i].ERROR= np.array(self.all_scans[i].DATA)[:,3]
                except: #some meassurements dont store the error !
                    self.all_scans[i].ERROR= np.array([np.nan for i in self.all_scans[i].DOSE])
            
                    
                
    def __get_axis_label(self, i):
        """ Implementation only viable for GANTRY_ANGLE==0, the axis definition is made such that x from read_mcc== x from read_3ddose"""        
        if self.all_scans[i].GANTRY_angle=="0.00":
            if self.all_scans[i].SCAN_curvetype == "CROSSPLANE_PROFILE":
                return "Y"
            
            elif self.all_scans[i].SCAN_curvetype == "INPLANE_PROFILE":
                return "X"
    
            else:
                return "undefined"
            
        else:
            return "undefined"
        
    def __load_profiles(self):
        self.pdd, self.x_profile, self.y_profile= np.array([]), np.array([]), np.array([])
        self.position = _xyz_array_Object({"x":[],"y":[],"z":[]})
        self.available_depths_x, self.available_depths_y = [], []
        
        if self.consider_reference:
            print("--- Reference will be considered in the normalized profiles ---")
            
        for i in self.all_scans:
            action = "loaded"                

            if self.all_scans[i].AXIS == "Z" and self.pdd.size==0:
                self.position.z = self.all_scans[i].POSITION
                self.pdd = self.all_scans[i].DOSE
                self.pdd_error = self.all_scans[i].ERROR
                self.pdd.error = self.all_scans[i].ERROR
                self.set_pdd_metrics()
                if self.consider_reference: #would be normalized on plateau so just normalize on max here !
                    try:    
                        self.pdd = normable_array(self.all_scans[i].DOSE/self.all_scans[i].ReferenceValue)                    
                    except:
                        print("--- Reference Values for Z-Axis are missing. ! Therefore they are not considered ! ---")
                
                print("Successfully set Percentage Depth Dose (PDD)!")
            
            elif self.all_scans[i].AXIS == "X":
                if self.x_profile.size==0:
                    self.position.x = self.all_scans[i].POSITION
                    self.x_profile = self.all_scans[i].DOSE                  
                    self.x_profile_error = self.all_scans[i].ERROR
                    self.x_profile.error = self.all_scans[i].ERROR
                    self.set_metrics_profile("X")
                    action = "set"
                    if self.consider_reference: #would be normalized on plateau so just normalize on max here !
                        try:    
                            self.x_profile = normable_array(self.all_scans[i].DOSE/self.all_scans[i].ReferenceValue)

                        except:
                            print("--- Reference Values for X-Axis are missing. ! Therefore they are not considered ! ---")

                        
                self.available_depths_x.append((f"Number of Scan {i}",float(self.all_scans[i].SCAN_DEPTH)/10))
                print(f"Successfully {action} X-Profile at depth Z={float(self.all_scans[i].SCAN_DEPTH)/10}cm!")
            
            elif self.all_scans[i].AXIS == "Y":
                if self.y_profile.size==0:
                    self.position.y = self.all_scans[i].POSITION
                    self.y_profile = self.all_scans[i].DOSE
                    self.y_profile_error = self.all_scans[i].ERROR
                    self.y_profile.error = self.all_scans[i].ERROR
                    self.set_metrics_profile("Y")
                    action = "set"
                    if self.consider_reference: #would be normalized on plateau so just normalize on max here !
                        try:    
                            self.y_profile = normable_array(self.all_scans[i].DOSE/self.all_scans[i].ReferenceValue)
                        except:
                            print("--- Reference Values for Y-Axis are missing. ! Therefore they are not considered ! ---")

                self.available_depths_y.append((f"Number of Scan {i}",float(self.all_scans[i].SCAN_DEPTH)/10))
                print(f"Successfully {action} Y-Profile at depth Z={float(self.all_scans[i].SCAN_DEPTH)/10}cm!")
            
            else:
                print(f"Warning---Axis of Scan {i} is unknown!")

        if self.average_profiles_at_init:
            self.average_profile()

        
    def average_profile(self):
        """ only works if profile is meassured symmetrically !"""
        if self.position.x.size:
            if (self.position.x.shape[0]%2==1):
                Nx = self.position.x.shape[0]//2 
                averaged_x = (self.x_profile.norm[:Nx] + self.x_profile.norm[:Nx:-1])/2
                self.x_profile.norm[:Nx] = averaged_x
                self.x_profile.norm[:Nx:-1] = averaged_x
                self.x_profile.norm = self.x_profile.norm/max(self.x_profile.norm)*100
            else:
                print("Warning --- averaging of X-Profile impossible (assymetric meassurements)")

        if self.position.y.size:
            if (self.position.y.shape[0]%2==1):
                Ny = self.position.y.shape[0]//2 
                averaged_y = (self.y_profile.norm[:Ny] + self.y_profile.norm[:Ny:-1])/2
                self.y_profile.norm[:Ny] = averaged_y
                self.y_profile.norm[:Ny:-1] = averaged_y
                self.y_profile.norm = self.y_profile.norm/max(self.y_profile.norm)*100
            else:
                print("Warning --- averaging of Y-Profile impossible (assymetric meassurements)")

    
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
        self.x_profile.error = self.all_scans[NUMSCAN].ERROR
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
        self.y_profile.error = self.all_scans[NUMSCAN].ERROR
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
        
    #Plot and show comparison
    if show_plt:
        axs[2].plot(axes_for_plot, gamma_array, label="γ(x)", lw=2, c="#42a1f5")
        axs[2].legend(loc="upper right", fontsize=20)

    #Step 4 calculate derivatives
    #--------------------------------------------------------------------------        
    if derivatives:
        dGamma_dx = np.zeros(len(gamma_array))
        dGamma_dD = np.zeros(len(gamma_array))
        #TODO! additionally use numba !    
        #--------------------------------------------------------------------------        
        for i in range(len(axes_reference)):
             if i!=0 and i!=len(axes_reference)-1:
                # O(delta**2) - Approximations
                dGamma_dx[ix+i] = (gamma_array[ix+i+1]-gamma_array[ix+i-1])/(2*(axes_reference[i+1]-axes_reference[i-1])) \
                                 if (axes_reference[i+1]-axes_reference[i-1]!=0) else 0     
                 
                dGamma_dD[ix+i] = (gamma_array[ix+i+1]-gamma_array[ix+i-1])/(2*(dose_reference[i+1]-dose_reference[i-1])) \
                                 if (dose_reference[i+1]-dose_reference[i-1]!=0) else 0    
        if show_plt:
            pass
            axs[3].plot(axes_for_plot, dGamma_dx, label="dγ/dx(x)", c="#4e42f5")
            axs[3].plot(axes_for_plot, dGamma_dD, label="dγ/dD(x)", c="red")
            axs[3].legend(loc="upper right", fontsize=20)        
        return gamma_array, dGamma_dx, dGamma_dD
            
    return gamma_array  

def rename_attr(name):
    name_map = {
        "D0_x": "D0",
        "D0_y": "D0",
        "H_x": "Hom",
        "H_y": "Hom",
        "S_x": "S",
        "S_y": "S",
        "dCAX_x": "ΔCAX",
        "dCAX_y": "ΔCAX",
        "penumbra_xl": "Left π",
        "penumbra_xr": "Right π",
        "penumbra_yl": "Left π",
        "penumbra_yr": "Right π",}
    return name_map[name] if name in name_map else name

def get_metric_string(DOSE_OBJ, AXIS, exclude={}):
    """return a string formated to fit the box"""
    # ---
    ATTRIBUTES = {"Z": ["D0", "Rmax", "D10", "D20", "Q"],
                  "X": ["D0_x", "penumbra_xl", "penumbra_xr", "H_x", "S_x", "dCAX_x"],
                  "Y": ["D0_y", "penumbra_yl", "penumbra_yr", "H_y", "S_y", "dCAX_y"]}
    UNITS = {
        "D0": "% ",
        "Rmax": "cm",
        "D10": "% ",
        "D20": "% ",
        "Q": "  ",
        "D0_x": "% ",
        "penumbra_xl": "cm",
        "penumbra_xr": "cm",
        "H_x": "  ",
        "S_x": "  ",
        "dCAX_x": "cm",
        "D0_y": "% ",
        "penumbra_yl": "cm",
        "penumbra_yr": "cm",
        "H_y": "  ",
        "S_y": "  ",
        "dCAX_y": "cm" } 
    def add_digits_to_string(substring):
        output = ""
        if "-0." in substring:
            output += " "
        elif len(substring) == 4:
            output +="  "
        elif len(substring)<6:
            output += " "
            if "0." in substring:
                output += " "                            
        
        output += substring
        return output

    if AXIS in exclude:
        ATTRIBUTES[AXIS] = [attr for attr in ATTRIBUTES[AXIS] if attr not in exclude[AXIS]] 

    #--- set amount of space such that "attr:" is equidistant (: at same position)
    dL = max([len(rename_attr(attr)) for attr in ATTRIBUTES[AXIS]])
    # --- gather requested  information in string and output
    metric_string = ""
    for attr in ATTRIBUTES[AXIS]:
        if AXIS.upper() == "Z":
            if attr == "D0":
                metric_string += f"{rename_attr(attr):<{dL-1}}: "
            else:
                metric_string += f"{rename_attr(attr):<{dL}}: "
            metric_string += add_digits_to_string(f"{getattr(DOSE_OBJ, attr):.2f}") 
            metric_string += f"{UNITS[attr]}\n"
        else:
            metric_string += f"{rename_attr(attr):<{dL}}: "
            metric_string += add_digits_to_string(f"{getattr(DOSE_OBJ, attr):.2f}")
            metric_string += f"{UNITS[attr]}\n"
    # --- remove last "\n"
    return metric_string[:-1]


def set_metrics(dose_objs, plot_axs, axes=["Z", "X", "Y"], difference=True, exclude={}, colors=["black", "blue"]):
    dh_x = 6-len(exclude["X"]) if "X" in exclude else 6
    dh_y = 6-len(exclude["Y"]) if "Y" in exclude else 6
    for  i, dose in enumerate(dose_objs[:2]):#compare maximum two dose objs !
        for j, ax in enumerate(axes):
            if difference:
                # Match statements require python 3,8 -> not available in the office therefore use this shananigans
                if ax == "Z":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="Z", exclude=exclude)
                    if len(axes)==1:
                        plot_axs[0].text(0.69+i*0.29, 0.94, metric_string, transform=plot_axs[0].transAxes, fontsize=19, bbox=props, ha="right", va="top", font="monospace")
                    else:
                        plot_axs[0,j].text(0.69+i*0.29, 0.94, metric_string, transform=plot_axs[0,j].transAxes, fontsize=19, bbox=props, ha="right", va="top", font="monospace")
                elif ax =="X":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="X", exclude=exclude)
                    if len(axes)==1:
                        plot_axs[0].text(0.98, 0.94-i*0.05*dh_x, metric_string, transform=plot_axs[0].transAxes, fontsize=19, bbox=props, ha="right", va="top", font="monospace")
                    else:
                        plot_axs[0,j].text(0.98, 0.94-i*0.05*dh_x, metric_string, transform=plot_axs[0,j].transAxes, fontsize=19, bbox=props, ha="right", va="top", font="monospace")
                elif ax =="Y":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="Y", exclude=exclude)
                    if len(axes)==1:
                        plot_axs[0].text(0.98, 0.94-i*0.05*dh_y, metric_string, transform=plot_axs[0].transAxes, fontsize=19, bbox=props, ha="right", va="top", font="monospace")
                    else:
                        plot_axs[0,j].text(0.98, 0.94-i*0.05*dh_y, metric_string, transform=plot_axs[0,j].transAxes, fontsize=19, bbox=props, ha="right", va="top", font="monospace")
            else:
                if ax =="Z":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="Z", exclude=exclude)
                    if len(axes)==1:
                        plot_axs.text(0.76+i*0.22, 0.94, metric_string, transform=plot_axs.transAxes, fontsize=15, bbox=props, ha="right", va="top", font="monospace")
                    else:
                        plot_axs[j].text(0.76+i*0.22, 0.94, metric_string, transform=plot_axs[j].transAxes, fontsize=15, bbox=props, ha="right", va="top", font="monospace")
                        
                elif ax =="X":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="X", exclude=exclude)
                    if len(axes)==1:
                        plot_axs.text(0.98, 0.94-i*0.06*(dh_x+1), metric_string, transform=plot_axs.transAxes, fontsize=15, bbox=props, ha="right", va="top", font="monospace")
                    else:
                        plot_axs[j].text(0.98, 0.94-i*0.06*(dh_x+1), metric_string, transform=plot_axs[j].transAxes, fontsize=15, bbox=props, ha="right", va="top", font="monospace")
                elif ax =="Y":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="Y", exclude=exclude)
                    if len(axes)==1:
                        plot_axs.text(0.98, 0.94-i*0.06*(dh_y+1), metric_string, transform=plot_axs.transAxes, fontsize=15, bbox=props, ha="right", va="top", font="monospace")
                    else:
                        plot_axs[j].text(0.98, 0.94-i*0.06*(dh_y+1), metric_string, transform=plot_axs[j].transAxes, fontsize=15, bbox=props, ha="right", va="top", font="monospace")


def compare_dose( dose_objects, labels=None, axes=["Z", "X", "Y"], difference=True, figsize=None, metrics=False, exclude= {}, colors = ["black", "red", "blue"], 
                  interpol="linear", diff_dx=0.1, scatter_first=True, show_labels=False, labels_loc="upper right", label_size=15, norm_on_max=False, plot_raw=False):
    """
        A function that return plot of selected axes.
        The limits can be set by calling the output !
        
        If difference is True and (dose_objects) len(2) it assumes the difference to the first element of dose_objects should be calculated for the others
        
        metrics = "average", "profiles"
        
    """ 
    if labels is None:
        labels = [i for i in range(len(dose_objects))]
        
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
        raise TypeError(f"Invalid input for Axis! THe following Axes are invalid {[ each for each in axes  if each not in ['X','Y','Z'] ]} Valid Inputs are [X, Y, Z] ")
   
    if len(colors) < len(dose_objects):
        import matplotlib.colors as mcolors
        if len(dose_objects) > 10:
            colors = [ rgb[1] for rgb in mcolors.CSS4_COLORS.items()]
        else:
            colors = [ rgb[1] for rgb in mcolors.TABLEAU_COLORS.items()]
            

    if len(labels) < len(dose_objects) and show_labels:
        raise ValueError("Missing Input labels (list with labels for plots has to be greated or equal the length of doseobjects to compare)")


    #--- set variables for layout
    cols, rows = (len(axes)), (int(difference)+1)

    #--- set figsize depending on requested graphs
    if not figsize:
        figsize = (8*(cols), (4*rows))
        figsize = (10*(cols), (5*rows))
    
    matplotlib.rcParams['axes.grid'] = True
    #--- add requested plots
    if difference:
        fig_, axs_ = plt.subplots(rows, cols, figsize=figsize, sharex="col", constrained_layout=True, gridspec_kw={'height_ratios': [4,1]})
    else:
        fig_, axs_ = plt.subplots(rows, cols, figsize=figsize, sharex="col", constrained_layout=True)


    #--- define helping functions !
    def plot_difference(fig_, axs_, position_dose_pair_arrays, labels, metrics=None, colors = ["black", "red", "blue"], difference=True, plot_raw=False):
        #--- group positions and pairs together and interpolate dose on a common grid !
        common_position, interpolated_dose = interpolate_arrays( [pair[0] for pair in position_dose_pair_arrays], [pair[1] for pair in position_dose_pair_arrays])
        #--- iterate through all the AXIS to plot !
        if plot_raw:
            for i in range(1, len(interpolated_dose)):
                if len(axes)==1:
                    axs_[1].plot(common_position, (interpolated_dose[i]-interpolated_dose[0])/interpolated_dose[0]*100, c=colors[i], lw=2)
                else:
                    axs_[1,col].plot(common_position, (interpolated_dose[i]-interpolated_dose[0])/interpolated_dose[0]*100, c=colors[i], lw=2)            
        else:
            for i in range(1, len(interpolated_dose)):
                if len(axes)==1:
                    axs_[1].plot(common_position, (interpolated_dose[i]-interpolated_dose[0])/interpolated_dose[0]*100, c=colors[i], lw=2)
                else:
                    axs_[1,col].plot(common_position, (interpolated_dose[i]-interpolated_dose[0])/interpolated_dose[0]*100, c=colors[i], lw=2)
            
        
    #--- define helping functions !
    def plot_dose_distribution(fig_, axs_, position_dose_pair_arrays, labels, metrics=None, colors = ["black", "red", "blue"], difference=True,scatter_first=True, plot_raw=False):
        #--- for each plot of the pairs 
        for i, pair in enumerate(position_dose_pair_arrays):
            # the first one is the reference -> make it appear different!
            if i==0 and scatter_first:
                if difference:
                    if len(axes)==1:
                        axs_[0].scatter(pair[0], pair[1], c=colors[i], label=labels[i], marker="o", s=45, alpha=0.75)
                    else:
                        axs_[0, col].scatter(pair[0], pair[1], c=colors[i], label=labels[i], marker="o", s=45, alpha=0.75)
                else:
                    if len(axes)==1:
                        axs_.scatter(pair[0], pair[1], c=colors[i], label=labels[i], marker="o", s=45, alpha=0.75)
                    else:
                        axs_[col].scatter(pair[0], pair[1], c=colors[i], label=labels[i], marker="o", s=45, alpha=0.75)
            else:                                       
                if difference:
                    if len(axes)==1:
                        axs_[0].plot(pair[0], pair[1], c=colors[i], lw=2, label=labels[i])
                    else:
                        axs_[0, col].plot(pair[0], pair[1], c=colors[i], lw=2, label=labels[i])
                else:
                    if len(axes)==1:
                        axs_.plot(pair[0], pair[1], c=colors[i], lw=2, label=labels[i])
                    else:
                        axs_[col].plot(pair[0], pair[1], c=colors[i], lw=2, label=labels[i])
            
        if show_labels:
            if difference:
                if len(axes) > 1:
                    axs_[0,0].legend(loc=labels_loc, fontsize=label_size)
                else:
                    axs_[0].legend(loc=labels_loc, fontsize=label_size)
            else:
                if len(axes) > 1:
                    axs_[0].legend(loc=labels_loc, fontsize=label_size)
                else:
                    axs_.legend(loc=labels_loc, fontsize=label_size)

        if difference:
            plot_difference(fig_, axs_, position_dose_pair_arrays, labels, metrics, colors, plot_raw=plot_raw)

    if plot_raw:
        for col in range(cols):
            #--- print original plots
            if axes[col] == "X":
                plot_dose_distribution(fig_, axs_, [[each.position.x, each.x_profile] for each in dose_objects], labels, difference=difference, colors=colors, scatter_first=scatter_first, plot_raw=True)
    
            if axes[col] == "Y":
                plot_dose_distribution(fig_, axs_, [[each.position.y, each.y_profile] for each in dose_objects], labels, difference=difference, colors=colors, scatter_first=scatter_first, plot_raw=True)
        
            if axes[col] == "Z":
                plot_dose_distribution(fig_, axs_, [[each.position.z, each.pdd] for each in dose_objects], labels, difference=difference, colors=colors, scatter_first=scatter_first, plot_raw=True)        

    else:
        for col in range(cols):
            #--- print original plots
            if axes[col] == "X":
                if norm_on_max:
                    plot_dose_distribution(fig_, axs_, [[each.position.x, each.x_profile/max(each.x_profile)*100] for each in dose_objects], labels, difference=difference, colors=colors, scatter_first=scatter_first)
                else:                
                    plot_dose_distribution(fig_, axs_, [[each.position.x, each.x_profile.norm] for each in dose_objects], labels, difference=difference, colors=colors, scatter_first=scatter_first)
    
            if axes[col] == "Y":
                if norm_on_max:
                    plot_dose_distribution(fig_, axs_, [[each.position.y, each.y_profile/max(each.y_profile)*100] for each in dose_objects], labels, difference=difference, colors=colors, scatter_first=scatter_first)
                else:
                    plot_dose_distribution(fig_, axs_, [[each.position.y, each.y_profile.norm] for each in dose_objects], labels, difference=difference, colors=colors, scatter_first=scatter_first)
        
            if axes[col] == "Z":
                if norm_on_max:
                    plot_dose_distribution(fig_, axs_, [[each.position.z, each.pdd/max(each.pdd)*100] for each in dose_objects], labels, difference=difference, colors=colors, scatter_first=scatter_first)
                else:
                    plot_dose_distribution(fig_, axs_, [[each.position.z, each.pdd.norm] for each in dose_objects], labels, difference=difference, colors=colors, scatter_first=scatter_first)


    #--- add x-labels
    axes_labels = {"X": "Position along X-Axis [cm]", "Y": "Position along Y-Axis [cm]", "Z": "Position along Z-Axis [cm]"}
    for i, each in enumerate(axes):
        if difference:
            if len(axes)==1:
                axs_[1].set_xlabel(axes_labels[each], fontsize=18)
            else:
                axs_[1,i].set_xlabel(axes_labels[each], fontsize=18)
        else:
            if len(axes)==1:
                axs_.set_xlabel(axes_labels[each], fontsize=18)
            else:
                axs_[i].set_xlabel(axes_labels[each], fontsize=18)

    #--- add y-labels
    if difference:
        if len(axes)==1:
            axs_[0].set_ylabel("relative Dose [%]", fontsize=18)
            axs_[1].set_ylabel("Δ [%]", fontsize=18)
        else:
            axs_[0,0].set_ylabel("relative Dose [%]", fontsize=18)
            axs_[1,0].set_ylabel("Δ [%]", fontsize=18)
    else:
        if len(axes)==1:
            axs_.set_ylabel("relative Dose [%]", fontsize=18)
        else:            
            axs_[0].set_ylabel("relative Dose [%]", fontsize=18)

    #---adjust labels
    matplotlib.rc('xtick', labelsize=18)
    matplotlib.rc('ytick', labelsize=18)

        
    if metrics:
        set_metrics(dose_objects[:2], axs_, axes=axes, difference=difference, colors=colors, exclude=exclude)
        
    #--- return fig and axs
    return fig_, axs_


def gamma_plot(reference, to_be_evaluated, dD=3, dx=0.3, swap_scatter=False, just_profiles=False, scatter_first=True, CUTOFF=20):
    if just_profiles:
        pass
    else:
        gamma_z = gamma(
                        axes_reference = reference.position.z, dose_reference  = reference.pdd.norm, 
                        axes_evaluation= to_be_evaluated.position.z, dose_evaluation = to_be_evaluated.pdd.norm, 
                        CUTOFF=0, dD=dD, dx=dx) #3% und 3mm
        Gamma_PR_Z = round(len(gamma_z[gamma_z<1])/len(gamma_z)*100,2)

    #X-Axis
    gamma_x =gamma(
                    axes_reference = reference.position.x, dose_reference  = reference.x_profile.norm, 
                    axes_evaluation= to_be_evaluated.position.x, dose_evaluation = to_be_evaluated.x_profile.norm, 
                    dD=dD, dx=dx, CUTOFF=CUTOFF) #3% und 3mm
    #Y-Axis
    gamma_y = gamma(
                    axes_reference = reference.position.y, dose_reference  = reference.y_profile.norm, 
                    axes_evaluation= to_be_evaluated.position.y, dose_evaluation = to_be_evaluated.y_profile.norm, 
                    dD=dD, dx=dx, CUTOFF=CUTOFF) #3% und 3mm


    Gamma_PR_x = round(len(gamma_x[gamma_x<1])/len(gamma_x)*100,2)
    Gamma_PR_y = round(len(gamma_y[gamma_y<1])/len(gamma_y)*100,2)

    if just_profiles:
        fig, axs = plt.subplots(2,2, figsize=(25, 9), sharex="col")
    else:
        fig, axs = plt.subplots(2,3, figsize=(25, 9), sharex="col")
        axs[0,2].set_title(f"γ(z) Passing Rate = {Gamma_PR_Z}%", font="monospace", fontsize=20)
        axs[0,2].plot(reference.position.z, gamma_z, c="black")
    fig.tight_layout()

    #Z-Axis
    if swap_scatter:
        axs[1,0].plot(reference.position.x, reference.x_profile.norm, c="red", lw=2)
        axs[1,1].plot(reference.position.y, reference.y_profile.norm, c="red", lw=2)
        if scatter_first:
            axs[1,0].scatter(to_be_evaluated.position.x, to_be_evaluated.x_profile.norm, c="black", marker="o", s=45, alpha=0.75)
            axs[1,1].scatter(to_be_evaluated.position.y, to_be_evaluated.y_profile.norm, c="black", marker="o", s=45, alpha=0.75)
        else:
            axs[1,0].plot(to_be_evaluated.position.x, to_be_evaluated.x_profile.norm, c="black", lw=2)
            axs[1,1].plot(to_be_evaluated.position.y, to_be_evaluated.y_profile.norm, c="black", lw=2)
            
        if just_profiles:
            pass
        else:
            axs[1,2].plot(reference.position.z, reference.pdd.norm, c="red", lw=2)
            axs[1,2].scatter(to_be_evaluated.position.z, to_be_evaluated.pdd.norm, c="black", marker="o", s=45, alpha=0.75)

    else:
        axs[1,0].plot(to_be_evaluated.position.x, to_be_evaluated.x_profile.norm, c="red", lw=2)
        axs[1,1].plot(to_be_evaluated.position.y, to_be_evaluated.y_profile.norm, c="red", lw=2)
        if scatter_first:
            axs[1,0].scatter(reference.position.x, reference.x_profile.norm, c="black", marker="o", s=45, alpha=0.75)
            axs[1,1].scatter(reference.position.y, reference.y_profile.norm, c="black", marker="o", s=45, alpha=0.75)
        else:
            axs[1,0].scatter(reference.position.x, reference.x_profile.norm, c="black", marker="o", s=45, alpha=0.75)
            axs[1,1].scatter(reference.position.y, reference.y_profile.norm, c="black", lw=2)
            
        if just_profiles:
            pass
        else:
            axs[1,2].scatter(reference.position.z, reference.pdd.norm, c="black", marker="o", s=45, alpha=0.75)
            axs[1,2].plot(to_be_evaluated.position.z, to_be_evaluated.pdd.norm, c="red", lw=2)
        
    #X-Axis
    axs[0,0].plot(reference.position.x, gamma_x, c="black")
    axs[0,0].set_title(f"γ(x) Passing Rate = {Gamma_PR_x}%", font="monospace", fontsize=20)
    #Y-Axis
    axs[0,1].plot(reference.position.y, gamma_y, c="black")
    axs[0,1].set_title(f"γ(y) Passing Rate = {Gamma_PR_y}%", font="monospace", fontsize=20)
    
    return fig, axs

class dose_chamber_log(dose_object):
    """
        A class that enables one to read dosevalues from .egslog-Files that follow a specific labeling scheme
        For more Info see Github 'https://github.com/Apelova/EGS_DOSE_TOOLS'.
    """
    def __init__(self, PATH, INFO=False, normalize_x_y_on_max=False, delete_deviants=False, error_limit=10):
        dose_object.__init__(self) # maintain Attributes of dose_object
        self.origin = PATH
        self.normed_on_max = normalize_x_y_on_max
        self.without_deviants = delete_deviants
        self.error_limit = error_limit

        if self.without_deviants:
            print(f"Warning --- Values with Error larger than {self.error_limit}% will be neglected!")

        if self.__test_path(): 
            #--- load data
            self.__read_values() #if input is list of paths corresponding to xyz load them all
            #--- metrics are automatically calculated after assigning the profiles and thus norm_on_plateu can be done ! 
            #--- normalize X/Y Profiles on start
            if not self.normed_on_max:
                 self.norm_plateau()
            #--- output information on read
            if INFO:
                print(self)
        else:
            raise TypeError("Missing Input for Variable Path!" )
        print(" ")

    def __str__(self):
        """Print Information about the egslog dose object """
        return f"""
##############################################################################
  Successfully read the egslog-File(s) and assigned it to the class attributes\n
  
  X-Profile is {"not" if (self.position.x.size==0) else ""} available {"within the intervall " + str([self.position.x[0], self.position.x[-1]])+" !" if (self.position.x.size!=0) else "."}\n

  Y-Profile is {"not" if (self.position.y.size==0) else ""} available {"within the intervall " + str([self.position.y[0], self.position.y[-1]])+" !" if (self.position.y.size!=0) else "."}\n

  PDD is {"not" if (self.position.z.size==0) else ""} available {"within the intervall " + str([self.position.z[0], self.position.z[-1]])+" !" if (self.position.z.size!=0) else "."}\n
  
  NOTE THAT THE INFORMATION ABOUT THE PROFILE DEPTH AND (X,Y)-COORDIANTES FOR THE PDD NEED TO BE CONCLUDED FROM THE .egsinp manually !
  
  You can access the Dose-values quickly through OBJECT_NAME.pdd, OBJECT_NAME.x_profile, OBJECT_NAME.y_profile.
  Alternatnatively you can use set_profiles or set_pdd to change the position of each Doseprofile or use get_plane to get a 2D-np.array representing a plane.
  If you want to access individual Elements use the class attribute  'dose_matrix' an access it with dose_matrix[z_index, y_index, x_index].
  
  To save a profile one can just call the pandas Datafram class method pd.DataFrame(profile).to_csv()  e.g pd.DataFrame(dose_3d_object.y_profile).to_csv(destination_path)      
##############################################################################
""" 
    def __test_path(self):
        is_valid = True
        
        if type(self.origin) == list:
            self.has_multiple_inputs = True
            for input_string in self.origin:
                if not( input_string.endswith("egslog") if type(input_string )==str else False):
                    return False
        else:
            is_valid = self.origin.endswith("egslog") if type(self.origin)==str else False
            
        return is_valid 

    def __read_values(self):
        """ calls read_file to get the information and then orders it into a xyz array !"""
        self.position = _xyz_array_Object({"x":[], "y":[], "z":[]})
        if self.has_multiple_inputs:
            for file_path in self.origin:
                self.__read_file(file_path)
        else:
            self.__read_file(self.origin)


    def __read_file(self, file_path):
        """ return arrays of position, dose, and error for one file and a char that entails which axis was simulated"""
        with open(file_path) as file:
            content = file.readlines()
            
        #--- skip to start of value list
        for i, line in enumerate(content):
            if "this MUST be zero!" in line:
                i+=1
                break
        
        axis = content[i].split("_dose")[0] # in chamber geometries are called either x_dose, y_dose or z_dose ! 
        #--- skip to start of value list
        position, dose, error = [], [], []
        
        if self.without_deviants:
            for line in content[i::]:
                if line == "\n":
                    break
                pos_i, dose_i, err_i = np.array(line.split())[[0, 1, 3]] 
                if (float(err_i) if "%" not in err_i else float(err_i[:-1])) < self.error_limit:
                    position.append(float(pos_i.split("dose_")[1])) 
                    dose.append(float(dose_i)) 
                    error.append(float(err_i) if "%" not in err_i else float(err_i[:-1])) 
                
        else:
            for line in content[i::]:
                if line == "\n":
                    break
                pos_i, dose_i, err_i = np.array(line.split())[[0, 1, 3]] 
                position.append(float(pos_i.split("dose_")[1])) 
                dose.append(float(dose_i)) 
                error.append(float(err_i) if "%" not in err_i else float(err_i[:-1])) 
        
        #--- assign the read in values to the corresponding class attribute
        self.__assign_values(axis, position, dose, error)

    def __assign_values(self, axis, pos, dose, error):
        if  axis=="x":
            self.position.x = np.array(pos)
            self.x_profile = normable_array(dose)
            self.x_profile_error = np.array(error)
            self.x_profile.error = np.array(error)
            self.set_metrics_profile("X")
        elif axis=="y":
            self.position.y = np.array(pos)
            self.y_profile = normable_array(dose)
            self.y_profile_error = np.array(error)
            self.y_profile.error = np.array(error)
            self.set_metrics_profile("Y")
        else:
            self.position.z = np.array(pos)
            self.pdd = normable_array(dose)
            self.pdd_error = np.array(error)
            self.pdd.error = np.array(error)
            self.set_pdd_metrics()
 
    def remove_value(self, axis, position):
        """ """
        if axis.upper() == "Z":
            ix = self.find_closest_index(self.position.z, position)
            print(f"Warning --- removed z-Value at {self.position.z[ix]}cm!")
            self.position.z = np.delete(self.position.z, ix) 
            self.pdd = normable_array(np.delete(self.pdd, ix)) #rescale pdd.norm ! 
            self.pdd_error = np.delete(self.pdd_error, ix) 
            self.pdd.error = np.delete(self.pdd.error, ix) 
            self.set_pdd_metrics()
        
        if axis.upper() == "X":
            ix = self.find_closest_index(self.position.x, position)
            print(f"Warning --- removed x-Value at {self.position.x[ix]}cm!")
            self.position.x = np.delete(self.position.x, ix) 
            self.x_profile = normable_array(np.delete(self.x_profile, ix)) #rescale pdd.norm ! 
            self.x_profile_error = np.delete(self.x_profile_error, ix) 
            self.x_profile.error = np.delete(self.x_profile.error, ix) 
            self.set_metrics_profile("X")
            
        if axis.upper() == "Y":
            ix = self.find_closest_index(self.position.y, position)
            print(f"Warning --- removed y-Value at {self.position.y[ix]}cm!")
            self.position.y = np.delete(self.position.y, ix) 
            self.y_profile = normable_array(np.delete(self.y_profile, ix)) #rescale pdd.norm ! 
            self.y_profile_error = np.delete(self.y_profile_error, ix) 
            self.y_profile.error = np.delete(self.y_profile.error, ix) 
            self.set_metrics_profile("Y")
            

def plot_chamber_sim(path, err_lim=10, label=None):
    #-- read file
    with open(path) as file:
        content = file.readlines()
        
    #--- skip to start of value list
    for i, line in enumerate(content):
        if "this MUST be zero!" in line:
            i+=1
            break

    #--- skip to start of value list
    position, dose, error = [], [], []
    for line in content[i::]:
        if line == "\n":
            break
        pos_i, dose_i, err_i = np.array(line.split())[[0, 1, 3]] 
        position.append(float(pos_i.split("dose_")[1])) 
        dose.append(float(dose_i)) 
        error.append(float(err_i) if "%" not in err_i else float(err_i[:-1])) 
        
    position_wod, dose_wod, error_wod = [], [], []
    #--- filter deviants
    for i in range(len(dose)):
        if error[i] < err_lim:
            position_wod.append(position[i])
            dose_wod.append(dose[i])
            error_wod.append(error[i])
        
    #--- Plot everything
    fig, axs = plt.subplots(3, 2, figsize=(20,10), sharey="row", tight_layout=True)
    #--- with deviants
    axs[0,0].scatter(position, dose)
    axs[1,0].scatter(position, error)
    axs[2,0].hist(error, bins=list(np.linspace(0, 20, 100)))
    #--- without deviants
    axs[0,1].scatter(position_wod, dose_wod)
    axs[1,1].scatter(position_wod, error_wod)
    axs[2,1].hist(error_wod, bins=list(np.linspace(0, 20, 100)))
    if label:
        fig.suptitle(label, fontsize=25)

def plot_dose_with_gamma(dose_exp, dose_sim, axes, dx=0.3, dD=3, CUTOFF=20):
    color_exp   = "#3e423f"
    color_sim   = "#317fde"
    color_gamma = "#de3131"
    
    fig, axs = plt.subplots(1,1, figsize=(12,8), dpi=100)
    #--- make the gamma_index scale 
    ax_gamma = axs.twinx()  # instantiate a second Axes that shares the same x-axis    
    ax_gamma.set_ylabel(f'γ({axes.lower()})', color=color_gamma, fontsize=20)  # we already handled the x-label with ax1
    ax_gamma.tick_params(axis="y", labelcolor=color_gamma)
    ax_gamma.set_ylim(0, 5)
    axs.grid(False)
    ax_gamma.grid(False)
    ax_gamma.plot([-100, 100], [1,1], c=color_gamma, ls="dashed")
    min_x, max_x = 0, 0
    if axes.upper()=="X":
        axs.scatter(dose_exp.position.x, dose_exp.x_profile.norm, color=color_exp, s=50, label="Experimental")
        axs.plot(dose_sim.position.x, dose_sim.x_profile.norm, color=color_sim, lw=3, label="Simulated")
        gamma_x =gamma(
                        axes_reference = dose_exp.position.x, dose_reference  = dose_exp.x_profile.norm, 
                        axes_evaluation= dose_sim.position.x, dose_evaluation = dose_sim.x_profile.norm, 
                        dD=dD, dx=dx, CUTOFF=CUTOFF) #3% und 3mm

        ax_gamma.plot(dose_exp.position.x, gamma_x, c=color_gamma, lw=3)
        Gamma_PR_X = round(len(gamma_x[gamma_x<1])/len(gamma_x)*100,2)
        fig.suptitle(f"γ(x) Passing Rate = {Gamma_PR_X}%", font="monospace", fontsize=20, y=0.925)
        min_x, max_x = min(dose_exp.position.x), max(dose_exp.position.x)

        

    elif axes.upper()=="Y":
        axs.scatter(dose_exp.position.y, dose_exp.y_profile.norm, color=color_exp, s=50, label="Experimental")
        axs.plot(dose_sim.position.y, dose_sim.y_profile.norm, color=color_sim, lw=3, label="Simulated")
        gamma_y =gamma(
                        axes_reference = dose_exp.position.y, dose_reference  = dose_exp.y_profile.norm, 
                        axes_evaluation= dose_sim.position.y, dose_evaluation = dose_sim.y_profile.norm, 
                        dD=dD, dx=dx, CUTOFF=CUTOFF) #3% und 3mm

        ax_gamma.plot(dose_exp.position.y, gamma_y, c=color_gamma, lw=3)
        Gamma_PR_Y = round(len(gamma_y[gamma_y<1])/len(gamma_y)*100,2)
        fig.suptitle(f"γ(y) Passing Rate = {Gamma_PR_Y}%", font="monospace", fontsize=20, y=0.925)
        min_x, max_x = min(dose_exp.position.y), max(dose_exp.position.y)

    if axes.upper()=="Z":
        axs.scatter(dose_exp.position.z, dose_exp.pdd.norm, color=color_exp, s=50, label="Experimental")
        axs.plot(dose_sim.position.z,    dose_sim.pdd.norm, color=color_sim, lw=3, label="Simulated")
        gamma_z =gamma(
                        axes_reference = dose_exp.position.z, dose_reference  = dose_exp.pdd.norm, 
                        axes_evaluation= dose_sim.position.z, dose_evaluation = dose_sim.pdd.norm, 
                        dD=dD, dx=dx, CUTOFF=CUTOFF) #3% und 3mm

        ax_gamma.plot(dose_exp.position.z, gamma_z, c=color_gamma, lw=3)
        Gamma_PR_Z = round(len(gamma_z[gamma_z<1])/len(gamma_z)*100,2)
        fig.suptitle(f"γ(z) Passing Rate = {Gamma_PR_Z}%", font="monospace", fontsize=20, y=0.925)
        min_x, max_x = 0, 0
        min_x, max_x = min(dose_exp.position.z), max(dose_exp.position.z)
        
    #--- fix labels so you dont see the dashed line
    axs.set_xlim(min_x, max_x)
    #--- labels and other cosmetics
    axs.set_xlabel(f"Position along {axes}-Axis [cm]", fontsize=15)
    axs.set_ylabel(f"relative Dose [%]", fontsize=15)
    #---
    axs.legend(loc="best", fontsize=12)
    
    return fig, axs 