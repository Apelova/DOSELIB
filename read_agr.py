"""
@author: Apelova 

Module Task: get Data from .agr-File as output from statdose.
"""
import matplotlib.pyplot as plt
import pandas as pd


def get_data_from_agr(path, norm_dose = True):
    if not path.endswith(".agr"):
        raise ValueError("Can't extract Data from this File. This routine can only read .agr Files from statdose.")
    
    #--- initialize Placeholder for Output
    data_dict = {}
    #--- used to test if the file is empty and raises an exception
    is_empty = True
    #--- This Dictionary Enables the correct Display of coordinate axis 
    axs = { "X": ["Y", "Z"],
               "Y": ["X", "Z"],
               "Z": ["X", "Y"],
             }
    #--- read and parse file from path
    with open(path) as file:
        line = file.readline()
        while line:
            #--- find lines at which meassurement start and extract each data in one data_set
            if "/cm" in line: 
                is_empty = False
                #--- get Information about which direction the Plot shows
                axis = str(line.split()[3][1])
                while line:
                    coordinates  = axis + "-Axis" 
                    if "legend string" in line: 
                        #--- initialize Dataframe that stores Data
                        temp_data = pd.DataFrame({
                                "position": [],
                                "dose":     [],
                                "error":    [],
                                })
                        #--- Get Information about start coordinated -> this Mumbo Jumbo fixes an inconsistency in the legend string of .agr 
                        vals = line.split("(")[-1].split(",")
                        vals[1] = vals[1].split("\n")
                        vals[1] = str(float(vals[1][0].split(")")[0])) if (")" in vals[1][0]) else vals[1][0].split('"')[0]                        
                        coordinates += f" {axs[axis][0]} = {vals[0]}  {axs[axis][1]} = {vals[1]} "
                        #--- skip lines until start of values:
                        file.readline()
                        file.readline()
                        file.readline()    
                        #--- extract the data converted to floats from ASCII-chars
                        line = file.readline()
                        while line != "&\n":
                            temp_data.loc[len(temp_data)] = [float(val) for val in line.split()]
                            line = file.readline()
                        #--- add curve_data to Dictionary and go to next curve
                        #--- get rid od NaN-Entries
                        data_dict[coordinates] = temp_data.fillna(0)
                    #---continue inner while loop
                    line = file.readline()
            #---continue outer while loop
            line = file.readline()
        
        #--- raise an Excpetion when file is empty:
        if is_empty:
            raise AttributeError("Files does not contain a valid Datasequence.")
        
    #--- norm dose if norm_dose
    if norm_dose:
        dose_max = max([ max(data_dict[curve].dose) for curve in data_dict])
        for curve in data_dict:
            data_dict[curve].dose = data_dict[curve].dose/dose_max*100
            data_dict[curve].error = data_dict[curve].error/dose_max*100 
    #--- return extracted data
    return data_dict


if __name__== "__main__":
    #--- Example 1 ()
    data = get_data_from_agr("EGSnrc/EGSnrc/egs_home/dosxyznrc/tdk_16MVp.agr")            
    for meas in data:
        plt.plot(data[meas].position, data[meas].dose, label=meas)
    plt.legend(loc="upper right")
    plt.show()
    #%%
    #--- Example 2
    data = get_data_from_agr("EGSnrc/EGSnrc/egs_home/dosxyznrc/tdk_16MVp_quer.agr")
    for meas in data:
        plt.plot(data[meas].position, data[meas].dose, label=meas)
    plt.legend(loc="upper right")
    plt.show()    
    #%%
    #--- Example 3
    data = get_data_from_agr("EGSnrc/EGSnrc/egs_home/dosxyznrc/pdd_z_diff_start.agr")            
    for meas in data:
        plt.plot(data[meas].position, data[meas].dose, label=meas)
    plt.legend(loc="upper right")
    plt.show()
    


    
    
    
    





















