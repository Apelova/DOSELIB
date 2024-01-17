"""
@author: apel04
email: marvin.apel.ext@ptb.de

path = "/home/marvin/Desktop/master"     # Notebook
path = "C:/Users/apel04/Desktop/master/" # Office Pc
"""
import matplotlib.pyplot as plt
from scipy import interpolate
from doselib import *
import pandas as pd
import numpy as np
import matplotlib

path = "C:/Users/apel04/Desktop/master"  # Office Pc
path = "/home/marvin/Desktop/master"    # Notebook

dose_dosy = dose_3d(
    path+"/Simulationen/DATA/TEST_PTB/DOSXYZ_NRC/ROUGH_GRID/PTB_6MV_old_dosxyz.3ddose")
dose_dosy.set_x_profile(Z=5.25, Y=0)
dose_dosy.set_y_profile(Z=5.25, X=0)

# --- chamber-3ddose-File!
dose_chamber = dose_3d(
    path+"/Simulationen/DATA/TEST_PTB/EGS_CHAMBER/ROUGH_GRID/PTB_6MeVp_10x10_old_dose_pdd.3ddose")
dose_chamber.add_profile(
    path+"/Simulationen/DATA/TEST_PTB/EGS_CHAMBER/ROUGH_GRID/PTB_6MeVp_10x10_old_dose_x_profile.3ddose", AXIS="X")
dose_chamber.add_profile(
    path+"/Simulationen/DATA/TEST_PTB/EGS_CHAMBER/ROUGH_GRID/PTB_6MeVp_10x10_old_dose_y_profile.3ddose", AXIS="Y")

# --- experimentell-Data!
#dose_exp = dose_mcc(path+"/Messungen/6MeV_10x10_Dose_Profiles/profiles_@_5cm_ff_and_fff_SSD_100/FF/combined_ff.mcc")
dose_exp = dose_mcc(
    path+"/Messungen/6MeV_10x10_Dose_Profiles/profiles_@_5_25cm_ff_SSD_98/all_profiles_scan.mcc")

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

# %%
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
            metric_string += f"{rename_attr(attr):<{dL}}: "
            metric_string += add_digits_to_string(f"{getattr(DOSE_OBJ, attr):.2f}") 
            metric_string += f"{UNITS[attr]}\n"
        else:
            metric_string += f"{rename_attr(attr):<{dL}}: "
            metric_string += add_digits_to_string(f"{getattr(DOSE_OBJ, attr):.2f}")
            metric_string += f"{UNITS[attr]}\n"
    # --- remove last "\n"
    return metric_string[:-1]

def set_metrics(dose_objs, plot_axs, axes=["Z", "X", "Y"], difference=True, exclude= {}):
    dh_x = 6-len(exclude["X"]) if "X" in exclude else 6
    dh_y = 6-len(exclude["Y"]) if "Y" in exclude else 6
    # Das hier umschreiben !!! so die reihenfolge mit axes stimmt !

    for  i, dose in enumerate(dose_objs[:2]):#compare maximum two dose objs !
        for j, ax in enumerate(axes):
            if difference:
                # Match statements require python 3,8 -> not available in the office therefore use this shananigans
                if ax == "Z":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="Z", exclude=exclude)
                    axs[0,j].text(0.69+i*0.29, 0.94, metric_string, transform=axs[0,j].transAxes, fontsize=19, bbox=props, ha="right", va="top", font="monospace")
                elif ax =="X":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="X", exclude=exclude)
                    axs[0,j].text(0.98, 0.94-i*0.05*dh_x, metric_string, transform=axs[0,j].transAxes, fontsize=19, bbox=props, ha="right", va="top", font="monospace")
                elif ax =="Y":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="Y", exclude=exclude)
                    axs[0,j].text(0.98, 0.94-i*0.05*dh_y, metric_string, transform=axs[0,j].transAxes, fontsize=19, bbox=props, ha="right", va="top", font="monospace")
            else:
                if ax =="Z":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="Z", exclude=exclude)
                    axs[j].text(0.76+i*0.22, 0.94, metric_string, transform=axs[j].transAxes, fontsize=15, bbox=props, ha="right", va="top", font="monospace")
                elif ax =="X":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="X", exclude=exclude)
                    axs[j].text(0.98, 0.94-i*0.06*(dh_x+1), metric_string, transform=axs[j].transAxes, fontsize=15, bbox=props, ha="right", va="top", font="monospace")
                elif ax =="Y":
                    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
                    metric_string = get_metric_string(dose, AXIS="Y", exclude=exclude)
                    axs[j].text(0.98, 0.94-i*0.06*(dh_y+1), metric_string, transform=axs[j].transAxes, fontsize=15, bbox=props, ha="right", va="top", font="monospace")
     
# TODO in the office
# Incorporate in main
# Test if it works in main !

exludio = {"Z":["D20", "D10"], "X":["D0_x", "dCAX_x"], "Y":[]}        
axerinios = ["Z", "X", "Y"]
differencio = True
colors_ = colors=["black", "blue"]
fig, axs = compare_dose([dose_exp, dose_chamber], labels=["Experimental","EGS_CHAMBER"], axes=axerinios, difference=differencio, interpol="quadratic", diff_dx=0.1, colors=colors_)
set_metrics([dose_exp, dose_chamber], axs, axes=axerinios, difference=differencio, exclude = exludio)








