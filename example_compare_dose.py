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
        "D0": "%",
        "Rmax": "cm",
        "D10": "%",
        "D20": "%",
        "Q": " ",
        "D0_x": "%",
        "penumbra_xl": "cm",
        "penumbra_xr": "cm",
        "H_x": " ",
        "S_x": " ",
        "dCAX_x": "cm",
        "D0_y": "%",
        "penumbra_yl": "cm",
        "penumbra_yr": "cm",
        "H_y": " ",
        "S_y": " ",
        "dCAX_y": "cm" } 

    if AXIS in exclude:
        ATTRIBUTES[AXIS] = [attr for attr in ATTRIBUTES[AXIS] if attr not in exclude[AXIS]] 

    # --- gather requested  information in string and output
    metric_string = ""
    for attr in ATTRIBUTES[AXIS]:
        if AXIS.upper() == "Z":
            metric_string += f"{rename_attr(attr):4}:{getattr(DOSE_OBJ, attr): 7.3f}{UNITS[attr]}\n"
            print(metric_string )
            #metric_string += f"{rename_attr(attr)}:{getattr(DOSE_OBJ, attr): .3f}{UNITS[attr]}\n"
        else:
            metric_string += f"{rename_attr(attr):7}:{round(getattr(DOSE_OBJ, attr),3): 3}{UNITS[attr]}\n"
            #metric_string += f"{rename_attr(attr)}:{getattr(DOSE_OBJ, attr): >.3f}{UNITS[attr]}\n"
    # --- remove last "\n"
    return metric_string[:-1]


def set_metrics(dose_objs, plot_axs, difference=True, exclude= {}):
    dh_x = 6-len(exclude["X"]) if "X" in exclude else 6
    dh_y = 6-len(exclude["Y"]) if "Y" in exclude else 6
    for  i, dose in enumerate(dose_objs[:2]):#compare maximum two dose objs !
        if difference:
            # Z-Axis
            props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
            metric_string = get_metric_string(dose, AXIS="Z", exclude=exclude)
            axs[0,0].text(0.7+i*0.27, 0.94, metric_string, transform=axs[0,0].transAxes, fontsize=20, bbox=props, ha="right", va="top")
            # X-Axis
            props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
            metric_string = get_metric_string(dose, AXIS="X", exclude=exclude)
            axs[0,1].text(0.98, 0.94-i*0.05*dh_x, metric_string, transform=axs[0,1].transAxes, fontsize=20, bbox=props, ha="right", va="top")
            # Y-Axis
            props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
            metric_string = get_metric_string(dose, AXIS="Y", exclude=exclude)
            axs[0,2].text(0.98, 0.94-i*0.05*dh_y, metric_string, transform=axs[0,2].transAxes, fontsize=20, bbox=props, ha="right", va="top")
        else:
            # Z-Axis
            props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
            metric_string = get_metric_string(dose, AXIS="Z", exclude=exclude)
            axs[0].text(0.45+i*0.28, 0.77, metric_string, transform=axs[0,0].transAxes, fontsize=20, bbox=props, ha="right", va="top")
            # X-Axis
            props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
            metric_string = get_metric_string(dose, AXIS="X", exclude=exclude)
            axs[1].text(0.98, 0.95-i*0.29, metric_string, transform=axs[0,1].transAxes, fontsize=20, bbox=props, ha="right", va="top")
            # Y-Axis
            props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
            metric_string = get_metric_string(dose, AXIS="Y", exclude=exclude)
            axs[2].text(0.98, 0.95-i*0.29, metric_string, transform=axs[0,2].transAxes, fontsize=20, bbox=props, ha="right", va="top")


colors_ = colors=["black", "blue"]
fig, axs = compare_dose([dose_exp, dose_chamber], labels=["Experimental","EGS_CHAMBER"], axes=[
                        "Z", "X", "Y"], difference=True, interpol="quadratic", diff_dx=0.1, colors=colors_)

set_metrics([dose_exp, dose_chamber], axs, difference=True, exclude = {"X": ["penumbra_xl", "penumbra_xr"], "Y": ["penumbra_yl", "penumbra_yr"], "Z": []})


# %%
#--- axis setting
axs[0, 0].legend(loc="lower left", fontsize=30)
axs[0, 0].set_xlim(0, 30)
axs[0, 1].set_xlim(-15, 15)
axs[0, 2].set_xlim(-15, 15)
# first diff
axs[1, 0].set_ylim(-10, 10)
axs[1, 0].set_yticks(np.arange(-10, 15, 5))
# second diff
axs[1, 1].set_ylim(-10, 10)
axs[1, 1].set_yticks(np.arange(-10, 15, 5))
# third diff
axs[1, 2].set_ylim(-20, 20)
axs[1, 2].set_yticks(np.arange(-20, 25, 10))

# %%
colors_ = colors=["black", "blue"]
fig, axs = compare_dose([dose_exp, dose_chamber], labels=["Experimental","EGS_CHAMBER"], axes=[
                        "X", "Y"], difference=True, interpol="quadratic", diff_dx=0.1, colors=colors_)
#--- axis setting
# first diff
# second diff
axs[1, 0].set_ylim(-10, 10)
axs[1, 0].set_yticks(np.arange(-10, 15, 5))
# third diff
axs[1, 1].set_ylim(-20, 20)
axs[1, 1].set_yticks(np.arange(-20, 25, 10))

# TODO!
#->auswählen der metriken
#-> in doselib integrieren -> für 2 dose_objetcs die funktio callen -> setzt figsize automatich  (avoids dynamic design)

# --- Metrics
for  i, dose in enumerate([dose_exp, dose_chamber]):    
    # X-Axis
    props = dict(boxstyle='square', facecolor=colors[i], edgecolor="black", alpha=0.25, lw=3)
    metric_string = get_metric_string(dose, AXIS="X")
    axs[0,0].text(0.98, 0.95-i*0.3, metric_string, transform=axs[0,0].transAxes, fontsize=20, bbox=props, ha="right", va="top")
    # Y-Axis
    props = dict(boxstyle='square', facecolor=colors[1], edgecolor="black", alpha=1, lw=0.25)
    metric_string = get_metric_string(dose, AXIS="Y", )
    axs[0,1].text(0.98, 0.95-i*0.3, metric_string, transform=axs[0,1].transAxes, fontsize=20, bbox=props, ha="right", va="top")



