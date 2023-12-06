# Doselib - A library for processing EGSnrc's Monte-Carlo Simulation output 

# 1. Install the package through the -whl-File by running:
	~ pip install /path/to/wheelfile.whl

# 2. Import the package by using either 
	from egs_doselib import * 
	or
	import egs_doselib as dl
	
# 3. Available Classes are:
 	- dose_3d (equivalent to statdose)
		- a tool which reads 3ddose files and enables quick and simple extraction of doseprofiles (x/y-profile and pdd) and slices (parallel to xy, xz or yz)
	- dose_mcc 
		 - a class with routines to read the .mcc files output from Mephysto

# 4. In The file Example.py the most important features of this package are shown.
