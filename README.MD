# Doselib - A library for processing EGSnrc's Monte-Carlo Simulation output 

For now this repository is just a remote storage for some code i wrote.
It will soon be ready to be used by others.

# 1. Test if required packages are installed:
	~ pip install wheel
	~ pip install setuptools
	~ pip install twine

# 2. Compile the build distribution (.whl) only neccessary if the file in /dist is corrupted 
 	~ python setup.py bdist_wheel

# 3. Install the package through the -whl-File by running:
	~ pip install /path/to/wheelfile.whl

# 4. Import the package by using 
	import doselib.functions as dl
	
# 5. Available methods are:
 	- get_data_from_agr, which reads the data ploted in the statdose method from the corresponding .agr-File 

