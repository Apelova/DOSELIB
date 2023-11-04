# Doselib - A library for processing EGSnrc's Monte-Carlo Simulation output 

Thanks  to KIA EISINGA's blog this repository can be turned into a python module by following a few steps
LINK: https://medium.com/analytics-vidhya/how-to-create-a-python-library-7d5aea80cc3f

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

