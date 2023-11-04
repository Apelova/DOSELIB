from setuptools import find_packages, setup

setup(
    name='doselib',
    version='0.2',
    author='Marvin Apel',
    url = "https://github.com/Apelova/EGS_DOSE_TOOLS",
    packages=find_packages(include = ["doselib"]),
    description='Library to extract Data or freatures from EGSnrc Monte-Carlo Simulation results',
    maintainer_email = "Marvin.Apel.99@gmail.com",
    setup_requires=['pandas']
)