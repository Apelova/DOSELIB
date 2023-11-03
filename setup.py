from setuptools import find_packages, setup

setup(
    name='doselib',
    packages=find_packages(include=["doselib"]),
    version='0.1',
    description='Library to extract Data or freatures from EGSnrc Monte-Carlo Simulation results',
    author='Apelova',
)