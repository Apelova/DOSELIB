from setuptools import find_packages, setup

setup(
    name='doselib',
    version='1.0',
    author='Marvin Apel aka Apelova',
    url = "https://github.com/Apelova/EGS_DOSE_TOOLS",
    description='Library to extract Data or freatures from EGSnrc Monte-Carlo Simulation results',
    maintainer_email = "Marvin.Apel.99@gmail.com",
    setup_requires=['pandas'],
    license='MIT',
    packages=['egs_doselib']
)
