# Polcomp: compartment model as a command-line interface

A command line utility to estimate parameters for Compartment Model. For more details, please refer to documentation and Dutta et al, Gen Res, 2023

Installation steps:
- cd /path/to/yourdirectory
- python3 -m venv
- source venv/activate
- python3 -m pip install compartmentModel

Quick start/user guide: 

Vignette: https://guertinlab.github.io/compartment_model/compModel_vignette/compModel_vignette/compartmentModel_vignette.html

Package: https://pypi.org/project/compartmentModel/
 

Package building:
- python3 -m pip install build 
- python3 -m build
- python3 -m pip install twine
- python3 -m twine check dist/* 
- python3 -m twine upload -r pypi dist/*

