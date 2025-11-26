# Polcomp: compartment model as a command-line interface

A command line utility to estimate parameters for Compartment Model. For more details, please refer to documentation and Dutta et al, Gen Res, 2023

Installation steps (virtual environment `venv` not needed but recommended):
- cd /path/to/your/directory
- python3 -m venv myenv ## create your environment
- source myenv/bin/activate ## activate your environment
- python3 -m pip install compartmentModel ## can directly install this without venv steps above

Quick start/user guide: 

Vignette: https://guertinlab.github.io/compartment_model/compModel_vignette/compModel_vignette/compartmentModel_vignette.html

Package: https://pypi.org/project/compartmentModel/
 

Package building (in case you are downloading/cloning the entire package source files):
- python3 -m pip install build 
- python3 -m build
- python3 -m pip install twine
- python3 -m twine check dist/* 
- python3 -m twine upload -r pypi dist/*

