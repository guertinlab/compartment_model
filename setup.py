from setuptools import setup

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
     name = 'compartmentModel',
     version='1.01',
     license='MIT',
     description='A command line utility to estimate parameters for Compartment Model. For more details see documentation and Dutta et al, Gen Res, 2023.',
     long_description=long_description,
     long_description_content_type='text/markdown',
     author='Lab of Michael J. Guertin',
     author_email='guertin@uchc.edu',
     url = 'https://github.com/guertinlab',
     py_modules=['compmodel.compmodel'],
     install_requires=['Click', 
                'pandas', 'configparser', 'PyYAML', 'tqdm', 'numpy','scipy','matplotlib'],
     entry_points='''
         [console_scripts]
     polcomp=compmodel.compmodel:cli
     '''
     )

