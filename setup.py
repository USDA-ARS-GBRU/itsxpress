"""
setup.py: python package setup for ITSxpress
"""
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools import setup
import os
import site
import sys
import subprocess

def custom_command():
    package_path = site.getsitepackages()[0]
    print(site.getsitepackages())
    
    # cmd  = "wget https://sourceforge.net/projects/bbmap/files/BBMap_38.60.tar.gz -O {package_path}/itsxpress/BBMap_38.60.tar.gz; tar -xvf {package_path}/itsxpress/BBMap_38.60.tar.gz -C {package_path}/itsxpress/".format(package_path=package_path)

    # cmd1 = "wget https://github.com/torognes/vsearch/releases/download/v2.13.6/vsearch-2.13.6-linux-x86_64.tar.gz -O {package_path}/itsxpress/vsearch-2.13.6-linux-x86_64.tar.gz; tar -xvf {package_path}/vsearch-2.13.6-linux-x86_64.tar.gz -C {package_path}/itsxpress//".format(package_path=package_path)

    # cmd2 = "wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz -O {package_path}/itsxpress/hmmer-3.1b2-linux-intel-x86_64.tar.gz; tar -xvf {package_path}/hmmer-3.1b2-linux-intel-x86_64.tar.gz -C {package_path}/itsxpress/".format(package_path=package_path)

    if sys.platform in ['darwin', 'linux']:
        parameters = ["mkdir",
                    "{package_path}/itsxpress".format(package_path=package_path)]
        subprocess.run(parameters)

        parameters = ["wget", 
                    "https://sourceforge.net/projects/bbmap/files/BBMap_38.60.tar.gz", 
                    "-O", "{package_path}/itsxpress/BBMap_38.60.tar.gz".format(package_path=package_path)]
        subprocess.run(parameters)
        parameters = ["tar", 
                    "-xvf" "{package_path}/itsxpress/BBMap_38.60.tar.gz".format(package_path=package_path), 
                    "-C",
                    "{package_path}/itsxpress/".format(package_path=package_path)]
        subprocess.run(parameters)

        parameters = ["wget", 
                    "https://github.com/torognes/vsearch/releases/download/v2.13.6/vsearch-2.13.6-linux-x86_64.tar.gz", 
                    "-O", "{package_path}/itsxpress/vsearch-2.13.6-linux-x86_64.tar.gz".format(package_path=package_path)]
        subprocess.run(parameters)
        parameters = ["tar", 
                    "-xvf" "{package_path}/itsxpress/vsearch-2.13.6-linux-x86_64.tar.gz".format(package_path=package_path), 
                    "-C",
                    "{package_path}/itsxpress/".format(package_path=package_path)]
        subprocess.run(parameters)

        parameters = ["wget", 
                    "http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz", 
                    "-O", "{package_path}/itsxpress/hmmer-3.1b2-linux-intel-x86_64.tar.gz".format(package_path=package_path)]
        subprocess.run(parameters)
        parameters = ["tar", 
                    "-xvf" "{package_path}/itsxpress/hmmer-3.1b2-linux-intel-x86_64.tar.gz".format(package_path=package_path), 
                    "-C",
                    "{package_path}/itsxpress/".format(package_path=package_path)]
        subprocess.run(parameters)
        

class PostInstallCommand(install):
    def run(self):
        install.run(self)
        print("custom install")     
        custom_command()

# def _post_install(setup):
#     def _post_actions():
#         custom_command()
#     _post_actions()
#     return setup



setup(
    name='itsxpress',
    version='1.8.0',
    packages=['itsxpress'],
    license='License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
    description="Rapidly trim sequences down to their Internally Transcribed Spacer (ITS) regions",
    long_description=open('README.rst').read(),
    classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.5',
                 'Development Status :: 3 - Alpha'],
    keywords='Amplicon sequencing fungal ITS',
    url='http://github.com/usda-ars-gbru/itsxpress',
    cmdclass={'install': PostInstallCommand},
    test_suite='nose.collector',
    author='Adam R. Rivers',
    author_email='adam.rivers@ars.usda.gov',
    install_requires=['biopython>=1.60'],
    python_requires='>3.5',
    tests_require=['nose'],
    include_package_data=True,
    entry_points={'console_scripts':['itsxpress=itsxpress.main:main'],
    'qiime2.plugins': ['itsxpress=itsxpress.plugin_setup:plugin']},
    zip_safe=False)


# custom_command()