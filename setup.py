"""
setup.py: python package setup for ITSxpress
"""
from setuptools.command.install import install
from setuptools.command.develop import develop
from setuptools import setup
import os

def custom_command():
    import sys
    cmd  = 'wget https://sourceforge.net/projects/bbmap/files/BBMap_38.60.tar.gz -O /usr/BBMap_38.60.tar.gz; tar -xvf /usr/BBMap_38.60.tar.gz -C /usr'

    cmd1 = 'wget https://github.com/torognes/vsearch/releases/download/v2.13.6/vsearch-2.13.6-linux-x86_64.tar.gz -O /usr/vsearch-2.13.6-linux-x86_64.tar.gz; tar -xvf /usr/vsearch-2.13.6-linux-x86_64.tar.gz -C /usr'

    cmd2 = 'wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz -O /usr/hmmer-3.1b2-linux-intel-x86_64.tar.gz; tar -xvf /usr/hmmer-3.1b2-linux-intel-x86_64.tar.gz -C /usr'

    if sys.platform in ['darwin', 'linux']:
        os.system(cmd)
        os.system(cmd1)
        os.system(cmd2)
        

class CustomInstallCommand(install):
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
    cmdclass={'install': CustomInstallCommand},
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
