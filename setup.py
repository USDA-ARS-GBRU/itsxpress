"""setup.py: python package setup for ITSxpress
"""

from setuptools import setup
import subprocess

if __name__ == "__main__":
    subprocess.run("conda install -y -c conda-forge hmmer==3.1b2",shell=True)
    subprocess.run("conda install -y -c bioconda bbmap==38.69",shell=True)
    subprocess.run("conda install -y -c bioconda vsearch==2.21.1",shell=True)
    setup()
