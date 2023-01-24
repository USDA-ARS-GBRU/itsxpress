"""setup.py: python package setup for ITSxpress

"""

from setuptools import setup
import subprocess

if __name__ == "__main__":
    subprocess.run("conda install -y -c conda-forge hmmer",shell=True)
    subprocess.run("conda install -y -c bioconda vsearch",shell=True)
    setup()

