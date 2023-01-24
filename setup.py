"""setup.py: python package setup for ITSxpress
"""

from setuptools import setup
import subprocess

if __name__ == "__main__":
    subprocess.run("conda install -y -c conda-forge hmmer==3.1b2",shell=True)
    subprocess.run("conda install -y -c bioconda bbmap==38.69",shell=True)
    subprocess.run("conda install -y -c bioconda vsearch==2.21.1",shell=True)
    #Check dependency versions and assert


    
    outp_hmmer = str(subprocess.run("conda list hmmer",shell=True,stdout=subprocess.PIPE))
    outp_bbmap = str(subprocess.run("conda list bbmap",shell=True,stdout=subprocess.PIPE))
    outp_vsearch = str(subprocess.run("conda list vsearch",shell=True,stdout=subprocess.PIPE))
    assert (outp_hmmer.find('3.1b2') != -1),"hmmer doesn't seem to be version 3.1b2, create a fresh conda environment and install itsxpress again"
    assert (outp_bbmap.find('38.69') != -1),"bbmap doesn't seem to be version 38.69, create a fresh conda environment and install itsxpress again"
    assert (outp_vsearch.find('2.21.1') != -1),"vsearch doesn't seem to be version 2.21.1, create a fresh conda environment and install itsxpress"
    setup()
