from setuptools import setup

setup(
    name='itsxpress',
    version='1.5.2',
    packages=['itsxpress'],
    license='License :: CC0 1.0 Universal (CC0 1.0) Public Domain Dedication',
    description=open('README.rst').read(),
    long_description="Rapidly trim sequences down to their Internally Transcribed Spacer (ITS) regions",
    classifiers=['Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.5',
                 'Development Status :: 3 - Alpha'],
    keywords='Amplicon sequencing fungal ITS',
    url='http://github.com/usda-ars-gbru/itsxpress',
    test_suite ='nose.collector',
    author='Adam R. Rivers',
    author_email='adam.rivers@ars.usda.gov',
    install_requires=['biopython>=1.60'],
    python_requires='>3.5',
    tests_require=['nose'],
    include_package_data=True,
    entry_points= {'console_scripts':['itsxpress=itsxpress.main:main']},
    zip_safe=False)
