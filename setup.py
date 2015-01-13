from setuptools import setup,find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here,'README.txt'), encoding = 'utf-8') as f:
    long_description = f.read()
    
setup(
    name='admm4block',
    
    version='0.3dev',
    description='ADMM4block solver for DNNSDP',
    long_description = long_description,
    
    url = 'https://github.com/bernatguillen/APC524-ADMM-4-block',
    
    author = 'Bernat Guillen, Michael Tarczon, Yuan Liu',
    author_email = 'bernatp@princeton.edu',
    license='GNU GENERAL PUBLIC LICENSE',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.3',
        'Programming Language :: Python :: 2.4',
        'Programming Language :: Python :: 2.5',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Mathematics'
        ],
    keywords = 'Convex Optimization SDP',
    test_suite =  'nose.collector',
    tests_require = ['nose'],
    packages=find_packages(exclude = ['contrib','docs','tests*']),
    scripts=['bin/points.py'],
    install_requires=['numpy'],

)