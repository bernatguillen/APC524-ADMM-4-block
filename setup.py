from distutils.core import setup

setup(
    name='admm4block',
    version='0.1dev',
    packages=['admm4block','admm4block.conic','admm4block.sdp'],
    license='GNU GENERAL PUBLIC LICENSE',
    long_description=open('README.md').read(),
)