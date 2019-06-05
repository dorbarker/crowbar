from crowbar import __version__
from setuptools import setup, find_packages

setup(
    name='crowbar',
    version=__version__,
    packages=find_packages(),
    install_requires=[
        'pandas>=0.22.0',
        'numpy>=1.16.3',
        'biopython>=1.73'
        ],
    author='Dillon Barker',
    author_email='dillon.barker@canada.ca',

    entry_points={
        'console_scripts': [
            'crowbar=crowbar.main:main',
            'crowbar-simulate=crowbar.test.simulate_recovery:main'
        ]
    }
)
