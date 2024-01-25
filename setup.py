#
# moonshot_project setuptools script
#
from setuptools import setup, find_packages


def get_version():
    """
    Get version number from the moonshot_project module.
    """
    import os
    import sys

    sys.path.append(os.path.abspath('moonshot_project'))
    version = "1.0.0"
    sys.path.pop()

    return version


def get_requirements():
    requirements = []
    with open("requirements.txt", "r") as file:
        for line in file:
            requirements.append(line)
    return requirements


setup(
    # Module name
    name='moonshot_project',

    # Version
    version=get_version(),

    description='Machine Learning model for the Covid Moonshot Project.',

    maintainer='Matthew Ghosh, Ronald Cvek',

    maintainer_email='matthew.ghosh@gtc.ox.ac.uk, ronald.cvek@stcatz.ox.ac.uk',

    url='https://github.com/mghosh00/MoonshotProject',

    # Packages to include
    packages=find_packages(include=('moonshot_project', 'moonshot_project.*')),

    # List of dependencies
    install_requires=get_requirements(),

    extras_require={
        'dev': [
            # Flake8 for code style checking
            'flake8>=3',
            'pytest',
            'pytest-cov',
        ],
    },
)
