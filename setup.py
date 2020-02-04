from setuptools import setup
from setuptools import find_packages


version_py = "fithic/_version.py"
exec(open(version_py).read())

setup(
    name = "fithic",
    version =__version__,
    description = 'Hi-C Analysis software created and maintained by the Ay Lab',
    long_description = "Fit-Hi-C is a Hi-C software analysis tool to compute statistical confidence estimates to Hi-C contact maps. It is good. Very good.",
    url = 'http://github.com/ay-lab/fithic/',
    entry_points = {
        "console_scripts": ['fithic = fithic.fithic:main']
        },
    python_requires = '~=3.6',
    author = 'Ferhat Ay',
    author_email = 'ferhatay@lji.org',
    license = 'MIT',
    packages = ['fithic'],
    install_requires = [
        'numpy',
        'matplotlib',
        'scipy',
        'scikit-learn',
        'sortedcontainers',
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    zip_safe = False,
  )
