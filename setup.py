from setuptools import setup
from setuptools import find_packages

setup(
    name = "fithic",
    version = '2.0.4',
    description = 'Hi-C Analysis software created and maintained by the Ay Lab',
    url = 'http://github.com/ay-lab/fithic/tree/development',
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
