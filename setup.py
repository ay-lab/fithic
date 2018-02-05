from setuptools import setup

setup(
    name = "fithic",
    version = '1.1.2',
    description = 'Hi-C Analysis software created and maintained by the Ay Lab',
    url = 'http://github.com/ay-lab/fithic',
    entry_points = {
        "console_scripts": ['fithic = fithic.fithic:main']
        },
    python_requires = '~=2.7',
    author = 'Ferhat Ay',
    author_email = 'ferhatay@lji.org',
    license = 'MIT',
    packages = ['fithic'],
    install_requires = [
        'numpy',
        'matplotlib',
        'scipy',
        'scikit-learn',
    ],
    zip_safe = False,
  )
