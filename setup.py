import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name = "synthedia",
    version = "0.0.1",
    description = "Generate synthetic DIA LC-MS/MS proteomics data",
    long_description = README,
    long_description_content_type = "text/markdown",
    url = "https://github.com/mgleeming/synthedia.git",
    author = "Michael Leeming",
    author_email = "leemingm@unimelb.edu.au",
    license = "BSD",
    classifiers = [
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages = ["synthedia"],
    scripts = ['synthedia/synthedia.py'],
    install_requires = [
        "pandas", "numpy", "pyopenms",
        "pyteomics", "numba", "matplotlib"
    ],
)

