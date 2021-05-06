from setuptools import setup, find_packages

version = {}
with open("ediacara/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="ediacara",
    version=version["__version__"],
    author="Peter Vegh",
    description="EGF's alignment comparator",
    long_description=open("pypi-readme.rst").read(),
    long_description_content_type="text/x-rst",
    keywords="biology",
    packages=find_packages(exclude="docs"),
    install_requires=[
        "matplotlib",
        "numpy",
        "pandas",
        "biopython",
        "dna_features_viewer",
    ],
)
