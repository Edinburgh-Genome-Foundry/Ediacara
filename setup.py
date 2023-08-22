from setuptools import setup, find_packages

version = {}
with open("ediacara/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="ediacara",
    version=version["__version__"],
    author="Peter Vegh",
    description="EGF's alignment reporter",
    long_description=open("pypi-readme.rst").read(),
    long_description_content_type="text/x-rst",
    license="GPLv3+",
    keywords="biology, sequencing",
    packages=find_packages(exclude="docs"),
    include_package_data=True,
    install_requires=[
        "matplotlib",
        "numpy",
        "pandas",
        "weighted_levenshtein",
        "biopython",
        "dna_features_viewer",
        "geneblocks",
        "cyvcf2",
        "pdf_reports",
        "portion",
    ],
)
