from setuptools import setup, find_packages

version = {}
with open("seqreport/version.py") as fp:
    exec(fp.read(), version)

setup(
    name="seqreport",
    version=version["__version__"],
    author="Peter Vegh",
    description="Simple reporting on a set of sequences for documentation purposes",
    long_description=open("pypi-readme.rst").read(),
    long_description_content_type="text/x-rst",
    license="MIT",
    keywords="biology dna",
    packages=find_packages(exclude="docs"),
    include_package_data=True,
    install_requires=[],
)
