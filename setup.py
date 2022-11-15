from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
	long_description = fh.read()

setup(
	name="LinearHillGillespie",
	version="0.5",
	author="Jacob Hilbert",
	author_email="jacob.hilbert.tree@gmail.com",
	description="Fortran implementation of Gillespie SSA",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/JacobHilbert/LinearHillSSA",
	packages=["LinearHillSSA"],
	install_requires=[
		"numpy >= 1.13.0",
		"scipy >= 1.8.0"
	],
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	license="LICENSE",
	zip_safe=True
)