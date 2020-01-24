import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="utils", 
    version="1.0",
    author="Wendy L. Williams",
    author_email="wllwen007@gmail.com",
    description="A set of python utilities",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wllwen007/utils",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    python_requires='>=3.6',
)
