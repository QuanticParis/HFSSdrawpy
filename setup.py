import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="HFSSdrawpy",  # Replace with your own username
    version="0.9",
    author="HQCteam",
    author_email="teamhqc@gmail.com",
    description="Drawing circuits in HFSS using python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/leghtas/HFSSdrawpy",
    packages=setuptools.find_packages(exclude=["tests", "tests.*"]),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "Pint>=0.10",
        "numpy",
        "sympy>=1.5.1",
        "Sphinx"
        'gdspy>=1.5.2 ; platform_system!="Windows"',
        'pywin32>=227 ; platform_system=="Windows"',
    ],
)
