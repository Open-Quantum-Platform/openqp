from setuptools import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="PyOpenQP",
    version="1.0",
    author="Jingbai Li",
    author_email="lijingbai@zspu.edu.cn",
    description="Python Wrapper for Open Quantum Platform",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    install_requires=[
        'numpy>=1.20.0',
        'scipy>=1.10.0',
        'libdlfind>=0.0.3',
        # 'dftd4>=3.5.0',
        'cffi>=1.16.0',
#        'mpi4py>=4.0.0',
    ],
    extras_require={
    },
    packages=find_packages(),
    include_package_data=False,
    package_data={"oqp": ["*.py"]},
    entry_points={
        'console_scripts': [
            'openqp=oqp.pyoqp:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
    keywords=["science", "MRSF-TDDFT"],
)
