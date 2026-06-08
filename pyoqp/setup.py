from setuptools import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="OpenQP",
    version="1.0",
    author="Jingbai Li",
    author_email="lijingbai@zspu.edu.cn",
    description="Python Wrapper for Open Quantum Platform",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    install_requires=[
        # numpy>=2.2 exposes an ffi.from_buffer dangling-pointer bug (zeroed
        # nuclear repulsion on SCF reload); cap until from_buffer lifetimes fixed.
        'numpy>=1.20.0,<2.2',
        'scipy>=1.10.0',
        'geometric>=1.0',
        # 'dftd4>=3.5.0',
        'cffi>=1.16.0',
        'mpi4py>=4.0.0',
        'basis_set_exchange',
        'openmm>=8.1.1',
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
