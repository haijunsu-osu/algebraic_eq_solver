#!/usr/bin/env python3
"""
Setup script for algebraic_eq_solver package
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="algebraic_eq_solver",
    version="1.0.0",
    author="AI Assistant",
    author_email="",
    description="Numerical solver for trigonometric systems with advanced singularity handling",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/haijunsu-osu/algebraic_eq_solver",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
        "test": [
            "pytest>=6.0",
            "pytest-cov",
        ],
    },
    entry_points={
        "console_scripts": [
            "trig-solver=algebraic_eq_solver.cli:main",
        ],
    },
)
