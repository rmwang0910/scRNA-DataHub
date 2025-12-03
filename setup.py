#!/usr/bin/env python3
"""Setup script for scRNA-DataHub"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_file = Path(__file__).parent / 'README.md'
long_description = readme_file.read_text(encoding='utf-8') if readme_file.exists() else ''

# Read requirements
requirements_file = Path(__file__).parent / 'requirements.txt'
requirements = []
if requirements_file.exists():
    requirements = [
        line.strip()
        for line in requirements_file.read_text().splitlines()
        if line.strip() and not line.startswith('#')
    ]

setup(
    name='scrna-datahub',
    version='1.0.0',
    description='Universal single-cell RNA-seq data reader and converter',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Wang Ruiming',
    author_email='your.email@example.com',
    url='https://github.com/yourusername/scRNA-DataHub',
    license='MIT',
    
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    
    install_requires=requirements,
    
    extras_require={
        'full': [
            'loompy>=3.0.6',
            'zarr>=2.10.0',
            'openpyxl>=3.0.0',
            'harmonypy>=0.0.9',
            'bbknn>=1.5.0',
            'scanorama>=1.7.0',
        ],
        'dev': [
            'pytest>=6.0.0',
            'pytest-cov>=2.12.0',
            'black>=21.0',
            'flake8>=3.9.0',
            'mypy>=0.900',
        ]
    },
    
    entry_points={
        'console_scripts': [
            'scrna-datahub=universal_reader:main',
        ],
    },
    
    python_requires='>=3.8',
    
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    
    keywords='single-cell RNA-seq data-processing bioinformatics',
    
    project_urls={
        'Bug Reports': 'https://github.com/yourusername/scRNA-DataHub/issues',
        'Source': 'https://github.com/yourusername/scRNA-DataHub',
        'Documentation': 'https://github.com/yourusername/scRNA-DataHub/docs',
    },
)

