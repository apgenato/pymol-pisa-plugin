from setuptools import setup, find_packages

setup(
    name='pymol-pisa-plugin',
    version='0.1.0',
    author='Your Name',
    author_email='your.email@example.com',
    description='A PyMOL plugin for analyzing interactions between receptor and ligand molecules.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/pymol-pisa-plugin',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'pymol',  # Specify the version if needed
        'numpy',  # Add any other dependencies required
        'matplotlib'  # Example of a potential dependency
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)