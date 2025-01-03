from setuptools import find_packages, setup
import re

DESCRIPTION = "The Text enHancement & Optimization for scienTific research with Python, or just ThothPy, allows you to create, modify and analyze all kinds of text files, with a special focus on (but not limited to) ab-initio calculations."

with open('thotpy/core.py', 'r') as f:
    content = f.read()
    version_match = re.search(r"version\s*=\s*'([^']+)'", content)
    if not version_match:
        raise RuntimeError("Unable to find version.")
    VERSION = version_match.group(1)

with open('README.md', 'r') as f:
    LONG_DESCRIPTION = f.read()

setup(
    name='thotpy', 
    version=VERSION,
    author='Pablo Gila-Herranz',
    author_email='pgila001@ikasle.ehu.eus',
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    packages=['thotpy'],
    install_requires=['pandas>=2.2', 'maatpy>=3.1.0'],
    extras_requires={
        'dev': ['pytest', 'twine']
        },
    python_requires='>=3',
    license='AGPL-3.0',
    keywords=['python', 'thot', 'thotpy', 'thoth', 'thothpy', 'text', 'inputmaker', 'DFT', 'Density Functional Theory', 'MD', 'Molecular Dynamics', 'ab-initio', 'Quantum ESPRESSO', 'Phonopy'],
    classifiers= [
        "Development Status :: 5 - Production/Stable",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: Other OS",
    ]
)
