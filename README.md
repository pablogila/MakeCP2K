# ThotPy v5.5.0

Welcome to the **T**ext en**H**ancement & **O**ptimization for scien**T**ific research with **Py**thon; or just **ThotPy**, as the modern incarnation of the ancient Egyptian god of writing and wisdom, [Thot](https://en.wikipedia.org/wiki/Thoth).  

This Python3 package allows you to easily create, edit and analyse all kinds of text files, with a special focus on *ab-initio* calculations. In particular, it contains interfaces for [Quantum ESPRESSO](https://www.quantum-espresso.org/), [Phonopy](https://phonopy.github.io/phonopy/) and [CASTEP](https://castep-docs.github.io/castep-docs/).  

Check the [full documentation online](https://pablogila.github.io/ThotPy/).  


# Installation

As always, it is recommended to install your packages in a virtual environment:  
```bash
python3 -m venv .venv
source .venv/bin/activate
```

## With pip

The fastest way to install ThotPy is to use pip:  
```bash
pip install thotpy
```

## From source

Optionally, you can install ThotPy from the [GitHub repository](https://github.com/pablogila/ThotPy/).  

Like its ancient Egyptian counterpart, ThotPy is *married* to [MaatPy](https://github.com/pablogila/MaatPy), another Python package with useful physico-chemical definitions and spectral analysis tools. Pip automatically installs MaatPy as a dependency, but it must be installed before you can build from source.  
To install the dependencies:  
```bash
pip install pandas maatpy
```

Then clone the repository or download the [latest stable release](https://github.com/pablogila/ThotPy/tags) as a ZIP, unzip it, and run inside the `ThotPy/` directory:  
```bash
pip install .
```


# Documentation

Check the [full ThotPy documentation online](https://pablogila.github.io/ThotPy/).  
An offline version of the documentation is available in `docs/thotpy.html`.  
Code examples are included in the `examples/` folder.  

## Submodules

ThotPy contains the following submodules for general text operations:  
- [file](https://pablogila.github.io/ThotPy/thotpy/file.html). Manipulate files.
- [find](https://pablogila.github.io/ThotPy/thotpy/find.html). Search for specific content in a text file.
- [text](https://pablogila.github.io/ThotPy/thotpy/text.html). Manipulate text files.
- [extract](https://pablogila.github.io/ThotPy/thotpy/extract.html). Extract data from raw text strings.
- [call](https://pablogila.github.io/ThotPy/thotpy/call.html). Run bash scripts and related.

Along with the [core](https://pablogila.github.io/ThotPy/thotpy/core.html) submodule with common utilities.

## Interfaces for ab-initio codes

The following interfaces for *ab-initio* codes are included:
- [qe](https://pablogila.github.io/ThotPy/thotpy/qe.html). Interface for [Quantum ESPRESSO](https://www.quantum-espresso.org/) calculations.
- [phonopy](https://pablogila.github.io/ThotPy/thotpy/phonopy.html). Interface for [Phonopy](https://phonopy.github.io/phonopy/) calculations.
- [castep](https://pablogila.github.io/ThotPy/thotpy/castep.html) Interface for [CASTEP](https://castep-docs.github.io/castep-docs/) calculations.


# Contributing

If you are interested in opening an issue or a pull request, please feel free to do so on [GitHub](https://github.com/pablogila/ThotPy/).  

For major changes, please get in touch first to discuss the details.  

## Code style

Please try to follow some general guidelines:  
- Use a code style consistent with the rest of the project.  
- Include docstrings to document new additions.  
- Include tests for new features or modifications.  
- Arrange function arguments by order of relevance. Most implemented functions follow something similar to `function(file, key/s, value/s, optional)`.  

## Testing with PyTest

If you are modifying the source code, you should run the automated tests of the `tests/` folder to check that everything works as intended.
To do so, first install PyTest in your environment,
```bash
pip install pytest
```

And then run PyTest inside the `ThotPy/` directory,
```bash
pytest -vv
```

## Compiling the documentation

The documentation can be compiled automatically to `docs/thotpy.html` with [pdoc](https://pdoc.dev/) and ThotPy itself, by running:
```shell
python3 makedocs.py
```


# License

Copyright (C) 2024  Pablo Gila-Herranz  
This program is free software: you can redistribute it and/or modify
it under the terms of the **GNU Affero General Public License** as published
by the Free Software Foundation, either version **3** of the License, or
(at your option) any later version.  
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the attached GNU Affero General Public License for more details.  

