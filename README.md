# Thoth v5.3.4

Welcome to the **T**ext enh**H**ancement & **O**ptimization for scien**T**ific researc**H**; or just **Thoth**, as the Egyptian god of writing, wisdom and magic.  

This Python3 package allows you to create, edit and analyze all kinds of text files, with a special focus on ab-initio calculations. In particular, it contains interfaces for [Quantum ESPRESSO](https://www.quantum-espresso.org/) and [Phonopy](https://phonopy.github.io/phonopy/).

Just as the Egyptian god, Thoth is *married* with [Maat](https://github.com/pablogila/Maat), another useful python package to analyze spectral data from your experiments. Check it out!  

Thoth was formally known as InputMaker.


# Installation

Thoth is installed as a regular Python package.
As always, it is recommended to install it in a virtual environment:  
```bash
python3 -m venv .venv
source .venv/bin/activate
```

Then, clone the repository from [GitHub](https://github.com/pablogila/Thoth/) or download it as a ZIP and run inside the `/Thoth/` directory:  
```bash
pip install .
```


# Documentation

Documentation is available locally on the `docs/thoth.html` folder.
An [online documentation](https://pablogila.github.io/InputMaker/) is also available.

## Submodules

This package contains the following submodules for general text edition:
- [file](https://pablogila.github.io/Thoth/thoth/file.html). Manipulate files.
- [find](https://pablogila.github.io/Thoth/thoth/find.html). Search for specific content in a text file.
- [text](https://pablogila.github.io/Thoth/thoth/text.html). Manipulate text files.
- [extract](https://pablogila.github.io/Thoth/thoth/extract.html). Extract data from raw text strings.
- [call](https://pablogila.github.io/Thoth/thoth/call.html). Run bash scripts and related.

Along with the [core](https://pablogila.github.io/Thoth/thoth/core.html) submodule with common utilities.

## Interfaces for ab-initio codes

The following interfaces for ab-initio codes are included:
- [qe](https://pablogila.github.io/Thoth/thoth/qe.html). Interface for [Quantum ESPRESSO](https://www.quantum-espresso.org/) calculations.
- [phonopy](https://pablogila.github.io/Thoth/thoth/phonopy.html). Interface for [Phonopy](https://phonopy.github.io/phonopy/) calculations.

## Compiling the documentation

The documentation can be compiled automatically using [pdoc](https://pdoc.dev/) and Thoth itself, by running:
```shell
python3 makedocs.py
```


## License

> TL;DR: Do what you want with this, as long as you share the source code of your modifications, also under GNU AGPLv3.  

Copyright (C) 2024  Pablo Gila-Herranz  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the attached GNU Affero General Public License for more details.  

