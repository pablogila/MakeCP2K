'''
# Description
Common utilities, autoloaded directly as `thotpy.value`.

# Index
- `version`
- `help()`

---
'''


version = 'v5.4.2'
'''
Package version, using semantic versioning to indicate breaking changes,
as in v<MAJOR>.<MINOR>.<PATCH>.
'''


def help() -> None:
    '''Print the help message.'''
    help_message ='''----------------------------------------------
# ThotPy - Text enHancement & Optimization for scienTific researcH

## Submodules
ThotPy contains the following submodules for general text edition:
- thotpy.file. Manipulate files.
- thotpy.find. Search for specific content in a text file.
- thotpy.text. Manipulate text files.
- thotpy.extract. Extract data from raw text strings.
- thotpy.call. Run bash scripts and related.
Along with the thotpy.core submodule with common utilities.

## Interfaces for ab-initio codes
The following interfaces for ab-initio codes are included:
- thotpy.qe. Interface for Quantum ESPRESSO calculations.
- thotpy.phonopy. Interface for Phonopy calculations.

## Full documentation
Check the full ThotPy documentation on https://pablogila.github.io/ThotPy/
----------------------------------------------'''
    print(help_message)
    return None

