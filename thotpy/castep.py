'''
# Description
Functions to work with [CASTEP](https://castep-docs.github.io/castep-docs/) calculation files.

# Index
- `read_castep()`

---
'''


from . import file
from . import find
from . import extract


def read_castep(filename) -> dict:
    '''
    Reads a CASTEP output file, specified in `filename`.
    Returns a dictionary with the following keys:
    `'Enthalpy'` (LBFGS: Final Enthalpy, in kJ/mol),
    `'Energy'` (Total energy corrected for finite basis set, in eV),
    `'Space group'`, `'Volume'` (Angstrom^3), `'Density'` (amu/Angstrom^3), `'Density_g'` (g/cm^3),
    `'A'`, `'B'`, `'C'` (Angstroms), `'Alpha'`, `'Beta'`, `'Gamma'` (Degrees).
    '''
    file_castep = file.get(filename)
    # Initial definitions
    enthalpy    = None
    energy      = None
    space_group = None
    volume      = None
    density     = None
    density_g   = None
    a           = None
    b           = None
    c           = None
    alpha       = None
    beta        = None
    gamma       = None
    # Find the output values in the file
    enthalpy_str    = find.lines('LBFGS: Final Enthalpy     =', file_castep, -1)
    energy_str      = find.lines('Total energy corrected for finite basis set =', file_castep, -1)
    space_group_str = find.lines('Space group of crystal =', file_castep, -1)
    volume_str      = find.lines('Current cell volume =', file_castep, -1)
    density_str     = find.lines('density =', file_castep, -1, 1)
    a_str           = find.lines('a =', file_castep, -1)
    b_str           = find.lines('b =', file_castep, -1)
    c_str           = find.lines('c =', file_castep, -1)
    
    if enthalpy_str:
        enthalpy = extract.number(enthalpy_str[0], 'LBFGS: Final Enthalpy')
    if energy_str:
        energy = extract.number(energy_str[0], 'Total energy corrected for finite basis set')
    if space_group_str:
        # Avoid little stupid errors
        space_group_str = space_group_str.replace(',','.')
        space_group = extract.string(space_group_str[0], 'Space group of crystal')
    if volume_str:
        volume = extract.number(volume_str[0], 'Current cell volume')
    if density_str:
        density = extract.number(density_str[0], 'density')
        density_g = extract.number(density_str[1], '')
    if a:
        a = extract.number(a_str, 'a')
        alpha = extract.number(a_str, 'alpha')
    if b:
        b = extract.number(b_str, 'b')
        beta = extract.number(b_str, 'beta')
    if c:
        c = extract.number(c_str, 'c')
        gamma = extract.number(c_str, 'gamma')
    # Return as a dict
    dictionary = {
        'Enthalpy'    : enthalpy,
        'Energy'      : energy,
        'Space group' : space_group,
        'Volume'      : volume,
        'Density'     : density,
        'Density_g'   : density_g,
        'A'           : a,
        'B'           : b,
        'C'           : c,
        'Alpha'       : alpha,
        'Beta'        : beta,
        'Gamma'       : gamma,
    }
    return dictionary

