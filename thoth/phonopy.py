'''
# Description
Functions to work with [Phonopy](https://phonopy.github.io/phonopy/) calculations, along with [Quantum ESPRESSO](https://www.quantum-espresso.org/).

# Index
Under construction, currently porting [PhonoMake](https://github.com/pablogila/phonomake) bash scripts...

---
'''

from .autoload import *
from . import file
from . import find
from . import text
from . import call
from . import extract

def make_scf(folder=None, input:str='relax.in', output:str='relax.out'):
    '''UNDER CONSTRUCTION'''
    if folder:
        call.here(folder)
    else:
        folder = call.here()
    relax_in = file.get(folder, input)
    relax_out = file.get(folder, output)
    scf_in = 'scf.in'
    scf_comment = f'! Automatic SCF input made with thoth.phonopy {version}'
    key1 = 'Begin final coordinates'
    key2 = 'End final coordinates'
    key_species = 'ATOMIC_SPECIES'
    key_species_end = r"(ATOMIC_POSITIONS|CELL_PARAMETERS)"
    key_cell = 'CELL_PARAMETERS'
    text_cell = r'CELL_PARAMETERS {alat}'
    file.from_template(relax_in, scf_in, scf_comment)
    species = find.between(key_species, key_species_end, scf_in, False, 1, True)
    species = '\n' + key_species + '\n' + species
    coords = find.between(key1, key2, relax_out, False, -1, False)
    coords_lines = coords.splitlines()
    coordinates = "\n".join(coords_lines[2:])
    text.delete_under(key_cell, scf_in, -1, -1)
    text.insert_at(species, scf_in, -1)
    text.insert_at(coordinates, scf_in, -1)
    alat_line = find.lines(key_cell, scf_in, -1, 0, False, False)[0]
    alat = extract.number(alat_line, 'alat')
    print(alat)
    text.replace_line(text_cell, key_cell, scf_in, -1, 0, 0, False)
