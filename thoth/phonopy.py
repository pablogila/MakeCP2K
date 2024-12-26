'''
# Description
Functions to work with [Phonopy](https://phonopy.github.io/phonopy/) calculations, along with [Quantum ESPRESSO](https://www.quantum-espresso.org/).

# Index
Under construction, currently porting [PhonoMake](https://github.com/pablogila/phonomake) bash scripts...

---
'''

from . import file
from . import text
from . import call
from .common import *

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
    key_species_end = r'?\s*\w{4}'
    key_cell = 'CELL_PARAMETERS'
    file.from_template(relax_in, scf_in, scf_comment)
    species = text.find_between(key_species, key_species_end, scf_in, -1, True)
    species.insert(0, key_species)
    coordinates = text.find_between(key1, key2, relax_out, -1, True)
    coordinates = coordinates[2:]
    text.delete_under(key_cell, scf_in, -1, -1)
    text.insert_at(species, scf_in, -1)
    