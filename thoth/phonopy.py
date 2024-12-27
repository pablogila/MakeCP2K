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
    '''
    Create a Quantum ESPRESSO `scf.in` file from a previous relax calculation.
    If no `folder` is provided, the current working directory is used.
    The `input` and `output` files by default are `relax.in` and `relax.out`, respectively;
    update the names if necessary.
    '''
    # Gather files and create the scf.in from the previous relax.in
    if folder:
        call.here(folder)
    else:
        folder = call.here()
    relax_in = file.get(folder, input)
    relax_out = file.get(folder, output)
    scf_in = 'scf.in'
    comment = f'! Automatic SCF input made with thoth.phonopy {version}'
    file.from_template(relax_in, scf_in, comment)
    # Key definitions
    key_system = r'(&SYSTEM|&system)'
    key_species = 'ATOMIC_SPECIES'
    key_species_end = r"(ATOMIC_POSITIONS|CELL_PARAMETERS)"
    key_cell = 'CELL_PARAMETERS'
    key_positions = 'ATOMIC_POSITIONS'
    key_alat = r'celldm(\d)\s*='
    key_A = r'[ABC]\s*='
    key_calculation = r'calculation\s*='
    calculation = "  calculation = 'scf'"
    # Save the ATOMIC_SPECIES from the relax input
    atomic_species = find.between(key_species, key_species_end, scf_in, False, 1, True)
    atomic_species = '\n' + key_species + '\n' + atomic_species
    # Delete the old CELL_PARAMETERS, ATOMIC_POSITIONS and ATOMIC_SPECIES
    text.delete_under(key_positions, scf_in, -1, -1)
    text.delete_under(key_species, scf_in, -1, -1)
    text.delete_under(key_cell, scf_in, -1, -1)
    # Get the final coordinates from the relax output
    coordinates = find.between('Begin final coordinates', 'End final coordinates', relax_out, False, -1, False)
    coordinates = coordinates.splitlines()
    coordinates = "\n".join(coordinates[2:])
    # Insert the new ATOMIC_SPECIES and coordinates. Atomic species will go first
    text.insert_at(atomic_species, scf_in, -1)
    text.insert_at(coordinates, scf_in, -1)
    # Get the lattice parameterts
    alat_line = find.lines(key_cell, scf_in, -1, 0, False, False)[0]
    alat = extract.number(alat_line, 'alat')
    alat_text = f'  celldm(1) = {alat}'
    # Update lattice parameters, ibrav and calculation
    text.replace_line(r'CELL_PARAMETERS {alat}', key_cell, scf_in, -1)
    text.replace_line('', key_alat, scf_in, 0, 0, 0, True)
    text.replace_line('', key_A, scf_in, 0, 0, 0, True)
    text.insert_under(alat_text, key_system, scf_in, 1, 0, True)
    text.replace_line('  ibrav = 0', r'ibrav\s*=', scf_in, 0, 0, 0, True)
    text.replace_line(calculation, key_calculation, scf_in, 0, 0, 0, True)

