'''
# Description
Functions to work with [Phonopy](https://phonopy.github.io/phonopy/) calculations, along with [Quantum ESPRESSO](https://www.quantum-espresso.org/).

# Index
Under construction, currently porting [PhonoMake](https://github.com/pablogila/phonomake) bash scripts...

- `make_scf()`  Working!
- `make_supercells()`  Working!
- `make_headers()` TO-DO
...

---
'''


import os
from .autoload import *
from . import file
from . import find
from . import text
from . import call
from . import extract


def make_scf(folder:str=None,
             input:str='relax.in',
             output:str='relax.out'
             ) -> str:
    '''
    Create a Quantum ESPRESSO `scf.in` file from a previous relax calculation.
    If no `folder` is provided, the current working directory is used.
    The `input` and `output` files by default are `relax.in` and `relax.out`, respectively;
    update the names if necessary.
    Returns the path of the created `scf.in` file.
    '''
    if folder:
        call.here(folder)
    else:
        folder = call.here()
    relax_in = file.get(folder, input)
    relax_out = file.get(folder, output)
    # Terminal feedback
    print(f'thoth.phonopy {version}\n'
          f'Creating Quantum ESPRESSO SCF input from previous relax calculation:\n'
          f'{relax_in}\n{relax_out}')
    # Create the scf.in from the previous relax.in
    scf_in = 'scf.in'
    comment = f'! Automatic SCF input made with thoth.phonopy {version}'
    file.from_template(relax_in, scf_in, comment)
    scf_in = file.get(folder, scf_in)
    # Save the ATOMIC_SPECIES from the relax input
    key_species = 'ATOMIC_SPECIES'
    key_species_end = r"(ATOMIC_POSITIONS|CELL_PARAMETERS)"
    atomic_species = find.between(key_species, key_species_end, scf_in, False, 1, True)
    atomic_species = '\n' + key_species + '\n' + atomic_species
    # Delete the old CELL_PARAMETERS, ATOMIC_POSITIONS and ATOMIC_SPECIES
    key_cell = 'CELL_PARAMETERS'
    text.delete_under('ATOMIC_POSITIONS', scf_in, -1, -1)
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
    text.replace_line('', r'celldm(\d)\s*=', scf_in, 0, 0, 0, True)
    text.replace_line('', r'[ABC]\s*=', scf_in, 0, 0, 0, True)
    text.insert_under(alat_text, r'(&SYSTEM|&system)', scf_in, 1, 0, True)
    text.replace_line('  ibrav = 0', r'ibrav\s*=', scf_in, 0, 0, 0, True)
    text.replace_line("  calculation = 'scf'", r'calculation\s*=', scf_in, 0, 0, 0, True)
    # Terminal feedback
    print(f'Created input SCF file at:\n'
          f'{scf_in}\n')
    return scf_in


def make_supercells(scf:str=None, dimension:str='2 2 2'):
    '''
    Creates supercells with Phonopy from a Quantum ESPRESSO `scf` input (`scf.in` by default).
    The `dimension` of the supercells has the format `X Y Z`, and is `2 2 2` by default.
    '''
    if scf:
        folder = os.path.dirname(scf)
        call.here(folder)
    else:
        folder = call.here()
    call.bash(f'phonopy --qe -d --dim="{dimension}" -c scf.in')


def make_headers(scf:str=None) -> None:
    '''
    '''
    if scf:
        folder = os.path.dirname(scf)
        call.here(folder)
    else:
        folder = call.here()
    # Check if the supercells exist
    pass   ########################  UNDER CONSTRUCTION


def run():
    '''
    THIS IS A DEMO
    '''
    print(f'\nWelcome to Thoth.phonopy {version}\n'
          'Enter below the desired workflow to run,\n'
          'press enter to run everything automatically!\n\n'
          '1. Make input SCF file\n'
          '2. Make supercells\n'
          '3. Add headers to supercells\n'
          '4. Run Phonopy calculation with SLURM\n'
          '5. Create force constants\n'
          '6. SCANCEL current jobs\n'
          '7. Fix unfinished jobs\n'
          '\nPress any other key to exit...\n\n')
    workflow = input('> ')
    if workflow == '1':
        make_scf()
    elif workflow == '2':
        make_supercells()
    elif workflow == '3':
        make_headers()
    elif workflow == '4':
        run_phonopy()
    elif workflow == '5':
        create_force_constants()
    elif workflow == '6':
        scancel_jobs()
    elif workflow == '7':
        fix_unfinished_jobs()
    elif workflow == '':
        make_scf()
        make_supercells()
        make_headers()
        run_phonopy()
        create_force_constants()
        scancel_jobs()
        fix_unfinished_jobs()
    else:
        print('Exiting...')
        return None
