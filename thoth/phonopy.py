'''
# Description
Functions to work with [Phonopy](https://phonopy.github.io/phonopy/) calculations,
along with [Quantum ESPRESSO](https://www.quantum-espresso.org/).

# Index
- `make()` Makes the inputs, by running the previous functions.
- `sbatch()` Sbatch'es the supercell calculations.
- `scf_from_relax()`
- `supercells_from_scf()`
- `scf_header_to_supercells()`
- `check_slurm_template()`

---
'''


import os
import re
from .autoload import *
from . import file
from . import find
from . import text
from . import call
from . import extract


def make(dimension:str='2 2 2',
         folder:str=None,
         relax_in:str='relax.in',
         relax_out:str='relax.out',
         slurm_template:str='scf.slurm'
         ) -> None:
    '''
    Starting on a given `folder` (CWD if none) from the `relax_in` and `relax_out` (default ones),
    creates the supercells of a `dimension` (`2 2 2` by default)
    needed for the Phonopy calculations with Quantum ESPRESSO.
    It runs sequentially `scf_from_relax()`, `supercells_from_scf()` and `scf_header_to_supercells()`.
    Finally, it checks the `slurm_template` with `check_slurm_template()`.
    '''
    print(f'\nthoth.phonopy {version}\n'
          'Creating all supercell inputs with Phonopy for Quantum ESPRESSO...\n')
    scf_from_relax(folder, relax_in, relax_out)
    supercells_from_scf(dimension, folder)
    scf_header_to_supercells(folder)
    print('\n------------------------------------------------------\n'
          'PLEASE CHECH BELOW THE CONTENT OF THE supercell-001.in\n'
          '------------------------------------------------------\n')
    call.bash('head -n 100 supercell-001.in')
    print('\n------------------------------------------------------\n'
          'PLEASE CHECH THE CONTENT OF THE supercell-001.in\n'
          'The first 100 lines of the input were printed above!\n'
          '------------------------------------------------------\n\n'
          'If it seems correct, run the calculations with thoth.phonopy.sbatch()\n')
    check_slurm_template(folder, slurm_template)
    return None


def sbatch(folder=None, slurm_template:str='scf.slurm') -> None:
    '''
    Launch all your supercell calculations to a cluster using a SLURM manager.
    Launched from a `folder` (CWD if empty), using a `slurm_template` (`scf.slurm` by default)
    '''
    print(f'\nthoth.phonopy {version}\n'
          'Sbatching all supercells...\n')
    key_input = 'INPUT_FILE'
    key_output = 'OUTPUT_FILE'
    key_jobname = 'JOB_NAME'
    id_pattern = re.compile(r'supercell-(\d\d\d).in')
    slurm_folder = 'slurms'
    folder = call.here(folder)
    # Get supercells and abort if not found
    supercells = file.get_list(folder, 'supercell-')
    if len(supercells) == 0:
        raise FileNotFoundError('Supercells were not found! Did you run thoth.phonopy.make() ?')
    call.bash(f"mkdir {slurm_folder}", folder, True, True)
    # Get the template
    slurm_file = check_slurm_template(folder, slurm_template)
    if not slurm_file:
        print(f'Aborting... Please correct {slurm_template}\n')
        return None
    for supercell in supercells:
        # Get the file ID
        match = re.search(id_pattern, supercell)
        calc_id = match.group(1)
        # Create slurm file for this supercell
        slurm_id = 'scf-' + calc_id + '.slurm'
        supercell = os.path.basename(supercell)
        supercell_out = supercell.replace('.in', '.out')
        fixing_dict = {
            key_jobname: calc_id,
            key_input: supercell,
            key_output: supercell_out
        }
        file.from_template(slurm_file, slurm_id, None, fixing_dict)
        #call.bash(f"echo {slurm_id}")  # This line may be useful for testing!
        call.bash(f"sbatch {slurm_id}", folder, True, False)
        call.bash(f"mv {slurm_id} {slurm_folder}", folder, False, True)  # Do not raise error if we can't move the file
    print(f'\nDone! Temporary slurm files were moved to /{slurm_folder}/\n')


def scf_from_relax(folder:str=None,
                   relax_in:str='relax.in',
                   relax_out:str='relax.out'
                   ) -> None:
    '''
    Create a Quantum ESPRESSO `scf.in` file from a previous relax calculation.
    If no `folder` is provided, the current working directory is used.
    The `relax_in` and `relax_out` files by default are `relax.in` and `relax.out`,
    update the names if necessary.
    '''
    folder = call.here(folder)
    relax_in = file.get(folder, relax_in)
    relax_out = file.get(folder, relax_out)
    # Terminal feedback
    print(f'\nthoth.phonopy {version}\n'
          f'Creating Quantum ESPRESSO SCF input from previous relax calculation:\n'
          f'{relax_in}\n{relax_out}\n')
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
    try:
        text.delete_under(key_species, scf_in, -1, -1)
    except:  # It might be deleted already!
        pass
    try:
        text.delete_under(key_cell, scf_in, -1, -1)
    except:  # It might be deleted already!
        pass
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
    text.replace_line(r'CELL_PARAMETERS alat', key_cell, scf_in, -1)
    text.replace_line('', r'celldm\(\d\)\s*=', scf_in, 0, 0, 0, True)
    text.replace_line('', r'[ABC]\s*=', scf_in, 0, 0, 0, True)
    text.insert_under(alat_text, r'(&SYSTEM|&system)', scf_in, 1, 0, True)
    text.replace_line('  ibrav = 0', r'ibrav\s*=', scf_in, 0, 0, 0, True)
    text.replace_line("  calculation = 'scf'", r'calculation\s*=', scf_in, 0, 0, 0, True)
    # Terminal feedback
    print(f'Created input SCF file at:'
          f'{scf_in}\n')
    return None


def supercells_from_scf(dimension:str='2 2 2',
                        folder:str=None,
                        scf:str='scf.in'
                        ) -> None:
    '''
    Creates supercells of a given `dimension` (`2 2 2` by default) inside a `folder`,
    from a Quantum ESPRESSO `scf` input (`scf.in` by default).
    '''
    print(f'\nthoth.phonopy {version}\n')
    folder = call.here(folder)
    scf_in = file.get(folder, scf, True)
    if scf_in is None:
        raise ValueError('No SCF input found in path!')
    call.bash(f'phonopy --qe -d --dim="{dimension}" -c {scf_in}')
    return None


def scf_header_to_supercells(folder:str=None,
                             scf:str='scf.in',
                             ) -> None:
    '''
    '''
    print(f'\nthoth.phonopy {version}\n'
          f'Adding headers to Phonopy supercells for Quantum ESPRESSO...\n')
    folder = call.here(folder)
    # Check if the header file, the scf.in, exists
    scf_file = file.get(folder, scf, True)
    if scf_file is None:
        raise ValueError('No header file found in path!')
    # Check if the supercells exist
    supercells = file.get_list(folder, 'supercell-')
    if supercells is None:
        raise ValueError('No supercells found in path!')
    # Check if the supercells contains '&CONTROL' and abort if so
    supercell_sample = supercells[0]
    is_control = find.lines(r'(&CONTROL|&control)', supercell_sample, 1, 0, False, True)
    if is_control:
        raise ValueError('Supercells already contain &CONTROL!')
    # Check if the keyword is in the scf file
    is_header = find.lines(r'ATOMIC_SPECIES', scf_file, 1, 0, False, False)
    if not is_header:
        raise ValueError('No ATOMIC_SPECIES found in header!')
    # Copy the scf to a temp file
    temp_scf = '_scf_temp.in'
    file.copy(scf_file, temp_scf)
    # Remove the top content from the temp file
    text.delete_under('ATOMIC_SPECIES', '_scf_temp.in', -1, -1, False)
    # Find the new number of atoms and replace the line
    updated_values = find.lines('ibrav', supercell_sample, 1)  # !    ibrav = 0, nat = 384, ntyp = 5
    if not updated_values:
        print("!!! Okay listen, this is weird. This code should never be running, "
              "but for some reson we couldn't find the updated values in the supercells. "
              "Please, introduce the NEW NUMBER OF ATOMS in the supercells manually (int):")
        nat = int(input('nat = '))
    else:
        nat = extract.number(updated_values[0], 'nat')
    text.replace_line(f'  nat = {nat}', r'nat\s*=', temp_scf, 1, 0, 0, True)
    # Remove the lattice parameters, since Phonopy already indicates units
    text.replace_line('', r'celldm\(\d\)\s*=', temp_scf, 0, 0, 0, True)
    text.replace_line('', r'[ABC]\s*=', temp_scf, 0, 0, 0, True)
    # Add the header to the supercells
    with open(temp_scf, 'r') as f:
        header = f.read()
    for supercell in supercells:
        text.insert_at(header, supercell, 0)
    # Remove the temp file
    os.remove('_scf_temp.in')
    print('Done!')
    return None


def check_slurm_template(folder=None, slurm_template:str='scf.slurm') -> str:
    '''
    Check a `slurm_template` inside `folder`.
    The current working directory is used if `folder` is not provided.
    If the file does not exist or is invalid, creates a `scf_EXAMPLE.slurm` file for reference.
    '''
    folder = call.here(folder)
    slurm_example = 'scf_EXAMPLE.slurm'
    new_slurm_file = os.path.join(folder, slurm_example)
    # Default slurm template
    content =f'''
# Automatic slurm created with thoth.phonopy {version}

#!/bin/bash
#SBATCH --partition=general
#SBATCH --qos=regular
#SBATCH --job-name=JOB_NAME
#SBATCH --ntasks=32
#SBATCH --time=1-00:00:00
#SBATCH --mem=128G
# #SBATCH --mail-user=YOUR@EMAIL
# #SBATCH --mail-type=END

module purge
module load QuantumESPRESSO/7.3-foss-2023a

mpirun pw.x -inp INPUT_FILE > OUTPUT_FILE
'''
    # If the slurm template does not exist, create one
    slurm_file = file.get(folder, slurm_template, True)
    if not slurm_file:
        with open(new_slurm_file, 'w') as f:
            f.write(content)
        print(f'!!! WARNING:  Slurm template missing, so an example was generated automatically:\n'
              f'{new_slurm_file}\n'
              f'PLEASE CHECK it, UPDATE it and RENAME it to {slurm_template}\n'
              'before running thoth.phonopy.sbatch()\n')
        return None
    # Check that the slurm file contains the INPUT_FILE, OUTPUT_FILE and JOB_NAME keywords
    key_input = find.lines('INPUT_FILE', slurm_file)
    key_output = find.lines('OUTPUT_FILE', slurm_file)
    key_jobname = find.lines('JOB_NAME', slurm_file)
    missing = []
    if not key_input:
        missing.append('INPUT_FILE')
    if not key_output:
        missing.append('OUTPUT_FILE')
    if not key_jobname:
        missing.append('JOB_NAME')
    if len(missing) > 0:
        with open(new_slurm_file, 'w') as f:
            f.write(content)
        print('!!! WARNING:  Some keywords were missing from your slurm template,\n'
              f'PLEASE CHECK the example at:\n'
              f'{new_slurm_file}\n'
              'before running thoth.phonopy.sbatch()\n'
              f'The following keywords were missing from your {slurm_template}:')
        for key in missing:
            print(key)
        print('')
        return None
    return slurm_file  # Ready to use!

