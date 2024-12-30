'''
# Description
Functions to work with [Quantum ESPRESSO](https://www.quantum-espresso.org/) calculation files.

# Index
- `read_in()`
- `read_out()`
- `read_dir()`
- `read_dirs()`

---
'''


import pandas as pd
import os
from .core import *
from . import file
from . import find
from . import text
from . import extract
from . import call


def read_in(filename) -> dict:
    '''
    Reads an input `filename` from Quantum ESPRESSO,
    returning a dictionary with the input values used.
    The keys are named after the name of the corresponding variable.
    '''
    filepath = file.get(filename)
    data = {}
    lines = find.lines('=', filepath)
    for line in lines:
        line.strip()
        var, value = line.split('=', 1)
        var = var.strip()
        value = value.strip()
        if var.startswith('!'):
            continue
        try:
            value_float = value.replace('d', 'e')
            value_float = value_float.replace('D', 'e')
            value_float = value_float.replace('E', 'e')
            value_float = float(value_float)
            value = value_float
        except ValueError:
            pass # Then it is a string
        data[var] = value
    # Get K_POINTS
    k_points = find.lines(r'(?!\s*!)(k_points|K_POINTS)', filepath, -1, 1, True, True)
    if k_points:
        k_points = k_points[1].strip()
        data['K_POINTS'] = k_points
    # Get ATOMIC_SPECIES
    key_species = r'(?!\s*!)(ATOMIC_SPECIES|atomic_species)'
    atomic_species = None
    if data['ntyp']:
        ntyp = data['ntyp']
        atomic_species_raw = find.lines(key_species, filepath, -1, int(ntyp+1), True)
        # Check that there was no empty line right after the keyword:
        if atomic_species_raw:
            atomic_species_cleaned = []
            for line in atomic_species_raw:
                line = line.strip()
                if not line == '' or not line.startswith('!'):
                    atomic_species_cleaned.append(line)
            atomic_species = atomic_species_cleaned[1:]
            if len(atomic_species) > ntyp:
                atomic_species = atomic_species[:ntyp]
    else:
        key_species_end = r"(?!\s*!)(ATOMIC_POSITIONS|CELL_PARAMETERS)"  # Assuming species go before 
        atomic_species = find.between(key_species, key_species_end, filepath, False, 1, True)
        atomic_species.split()
    if atomic_species:
        data['ATOMIC_SPECIES'] = atomic_species
    # Get CELL_PARAMETERS. Let's take some extra lines just in case there were empty or commented lines in between.
    cell_parameters_raw = find.lines(r'(?!\s*!)(cell_parameters|CELL_PARAMETERS)', filepath, -1, 4, True, True)
    if cell_parameters_raw:
        cell_parameters_cleaned = []
        for line in cell_parameters_raw:
            line = line.strip()
            if not line == '' or not line.startswith('!'):
                cell_parameters_cleaned.append(line)
        if len(cell_parameters_cleaned) > 4:
            cell_parameters_cleaned = cell_parameters_cleaned[:4]
        # extract a possible alat from CELL_PARAMETERS
        alat = extract.number(cell_parameters_cleaned[0])
        if alat:  # This overwrites any possible celldm(1) previously defined!
            data['celldm(1)'] = alat
        cell_parameters = cell_parameters_cleaned[1:]
        data['CELL_PARAMETERS'] = cell_parameters
    # Get ATOMIC_POSITIONS. We assume nat is correct.
    if data['nat']:
        nat = data['nat']
        atomic_positions_raw = find.lines(r'(?!\s*!)(atomic_positions|ATOMIC_POSITIONS)', filepath, -1, int(nat+1), True, True)
        if atomic_positions_raw:
            atomic_positions_cleaned = []
            for line in atomic_positions_raw:
                line.strip()
                if not line == '' or not line.startswith('!'):
                    atomic_positions_cleaned.append(line)
            atomic_positions = atomic_positions_cleaned[1:]
            if len(atomic_positions) > nat:
                atomic_positions = atomic_positions[:nat]
            data['ATOMIC_POSITIONS'] = atomic_positions
    return data


def read_out(filename) -> dict:
    '''
    Reads an output `filename` from Quantum ESPRESSO,
    returning a dict with the following keys:
    `'Energy'` (float), `'Total force'` (float), `'Total SCF correction'` (float),
    `'Runtime'` (str), `'JOB DONE'` (bool), `'BFGS converged'` (bool), `'BFGS failed'` (bool),
    `'Maxiter reached'` (bool), `'Error'` (str), `'Success'` (bool), `'CELL_PARAMETERS_out'` (list of str), `'ATOMIC_POSITIONS_out'` (list of str), `'alat'` (float).
    '''
    filepath = file.get(filename)

    energy_key           = '!    total energy'
    force_key            = 'Total force'
    scf_key              = 'Total SCF correction'
    time_key             = 'PWSCF'
    time_stop_key        = 'CPU'
    job_done_key         = 'JOB DONE.'
    bfgs_converged_key   = 'bfgs converged'
    bfgs_failed_key      = 'bfgs failed'
    maxiter_reached_key  = 'Maximum number of iterations reached'
    error_key            = 'Error in routine'
    cell_parameters_key  = 'CELL_PARAMETERS'
    atomic_positions_key = 'ATOMIC_POSITIONS'

    energy_line          = find.lines(energy_key, filepath, -1)
    force_line           = find.lines(force_key, filepath, -1)
    time_line            = find.lines(time_key, filepath, -1)
    job_done_line        = find.lines(job_done_key, filepath, -1)
    bfgs_converged_line  = find.lines(bfgs_converged_key, filepath, -1)
    bfgs_failed_line     = find.lines(bfgs_failed_key, filepath, -1)
    maxiter_reached_line = find.lines(maxiter_reached_key, filepath, -1)
    error_line           = find.lines(error_key, filepath, -1, 1, True)

    energy: float = None
    force: float = None
    scf: float = None
    time: str = None
    job_done: bool = False
    bfgs_converged: bool = False
    bfgs_failed: bool = False
    maxiter_reached: bool = False
    error: str = ''
    success: bool = False

    if energy_line:
        energy = extract.number(energy_line[0], energy_key)
    if force_line:
        force = extract.number(force_line[0], force_key)
        scf = extract.number(force_line[0], scf_key)
    if time_line:
        time = extract.string(time_line[0], time_key, time_stop_key)
    if job_done_line:
        job_done = True
    if bfgs_converged_line:
        bfgs_converged = True
    if bfgs_failed_line:
        bfgs_failed = True
    if maxiter_reached_line:
        maxiter_reached = True
    if error_line:
        error = error_line[1].strip()
    if job_done and not bfgs_failed and not maxiter_reached and not error:
        success = True

    # CELL_PARAMETERS and ATOMIC_POSITIONS
    cell_parameters = None
    atomic_positions = None
    alat = None
    coordinates_raw = find.between('Begin final coordinates', 'End final coordinates', filepath, False, -1, False)
    if coordinates_raw:
        coordinates_raw = coordinates_raw.splitlines()
        append_cell = False
        append_positions = False
        cell_parameters = []
        atomic_positions = []
        for line in coordinates_raw:
            line = line.strip()
            if cell_parameters_key in line:
                append_cell = True
                append_positions = False
                alat = extract.number(line, 'alat')
            elif atomic_positions_key in line:
                append_cell = False
                append_positions = True
            if line == '' or line.startswith('!'):
                continue
            if append_cell:
                cell_parameters.append(line)
            elif append_positions:
                atomic_positions.append(line)
        cell_parameters = cell_parameters[1:]
        atomic_positions = atomic_positions[1:]

    output = {
        'Energy'                : energy,
        'Total force'           : force,
        'Total SCF correction'  : scf,
        'Runtime'               : time,
        'JOB DONE'              : job_done,
        'BFGS converged'        : bfgs_converged,
        'BFGS failed'           : bfgs_failed,
        'Maxiter reached'       : maxiter_reached,
        'Error'                 : error,
        'Success'               : success,
        'CELL_PARAMETERS_out'   : cell_parameters,
        'ATOMIC_POSITIONS_out'  : atomic_positions,
        'alat'                  : alat,
    }
    return output


def read_dir(folder,
             input_str:str='.in',
             output_str:str='.out'
             ) -> dict:
    '''
    Takes a `folder` containing a Quantum ESPRESSO calculation,
    and returns a dictionary containing the input parameters and output results.
    Input and output files are determined automatically,
    but must be specified with `input_str` and `output_str` if more than one file ends with `.in` or `.out`.
    To extract values only from the input or only from the output, check `read_in()` and `read_out()`.
    '''
    input_file = file.get(folder, input_str)
    output_file = file.get(folder, output_str)
    if not input_file and not output_file:
        return None
    if input_file:
        dict_in = read_in(input_file)
        if not output_file:
            return dict_in
    if output_file:
        dict_out = read_out(output_file)
        if not input_file:
            return dict_out
    # Merge both dictionaries
    merged_dict = {**dict_in, **dict_out}
    return merged_dict


def read_dirs(directory,
              input_str:str='.in',
              output_str:str='.out',
              calc_splitter='_',
              calc_type_index=0,
              calc_id_index=1
              ) -> None:
    '''
    Calls recursively `read_dir()`, reading Quantum ESPRESSO calculations
    from all the subfolders inside the given `directory`.
    The results are saved to CSV files inside the current directory.
    Input and output files are determined automatically, but must be specified with
    `input_str` and `output_str` if more than one file ends with `.in` or `.out`.

    To properly group the calculations per type, saving separated CSVs for each calculation type,
    you can modify `calc_splitter` ('_' by default), `calc_type_index` (0) and `calc_id_index` (1).
    With these default values, a subfolder named './CalculationType_CalculationID_AdditionalText/'
    will be interpreted as follows:
    - Calculation type: 'CalculationType' (The output CSV will be named after this)
    - CalculationID: 'CalculationID' (Stored in the 'ID' column of the resulting dataframe)

    If everything fails, the subfolder name will be used.
    '''
    print(f'Reading all Quantum ESPRESSO calculations from {directory} ...')
    folders = file.get_list(directory)
    if not folders:
        raise FileNotFoundError('The directory is empty!')
    # Separate calculations by their title in an array
    calc_types = []
    folders.sort()
    for folder in folders:
        if not os.path.isdir(folder):
            folders.remove(folder)
            continue
        folder_name = os.path.basename(folder)
        try:
            calc_name = folder_name.split(calc_splitter)[calc_type_index]
        except:
            calc_name = folder_name
        if not calc_name in calc_types:
            calc_types.append(calc_name)
    len_folders = len(folders)
    total_success_counter = 0
    for calc in calc_types:
        len_calcs = 0
        success_counter = 0
        results = pd.DataFrame()
        for folder in folders:
            if not calc in folder:
                continue
            len_calcs += 1
            folder_name = os.path.basename(folder)
            try:
                calc_id = folder_name.split(calc_splitter)[calc_id_index]
            except:
                calc_id = folder_name
            df = pd.DataFrame.from_dict(read_dir(folder, input_str, output_str))
            if df is None:
                continue
            # Join input and output in the same dataframe
            df.insert(0, 'ID', calc_id)
            df = df.dropna(axis=1, how='all')
            results = pd.concat([results, df], axis=0, ignore_index=True)
            if df['Success'][0]:
                success_counter += 1
                total_success_counter += 1
        results.to_csv(os.path.join(directory, calc+'.csv'))
        print(f'Saved to CSV: {calc} ({success_counter} successful calculations out of {len_calcs})')
    print(f'Total successful calculations: {total_success_counter} out of {len_folders}')


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

