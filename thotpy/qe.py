'''
# Description
Functions to work with [Quantum ESPRESSO](https://www.quantum-espresso.org/) calculation files.

# Index
- `pw_description`
- `read_in()`
- `read_out()`
- `read_dir()`
- `read_dirs()`
- `set_value()`
- `scf_from_relax()`

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


pw_description = {
    '&CONTROL' : ['calculation', 'title', 'verbosity', 'restart_mode', 'wf_collect', 'nstep', 'iprint', 'tstress', 'tprnfor', 'dt', 'outdir', 'wfcdir', 'prefix', 'lkpoint_dir', 'max_seconds', 'etot_conv_thr', 'forc_conv_thr', 'disk_io', 'pseudo_dir', 'tefield', 'dipfield', 'lelfield', 'nberrycyc', 'lorbm', 'lberry', 'gdir', 'nppstr', 'gate', 'twochem', 'lfcp', 'trism'],
    #
    '&SYSTEM' : ['ibrav', 'celldm(1)', 'celldm(2)', 'celldm(3)', 'celldm(4)', 'celldm(5)', 'celldm(6)', 'A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC', 'nat', 'ntyp', 'nbnd', 'nbnd_cond', 'tot_charge', 'starting_charge', 'tot_magnetization', 'starting_magnetization', 'ecutwfc', 'ecutrho', 'ecutfock', 'nr1', 'nr2', 'nr3', 'nr1s', 'nr2s', 'nr3s', 'nosym', 'nosym_evc', 'noinv', 'no_t_rev', 'force_symmorphic', 'use_all_frac', 'occupations', 'one_atom_occupations', 'starting_spin_angle', 'degauss_cond', 'nelec_cond', 'degauss', 'smearing', 'nspin', 'sic_gamma', 'pol_type', 'sic_energy', 'sci_vb', 'sci_cb', 'noncolin', 'ecfixed', 'qcutz', 'q2sigma', 'input_dft', 'ace', 'exx_fraction', 'screening_parameter', 'exxdiv_treatment', 'x_gamma_extrapolation', 'ecutvcut' 'nqx1', 'nqx2', 'nqx3', 'localization_thr', 'Hubbard_occ', 'Hubbard_alpha', 'Hubbard_beta', 'starting_ns_eigenvalue', 'dmft', 'dmft_prefix', 'ensemble_energies', 'edir', 'emaxpos', 'eopreg', 'eamp', 'angle1', 'angle2', 'lforcet', 'constrained_magnetization', 'fixed_magnetization', 'lambda', 'report', 'lspinorb', 'assume_isolated', 'esm_bc', 'esm_w', 'esm_efield', 'esm_nfit', 'lgcscf', 'gcscf_mu', 'gcscf_conv_thr', 'gcscf_beta', 'vdw_corr', 'london', 'london_s6', 'london_c6', 'london_rvdw', 'london_rcut', 'dftd3_version', 'dftd3_threebody', 'ts_vdw_econv_thr', 'ts_vdw_isolated', 'xdm', 'xdm_a1', 'xdm_a2', 'space_group', 'uniqueb', 'origin_choice', 'rhombohedral', 'zgate', 'relaxz', 'block', 'block_1', 'block_2', 'block_height', 'nextffield'],
    #
    '&ELECTRONS' : ['electron_maxstep', 'exx_maxstep', 'scf_must_converge', 'conv_thr', 'adaptive_thr', 'conv_thr_init', 'conv_thr_multi', 'mixing_mode', 'mixing_beta', 'mixing_ndim', 'mixing_fixed_ns', 'diagonalization', 'diago_thr_init', 'diago_cg_maxiter', 'diago_ppcg_maxiter', 'diago_david_ndim', 'diago_rmm_ndim', 'diago_rmm_conv', 'diago_gs_nblock', 'diago_full_acc', 'efield', 'efield_cart', 'efield_phase', 'startingpot', 'startingwfc', 'tqr', 'real_space'],
    #
    '&IONS' : ['ion_positions', 'ion_velocities', 'ion_dynamics', 'pot_extrapolation', 'wfc_extrapolation', 'remove_rigid_rot', 'ion_temperature', 'tempw', 'tolp', 'delta_t', 'nraise', 'refold_pos', 'upscale', 'bfgs_ndim', 'trust_radius_max', 'trust_radius_min', 'trust_radius_ini', 'w_1', 'w_2', 'fire_alpha_init', 'fire_falpha', 'fire_nmin', 'fire_f_inc', 'fire_f_dec', 'fire_dtmax'],
    #
    '&CELL' : ['cell_dynamics', 'press', 'wmass', 'cell_factor', 'press_conv_thr' 'cell_dofree'],
    #
    '&FCP' : ['fcp_mu', 'fcp_dynamics', 'fcp_conv_thr', 'fcp_ndiis', 'fcp_mass','fcp_velocity', 'fcp_temperature', 'fcp_tempw', 'fcp_tolp ', 'fcp_delta_t', 'fcp_nraise', 'freeze_all_atoms'],
    #
    '&RISM' : ['nsolv', 'closure', 'tempv', 'ecutsolv', 'solute_lj', 'solute_epsilon', 'solute_sigma', 'starting1d', 'starting3d', 'smear1d', 'smear3d', 'rism1d_maxstep', 'rism3d_maxstep', 'rism1d_conv_thr', 'rism3d_conv_thr', 'mdiis1d_size', 'mdiis3d_size', 'mdiis1d_step', 'mdiis3d_step', 'rism1d_bond_width', 'rism1d_dielectric', 'rism1d_molesize', 'rism1d_nproc', 'rism3d_conv_level', 'rism3d_planar_average', 'laue_nfit', 'laue_expand_right', 'laue_expand_left', 'laue_starting_right', 'laue_starting_left', 'laue_buffer_right', 'laue_buffer_left', 'laue_both_hands', 'laue_wall', 'laue_wall_z', 'laue_wall_rho', 'laue_wall_epsilon', 'laue_wall_sigma', 'laue_wall_lj6'],
    #
    'ATOMIC_SPECIES' : ['X', 'Mass_X', 'PseudoPot_X'],
    #
    'ATOMIC_POSITIONS' : ['X', 'x', 'y', 'z', 'if_pos(1)', 'if_pos(2)', 'if_pos(3)'],
    #
    'K_POINTS' : ['nks', 'xk_x', 'xk_y', 'xk_z', 'wk', 'nk1', 'nk2', 'nk3', 'sk1', 'sk2', 'sk3'],
    #
    'ADDITIONAL_K_POINTS' : ['nks_add', 'k_x', 'k_y', 'k_z', 'wk_'],
    #
    'CELL_PARAMETERS': ['v1', 'v2', 'v3'],
    #
    'CONSTRAINTS' : ['nconstr', 'constr_tol', 'constr_type', 'constr(1)', 'constr(2)', 'constr(3)', 'constr(4)', 'constr_target'],
    #
    'OCCUPATIONS': ['f_inp1', 'f_inp2'],
    #
    'ATOMIC_VELOCITIES' : ['V', 'vx', 'vy', 'vz'],
    #
    'ATOMIC_FORCES' : ['X', 'fx', 'fy', 'fz'],
    #
    'SOLVENTS' : ['X', 'Density', 'Molecule', 'X', 'Density_Left', 'Density_Right', 'Molecule'],
    #
    'HUBBARD' : ['label(1)-manifold(1)', 'u_val(1)', 'label(1)-manifold(1)', 'j0_val(1)', 'paramType(1)', 'label(1)-manifold(1)', 'paramValue(1)', 'label(I)-manifold(I)', 'u_val(I)', 'label(I)-manifold(I)', 'j0_val(I)', 'label(I)-manifold(I)', 'label(J)-manifold(J)', 'I', 'J', 'v_val(I,J)'],
}
'''
Dictionary with every possible namelist as keys, and the corresponding variables as values.
'''


def read_in(filename) -> dict:
    '''
    Reads an input `filename` from Quantum ESPRESSO,
    returning a dictionary with the input values used.
    The keys are named after the name of the corresponding variable.
    '''
    must_be_int = ['max_seconds', 'nstep', 'ibrav', 'nat', 'ntyp', 'dftd3_version', 'electron_maxstep']
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
            if var in must_be_int: # Keep ints as int
                value = int(value)
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
        atomic_species_raw = find.lines(key_species, filepath, -1, int(ntyp+1), True, True)
        # Check that there was no empty line right after the keyword:
        if atomic_species_raw:
            atomic_species_cleaned = []
            for line in atomic_species_raw:
                line = line.strip()
                if not line == '' or not line.startswith('!'):
                    atomic_species_cleaned.append(line)
            atomic_species = atomic_species_cleaned[1:]
            if len(atomic_species) > ntyp:
                atomic_species = atomic_species[:int(ntyp)]
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


def read_dir(
        folder,
        in_str:str='.in',
        out_str:str='.out'
    ) -> dict:
    '''
    Takes a `folder` containing a Quantum ESPRESSO calculation,
    and returns a dictionary containing the input parameters and output results.
    Input and output files are determined automatically,
    but must be specified with `in_str` and `out_str` if more than one file ends with `.in` or `.out`.
    '''
    input_file = file.get(folder, in_str)
    output_file = file.get(folder, out_str)
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


def read_dirs(
        directory,
        in_str:str='.in',
        out_str:str='.out',
        calc_splitter='_',
        calc_type_index=0,
        calc_id_index=1
    ) -> None:
    '''
    Calls recursively `read_dir()`, reading Quantum ESPRESSO calculations
    from all the subfolders inside the given `directory`.
    The results are saved to CSV files inside the current directory.
    Input and output files are determined automatically, but must be specified with
    `in_str` and `out_str` if more than one file ends with `.in` or `.out`.

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
            df = pd.DataFrame.from_dict(read_dir(folder, in_str, out_str))
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


def set_value(
        value,
        key:str,
        filename
    ) -> None:
    '''
    Replace the `value` of a `key` parameter in an input file with `filename`.
    If `value=''`, the parameter gets deleted.\n
    Remember to include the upper commas `'` on values that use them.\n
    Note that you must update some values before replacing others:
    'nat' before 'ATOMIC_POSITIONS', 'ntyp' before 'ATOMIC_SPECIES',
    and lattice parameters before 'CELL_PARAMETERS.
    '''
    key_uncommented = key
    key_uncommented.replace('(', r'\\(')
    key_uncommented.replace(')', r'\\)')
    key_uncommented = rf'(?!\s*!){key}'
    filepath = file.get(filename)
    input_old = read_in(filepath)
    # Check if the value is already in the file
    if not key in input_old.keys():
        _add_value(value, key, filename)
        return None
    # Check for the special values, else replace it as a regular value. ATOMIC_POSITIONS ?
    if key in ['ATOMIC_POSITIONS', 'ATOMIC_POSITIONS_old']:    
        nat = input_old['nat']
        if isinstance(value, list):
            if 'ATOMIC_SPECIES' in value[0]:
                value = value[1:]
            if len(value) != int(nat):
                raise ValueError('Update nat before updating ATOMIC_POSITIONS!')
            value = '\n'.join(value)
        if value == '':  # Remove from the file
            text.replace_line('', r'(?!\s*!)ATOMIC_POSITIONS', filepath, -1, 0, int(nat), True)
        else:
            text.replace_line(value, r'(?!\s*!)ATOMIC_POSITIONS', filepath, -1, 1, int(nat-1), True)
    # CELL_PARAMETERS ?
    elif key in ['CELL_PARAMETERS', 'CELL_PARAMETERS_old']:
        if isinstance(value, list):
            if len(value) == 4:
                if 'angstrom' in value[0] or 'bohr' in value[0]:
                    text.replace_line('', r'(?!\s*!)celldm\(\d\)\s*=', filepath, 1, 0, 0, True)
                    text.replace_line('', r'(?!\s*!)[ABC]\s*=', filepath, 1, 0, 0, True)
                    text.replace_line('', r'(?!\s*!)cos[ABC]\s*=', filepath, 1, 0, 0, True)
                elif not 'alat' in value[0]:
                    raise ValueError(f'Your CELL_PARAMETERS are invalid, please check them. Hint: card options must always be specified (angstrom, bohr, or alat). Your current CELL_PARAMETERS are:\n{value}')
                value = '\n'.join(value)
                text.replace_line(value, r'(?!\s*!)CELL_PARAMETERS', filepath, -1, 0, 3, True)
            elif len(value) == 3:
                value = '\n'.join(value)
                text.replace_line(value, r'(?!\s*!)CELL_PARAMETERS', filepath, -1, 1, 2, True)
                # We assume that lattice parameters are now in alat, so we update the title of the card
                text.replace_line('CELL_PARAMETERS alat', r'(?!\s*!)CELL_PARAMETERS', filepath, -1, 0, 0, True)
            else:
                raise ValueError('CELL_PARAMETERS must be a set of three vectors!')
        elif value == '':  # Remove from the file
            text.replace_line('', r'(?!\s*!)CELL_PARAMETERS', filepath, -1, 0, 3, True)
        else:  # Assume it was only three lines
            text.replace_line(value, r'(?!\s*!)CELL_PARAMETERS', filepath, -1, 1, 2, True)
            # We assume that lattice parameters are now in alat, so we update the title of the card
            text.replace_line('CELL_PARAMETERS alat', r'(?!\s*!)CELL_PARAMETERS', filepath, -1, 0, 0, True)
    # ATOMIC_SPECIES ?
    elif key == 'ATOMIC_SPECIES':
        ntyp = input_old['ntyp']
        if isinstance(value, list):
            if 'ATOMIC_SPECIES' in value[0]:
                value = value[1:]
            if len(value) != ntyp:
                raise ValueError('Update ntyp before updating ATOMIC_SPECIES!')
            value = '\n'.join(value)
        if value == '':  # Remove from the file
            text.replace_line('', r'(?!\s*!)ATOMIC_SPECIES', filepath, -1, 0, int(ntyp), True)
        else:
            text.replace_line(value, r'(?!\s*!)ATOMIC_SPECIES', filepath, -1, 1, int(ntyp-1), True)
    # K_POINTS ?
    elif key == 'K_POINTS':
        if value == '':  # Remove from the file
            text.replace_line('', key_uncommented, filepath, -1, 0, 1, True)
        else:
            text.replace_line(value, key_uncommented, filepath, -1, 1, 0, True)
    elif value == '':
        text.replace_line('', key_uncommented, filepath, 1, 0, 0, True)
    else:
        text.replace_line(f"  {key} = {str(value)}", key_uncommented, filepath, 1, 0, 0, True)
        # If the key is a lattice parameter, remove previous lattice parameter definitions
        if key in ['A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC']:
            text.replace_line('', r'(?!\s*!)celldm\(\d\)\s*=', filepath, 1, 0, 0, True)
            # Since lattice parameters are in alat, we update the title of the corresponding card
            text.replace_line('CELL_PARAMETERS alat', 'CELL_PARAMETERS', filepath, -1)
        elif 'celldm(' in key:
            text.replace_line('', r'(?!\s*!)[ABC]\s*=', filepath, 1, 0, 0, True)
            text.replace_line('', r'(?!\s*!)cos[ABC]\s*=', filepath, 1, 0, 0, True)
            # Since lattice parameters are in alat, we update the title of the corresponding card
            text.replace_line('CELL_PARAMETERS alat', 'CELL_PARAMETERS', filepath, -1)
    return None


def _add_value(
        value,
        key:str,
        filename
    ) -> None:
    '''
    Adds an input `value` for a `key_uncommented` that was not present before in the `filename`.
    Note that namelists must be in capital letters in yor file. Namelists must be introduced by hand.
    '''
    if value == '':  # The value was not there in the first place!
        return None
    filepath = file.get(filename)
    # CELL_PARAMETERS ?
    if key in ['CELL_PARAMETERS', 'CELL_PARAMETERS_old']:
        if isinstance(value, list):
            if len(value) == 4:
                if 'angstrom' in value[0] or 'bohr' in value[0]:
                    text.replace_line('', r'(?!\s*!)celldm\(\d\)\s*=', filepath, 1, 0, 0, True)
                    text.replace_line('', r'(?!\s*!)[ABC]\s*=', filepath, 1, 0, 0, True)
                    text.replace_line('', r'(?!\s*!)cos[ABC]\s*=', filepath, 1, 0, 0, True)
                elif not 'alat' in value[0]:
                    raise ValueError(f'Your CELL_PARAMETERS are invalid, please check them. Hint: card options must always be specified (angstrom, bohr, or alat). Your current CELL_PARAMETERS are:\n{value}')
                value = '\n'.join(value)
                text.insert_at(value, filepath, -1)
            elif len(value) == 3:
                value = '\n'.join(value)
                text.insert_at(f'CELL_PARAMETERS alat\n{value}', filepath, -1)
            else:
                raise ValueError('CELL_PARAMETERS must be a set of three vectors!')
        else:  # Assume it was only three lines
            text.insert_at(f'CELL_PARAMETERS alat\n{value}', filepath, -1)
        return None
    # ATOMIC_SPECIES?
    elif key == 'ATOMIC_SPECIES':    
        if isinstance(value, list):
            if not 'ATOMIC_SPECIES' in value[0]:
                value = value.insert(0, 'ATOMIC_SPECIES')
            value = '\n'.join(value)
        elif not value.startswith('ATOMIC_SPECIES'):
            value = 'ATOMIC_SPECIES\n' + value
        text.insert_at(value, filepath, -1)
        return None
    # ATOMIC_POSITIONS ?
    elif key in ['ATOMIC_POSITIONS', 'ATOMIC_POSITIONS_old']:    
        if isinstance(value, list):
            if not 'ATOMIC_POSITIONS' in value[0]:
                value = value.insert(0, 'ATOMIC_POSITIONS')
            value = '\n'.join(value)
        elif not value.startswith('ATOMIC_POSITIONS'):
            value = 'ATOMIC_POSITIONS\n' + value
        text.insert_at(value, filepath, -1)
        return None
    # K_POINTS ?
    elif key == 'K_POINTS':
        text.insert_at(f'K_POINTS\n{value}', filepath, -1)
        return None
    # Try with regular parameters
    done = False
    for section in pw_description.keys():
        if key in pw_description[section]:
            is_section_on_file = find.lines(section, filepath)
            if not is_section_on_file:
                _add_section(section, filepath)
            text.insert_under(f'  {key} = {str(value)}', section, filepath, 1)
            done = True
            break
    if not done:
        raise ValueError(f'Could not update the following variable: {key}. Are namelists in CAPITAL letters?')
    if key in ['A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC']:
        text.replace_line('', r'(?!\s*!)celldm\(\d\)\s*=', filepath, 1, 0, 0, True)
    elif 'celldm(' in key:
        text.replace_line('', r'(?!\s*!)[ABC]\s*=', filepath, 1, 0, 0, True)
        text.replace_line('', r'(?!\s*!)cos[ABC]\s*=', filepath, 1, 0, 0, True)
    return None


def _add_section(
        section:str,
        filename
    ) -> None:
    '''
    Adds a `section` namelist to the file with `filename`.
    The section must be in CAPITAL LETTERS, as in `&CONTROL`.
    '''
    filepath = file.get(filename)
    namelists = pw_description.keys()
    if not section in namelists:
        raise ValueError(f'{section} is not a valid namelist!')
    namelists_reversed = namelists.reverse()
    next_namelist = None
    for namelist in namelists_reversed:
        if namelist == section:
            break
        next_namelist = namelist
    next_namelist_uncommented = rf'(?!\s*!){next_namelist}'
    text.insert_under(f'{section}\n/', next_namelist_uncommented, filepath, 1, -1, True)
    return None


def scf_from_relax(
        folder:str=None,
        relax_in:str='relax.in',
        relax_out:str='relax.out'
    ) -> None:
    '''
    Create a Quantum ESPRESSO `scf.in` file from a previous relax calculation.
    If no `folder` is provided, the current working directory is used.
    The `relax_in` and `relax_out` files by default are `relax.in` and `relax.out`,
    update the names if necessary.
    '''
    # Terminal feedback
    print(f'\nthotpy.qe {version}\n'
          f'Creating Quantum ESPRESSO SCF input from previous relax calculation:\n'
          f'{relax_in}\n{relax_out}\n')
    folder = call.here(folder)
    relax_in = file.get(folder, relax_in)
    relax_out = file.get(folder, relax_out)
    data = read_dir(folder, relax_in, relax_out)
    # Create the scf.in from the previous relax.in
    scf_in = 'scf.in'
    comment = f'! Automatic SCF input made with thotpy.qe {version}. https://github.com/pablogila/ThotPy'
    file.from_template(relax_in, scf_in, comment)
    scf_in = file.get(folder, scf_in)
    # Replace CELL_PARAMETERS, ATOMIC_POSITIONS, ATOMIC_SPECIES, alat, ibrav and calculation
    atomic_species = data['ATOMIC_SPECIES']
    cell_parameters = data['CELL_PARAMETERS_out']
    atomic_positions = data['ATOMIC_POSITIONS_out']
    alat = data['alat']
    set_value(atomic_species, 'ATOMIC_SPECIES', scf_in)
    set_value(cell_parameters, 'CELL_PARAMETERS', scf_in)
    set_value(atomic_positions, 'ATOMIC_POSITIONS', scf_in)
    set_value(alat, 'celldm(1)', scf_in)
    set_value(0, 'ibrav', scf_in)
    set_value("'scf'", 'calculation', scf_in)
    # Terminal feedback
    print(f'Created input SCF file at:'
          f'{scf_in}\n')
    return None

