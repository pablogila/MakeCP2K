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
- `add_atom()`
- `normalize_cell_parameters()`
- `normalize_atomic_positions()`
- `normalize_atomic_species()`
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
import maatpy as mt


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


def read_in(filepath) -> dict:
    '''
    Reads an input `filepath` from Quantum ESPRESSO,
    returning a dictionary with the input values used.
    The keys are named after the name of the corresponding variable.
    '''
    must_be_int = ['max_seconds', 'nstep', 'ibrav', 'nat', 'ntyp', 'dftd3_version', 'electron_maxstep']
    file_path = file.get(filepath)
    data = {}
    lines = find.lines(file_path, '=')
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
    # K_POINTS
    k_points = find.lines(file_path, r'(?!\s*!)(k_points|K_POINTS)', -1, 1, True, True)
    if k_points:
        k_points = k_points[1].strip()
        data['K_POINTS'] = k_points
    # ATOMIC_SPECIES
    key_species = r'(?!\s*!)(ATOMIC_SPECIES|atomic_species)'
    atomic_species = None
    if data['ntyp']:
        ntyp = data['ntyp']
        atomic_species_raw = find.lines(file_path, key_species, -1, int(ntyp+1), True, True)
        atomic_species = normalize_atomic_species(atomic_species_raw)
        if atomic_species:  # Just a check for clueless people
            if len(atomic_species) != ntyp:
                print(f'WARNING: ntyp={ntyp}, len(ATOMIC_SCPECIES)={len(atomic_species)}')
    else: # We assume species go before
        key_species_end = r"(?!\s*!)(ATOMIC_POSITIONS|atomic_positions|CELL_PARAMETERS|cell_parameters)" 
        atomic_species = find.between(file_path, key_species, key_species_end, False, 1, True)
        atomic_species = normalize_atomic_species(atomic_species_raw)
    if atomic_species:
        data['ATOMIC_SPECIES'] = atomic_species
    # CELL_PARAMETERS. Let's take some extra lines just in case there were empty or commented lines in between.
    cell_parameters_raw = find.lines(file_path, r'(?!\s*!)(cell_parameters|CELL_PARAMETERS)', -1, 4, True, True)
    if cell_parameters_raw:
        cell_parameters = normalize_cell_parameters(cell_parameters_raw)
        if cell_parameters:
            # extract a possible alat from CELL_PARAMETERS
            alat = extract.number(cell_parameters[0])
            if alat:  # This overwrites any possible celldm(1) previously defined!
                data['celldm(1)'] = alat
            cell_parameters[0] = 'CELL_PARAMETERS alat'
            data['CELL_PARAMETERS'] = cell_parameters
    # ATOMIC_POSITIONS. We assume nat is correct.
    if data['nat']:
        nat = data['nat']
        atomic_positions_raw = find.lines(file_path, r'(?!\s*!)(atomic_positions|ATOMIC_POSITIONS)', -1, int(nat+1), True, True)
        if atomic_positions_raw:
            atomic_positions = normalize_atomic_positions(atomic_positions_raw)
            if len(atomic_positions) != (nat + 1):
                print("WARNING:  len(nat) != len (ATOMIC_POSITIONS)")
            data['ATOMIC_POSITIONS'] = atomic_positions
    else:
        print("WARNING: 'nat' is missing, so no ATOMIC_POSITIONS were obtained!")
    return data


def read_out(filepath) -> dict:
    '''
    Reads an output `filepath` from Quantum ESPRESSO,
    returning a dict with the following keys:\n
    `'Energy'` (Ry), `'Total force'` (float), `'Total SCF correction'` (float),
    `'Runtime'` (str), `'JOB DONE'` (bool), `'BFGS converged'` (bool), `'BFGS failed'` (bool),
    `'Maxiter reached'` (bool), `'Error'` (str), `'Success'` (bool), `'CELL_PARAMETERS_out'` (list of str), `'ATOMIC_POSITIONS_out'` (list of str), `'Alat'` (bohr), `'Volume'` (a.u.^3), `'Density'` (g/cm^3).\n
    Note that these output keys start with a **C**apital letter.
    '''
    file_path = file.get(filepath)

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

    energy_line          = find.lines(file_path, energy_key, -1)
    force_line           = find.lines(file_path, force_key, -1)
    time_line            = find.lines(file_path, time_key, -1)
    job_done_line        = find.lines(file_path, job_done_key, -1)
    bfgs_converged_line  = find.lines(file_path, bfgs_converged_key, -1)
    bfgs_failed_line     = find.lines(file_path, bfgs_failed_key, -1)
    maxiter_reached_line = find.lines(file_path, maxiter_reached_key, -1)
    error_line           = find.lines(file_path, error_key, -1, 1, True)

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
    volume = None
    density = None
    coordinates_raw = find.between(file_path, 'Begin final coordinates', 'End final coordinates', False, -1, False)
    if coordinates_raw:
        coordinates_raw = coordinates_raw.splitlines()
        append_cell = False
        append_positions = False
        cell_parameters_raw = []
        atomic_positions_raw = []
        for line in coordinates_raw:
            line = line.strip()
            if cell_parameters_key in line:
                append_cell = True
                append_positions = False
            elif atomic_positions_key in line:
                append_cell = False
                append_positions = True
            elif 'volume' in line:
                volume = extract.number(line, 'volume')
            elif 'density' in line:
                density = extract.number(line, 'density')
            if line == '' or line.startswith('!'):
                continue
            if append_cell:
                cell_parameters_raw.append(line)
            elif append_positions:
                atomic_positions_raw.append(line)
        cell_parameters = normalize_cell_parameters(cell_parameters_raw)
        atomic_positions = normalize_atomic_positions(atomic_positions_raw)
        if cell_parameters:
            if 'alat' in cell_parameters[0]:
                alat = extract.number(cell_parameters[0], 'alat')

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
        'Alat'                  : alat,
        'Volume'                : volume,
        'Density'               : density,
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
        filepath,
        key:str,
        value
    ) -> None:
    '''
    Replace the `value` of a `key` parameter in an input `filepath`.
    If `value=''`, the parameter gets deleted.\n
    Remember to include the upper commas `'` on values that use them.\n
    Updating 'ATOMIC_POSITIONS' updates 'nat' automatically,
    and updating 'ATOMIC_SPECIES' updates 'ntyp'.
    '''
    key_uncommented = key
    key_uncommented = key_uncommented.replace('(', r'\(')
    key_uncommented = key_uncommented.replace(')', r'\)')
    key_uncommented = rf'(?!\s*!){key_uncommented}'
    file_path = file.get(filepath)
    input_old = read_in(file_path)
    # Do we have to include the value from scratch?
    if not key in input_old.keys():
        _add_value(file_path, key, value)
        return None
    # K_POINTS ?
    if key == 'K_POINTS':
        if value == '':  # Remove from the file
            text.replace_line(file_path, key_uncommented, '', -1, 0, 1, True)
        else:
            text.replace_line(file_path, key_uncommented, value, -1, 1, 0, True)
        return None
    # ATOMIC_SPECIES ?
    elif key == 'ATOMIC_SPECIES':
        ntyp = input_old['ntyp']
        atomic_species = normalize_atomic_species(value)
        atomic_species_str = '\n'.join(atomic_species)
        if value == '':  # Remove from the file
            text.replace_line(file_path, r'(?!\s*!)ATOMIC_SPECIES', '', -1, 0, int(ntyp), True)
        else:
            text.replace_line(file_path, r'(?!\s*!)ATOMIC_SPECIES', atomic_species_str, -1, 1, int(ntyp-1), True)
        new_ntyp = len(atomic_species)
        if new_ntyp != ntyp:
            text.replace_line(file_path, r'(?!\s*!)ntyp\s*=', f'  ntyp = {new_ntyp}', 1, 0, 0, True)
        return None
    # CELL_PARAMETERS ?
    elif key in ['CELL_PARAMETERS', 'CELL_PARAMETERS_out']:
        cell_parameters = normalize_cell_parameters(value)
        cell_parameters_str = '\n'.join(cell_parameters)
        if value == '':  # Remove from the file
            text.replace_line(file_path, r'(?!\s*!)CELL_PARAMETERS', '', -1, 0, 3, True)
        else:
            text.replace_line(file_path, r'(?!\s*!)CELL_PARAMETERS', cell_parameters_str, -1, 0, 3, True)
            # Now, if there were units there, overwrite previous definitions
            if 'angstrom' in cell_parameters[0] or 'bohr' in cell_parameters[0]:
                text.replace_line(file_path, r'(?!\s*!)celldm\(\d\)\s*=', '', 1, 0, 0, True)
                text.replace_line(file_path, r'(?!\s*!)[ABC]\s*=', '', 1, 0, 0, True)
                text.replace_line(file_path, r'(?!\s*!)cos[ABC]\s*=', '', 1, 0, 0, True)
            elif 'alat' in cell_parameters[0]:
                alat = extract.number(cell_parameters[0])
                if alat:
                    text.replace_line(file_path, r'(?!\s*!)CELL_PARAMETERS', 'CELL_PARAMETERS alat', -1, 0, 0, True)
                    text.replace_line(file_path, r'(?!\s*!)celldm\(\d\)\s*=', '', 1, 0, 0, True)
                    text.insert_under(file_path, r'(?!\s*!)&SYSTEM', f'  celldm(1) = {alat}', 1, 0, True)
        return None
    # ATOMIC_POSITIONS ?
    elif key in ['ATOMIC_POSITIONS', 'ATOMIC_POSITIONS_out']:    
        nat = input_old['nat']
        new_nat = None
        atomic_positions = normalize_atomic_positions(value)
        new_nat = len(atomic_positions) - 1
        atomic_positions_str = '\n'.join(atomic_positions)
        if value == '':  # Remove from the file
            text.replace_line(file_path, r'(?!\s*!)ATOMIC_POSITIONS', '', -1, 0, int(nat), True)
        else:
            text.replace_line(file_path, r'(?!\s*!)ATOMIC_POSITIONS', atomic_positions_str, -1, 0, int(nat), True)
        if new_nat != nat:
            text.replace_line(file_path, r'(?!\s*!)nat\s*=', f'  nat = {new_nat}', 1, 0, 0, True)
        return None
    # The value seems single-lined, do we want to delete it?
    elif value == '':
        text.replace_line(file_path, key_uncommented, '', 1, 0, 0, True)
        return None
    # Update a single-line value
    else:
        must_be_int = ['max_seconds', 'nstep', 'ibrav', 'nat', 'ntyp', 'dftd3_version', 'electron_maxstep']
        if key in must_be_int:
            value = int(value)
        text.replace_line(file_path, key_uncommented, f"  {key} = {str(value)}", 1, 0, 0, True)
        # If the key is a lattice parameter, remove previous lattice parameter definitions
        if key in ['A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC']:
            text.replace_line(file_path, r'(?!\s*!)celldm\(\d\)\s*=', '', 1, 0, 0, True)
            text.replace_line(file_path, r'(?!\s*!)CELL_PARAMETERS', 'CELL_PARAMETERS alat', -1)
        elif 'celldm(' in key:
            text.replace_line(file_path, r'(?!\s*!)[ABC]\s*=', '', 1, 0, 0, True)
            text.replace_line(file_path, r'(?!\s*!)cos[ABC]\s*=', '', 1, 0, 0, True)
            text.replace_line(file_path, r'(?!\s*!)CELL_PARAMETERS', 'CELL_PARAMETERS alat', -1)
    return None


def _add_value(
        filepath,
        key:str,
        value
    ) -> None:
    '''
    Adds an input `value` for a `key_uncommented` that was not present before in the `filepath`.
    Note that namelists must be in capital letters in yor file. Namelists must be introduced by hand.
    '''
    if value == '':  # We are trying to delete a value that is not present!
        return None
    file_path = file.get(filepath)
    old_values = read_in(file_path)
    # K_POINTS ?
    if key == 'K_POINTS':
        k_points = f'K_POINTS\n{value}'
        text.insert_at(file_path, k_points, -1)
        return None
    # ATOMIC_SPECIES ?
    elif key == 'ATOMIC_SPECIES':    
        atomic_species = normalize_atomic_species(value)
        new_ntyp = len(atomic_species)
        atomic_species_str = '\n'.join(atomic_species)
        atomic_species_str = 'ATOMIC_SPECIES\n' + atomic_species_str
        text.insert_at(file_path, atomic_species_str, -1)
        if old_values['ntyp'] != new_ntyp:
            text.replace_line(file_path, r'(?!\s*!)ntyp\s*=', f'ntyp = {new_ntyp}', 1, 0, 0, True)
        return None
    # CELL_PARAMETERS ?
    elif key in ['CELL_PARAMETERS', 'CELL_PARAMETERS_old']:
        cell_parameters = normalize_cell_parameters(value)
        cell_parameters_str = '\n'.join(cell_parameters)
        text.insert_at(file_path, cell_parameters_str, -1)
        if 'angstrom' in cell_parameters[0] or 'bohr' in cell_parameters[0]:
            text.replace_line(file_path, r'(?!\s*!)celldm\(\d\)\s*=', '', 1, 0, 0, True)
            text.replace_line(file_path, r'(?!\s*!)[ABC]\s*=', '', 1, 0, 0, True)
            text.replace_line(file_path, r'(?!\s*!)cos[ABC]\s*=', '', 1, 0, 0, True)
        elif 'alat' in cell_parameters[0]:
            alat = extract.number(cell_parameters[0])
            if alat:
                text.replace_line(file_path, r'(?!\s*!)CELL_PARAMETERS', 'CELL_PARAMETERS alat', -1, 0, 0, True)
                text.replace_line(file_path, r'(?!\s*!)celldm\(\d\)\s*=', f'celldm(1) = {alat}', 1, 0, 0, True)
        return None
    # ATOMIC_POSITIONS ?
    elif key in ['ATOMIC_POSITIONS', 'ATOMIC_POSITIONS_old']:
        atomic_positions = normalize_atomic_positions(value)
        new_nat = len(atomic_positions) - 1
        atomic_positions_str = '\n'.join(atomic_positions)
        text.insert_at(file_path, atomic_positions_str, -1)
        if old_values['nat'] != new_nat:
            text.replace_line(file_path, r'(?!\s*!)nat\s*=', f'nat = {new_nat}', 1, 0, 0, True)
        return None
    # Try with regular parameters
    done = False
    for section in pw_description.keys():
        if key in pw_description[section]:
            is_section_on_file = find.lines(file_path, section)
            if not is_section_on_file:
                _add_section(file_path, section)
            text.insert_under(file_path, section, f'  {key} = {str(value)}', 1)
            done = True
            break
    if not done:
        raise ValueError(f'Could not update the following variable: {key}. Are namelists in CAPITAL letters?')
    if key in ['A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC']:
        text.replace_line(file_path, r'(?!\s*!)celldm\(\d\)\s*=', '', 1, 0, 0, True)
        text.replace_line(file_path, r'(?!\s*!)CELL_PARAMETERS', 'CELL_PARAMETERS alat', -1, 0, 0, True)
    elif 'celldm(' in key:
        text.replace_line(file_path, r'(?!\s*!)[ABC]\s*=', '', 1, 0, 0, True)
        text.replace_line(file_path, r'(?!\s*!)cos[ABC]\s*=', '', 1, 0, 0, True)
        text.replace_line(file_path, r'(?!\s*!)CELL_PARAMETERS', 'CELL_PARAMETERS alat', -1, 0, 0, True)
    return None


def _add_section(
        filepath,
        section:str
    ) -> None:
    '''
    Adds a `section` namelist to the `filepath`.
    The section must be in CAPITAL LETTERS, as in `&CONTROL`.
    '''
    file_path = file.get(filepath)
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
    text.insert_under(file_path, next_namelist_uncommented, f'{section}\n/', 1, -1, True)
    return None


def add_atom(filepath, position) -> None:
    '''
    Adds an atom in a given `filepath` at a specified `position`.
    Position must be a string or a list, as follows:\n
    `"specie:str float float float"` or `[specie:str, float, float, float]`\n
    This method updates automatically 'ntyp' and 'nat'.
    '''
    if isinstance(position, list):
        if not len(position) == 4 or not isinstance(position[0], str):
            raise ValueError('If your atomic position is a list, it must contain the atom type and the three coords, as in [str, str/float, str/float, str/float]')
        new_atom = position[0] + '   ' + str(position[1]) + '   ' + str(position[2]) + '   ' + str(position[3])
    elif isinstance(position, str):
        new_atom = position.strip()
    else:
        raise ValueError(f'The specified position must be a list of size 4 (atom type and three coordinates) or an equivalent string! Your position was:\n{coords}')
    # Let's check that our new_atom makes sense
    atom = extract.element(new_atom)
    coords = extract.coords(new_atom)
    if not atom:
        raise ValueError(f'The specified position does not contain an atom at the beginning! Your position was:\n{coords}')
    if len(coords) < 3:
        raise ValueError(f'Your position has len(coordinates) < 3, please check it.\nYour position was: {position}\nCoordinates detected: {coords}')
    if len(coords) > 3:
        coords = coords[:3]
    new_atom = atom + '   ' + str(coords[0]) + '   ' + str(coords[1]) + '   ' + str(coords[2])
    # Get the values from the file
    values = read_in(filepath)
    nat = values['nat']
    ntyp = values['ntyp']
    atomic_positions = values['ATOMIC_POSITIONS']
    nat += 1
    atomic_positions.append(new_atom)
    set_value(filepath, 'ATOMIC_POSITIONS', atomic_positions)
    set_value(filepath, 'nat', nat)
    # We might have to update ATOMIC_SPECIES!
    atomic_species = values['ATOMIC_SPECIES']
    is_atom_missing = True
    for specie in atomic_species:
        if atom == extract.element(specie):
            is_atom_missing = False
            break
    if is_atom_missing:
        mass = mt.atom[atom].mass
        atomic_species.append(f'{atom}   {mass}   {atom}.upf')
        ntyp += 1
        set_value(filepath=filepath, key='ATOMIC_SPECIES', value=atomic_species)
        set_value(filepath=filepath, key='ntyp', value=ntyp)
    return None


def normalize_cell_parameters(params) -> list:
    '''
    Takes a `params` string or a list of strings with the cell parameters
    and possibly some additional rogue lines, and returns a list os size 4,
    with the "CELL_PARAMETERS {alat|bohr|angstrom}" on list[0],
    followed by the three coordinates.
    '''
    if params == None:
        return None
    if isinstance(params, str):  # Convert to a list
        params = params.splitlines()
    if not isinstance(params, list):
        raise ValueError(f'The provided cell parameters must be a list or a string! Yours was:\n{params}')
    # Clean it
    cell_key = 'cell_parameters'
    stop_keys = ['atomic_species','atomic_positions']
    header = ''
    cell_parameters = []
    for line in params:
        line = line.strip()
        if line == '' or line.startswith('!'):
            continue
        if any(key in line.lower() for key in stop_keys):
            break
        if cell_key in line.lower():
            header = line
            continue
        coords = extract.coords(line)
        if len(coords) < 3:
            raise ValueError(f'Each CELL_PARAMETER must have three coordinates! Yours was:\n{params}\nDetected coordinates were:\n{coords}')
        if len(coords) > 3:
            coords = coords[:3]
        new_line = f"  {coords[0]:.15f}   {coords[1]:.15f}   {coords[2]:.15f}"
        cell_parameters.append(new_line)
    # Check the header
    if 'bohr' in header.lower():
        header = 'CELL_PARAMETERS bohr'
    elif 'angstrom' in header.lower():
        header = 'CELL_PARAMETERS angstrom'
    elif 'alat' in header.lower():
        alat = extract.number(header, 'alat')
        header = 'CELL_PARAMETERS alat'
        if alat:
            header = f'CELL_PARAMETERS alat= {alat}'
    elif not header:
        header = 'CELL_PARAMETERS alat'
    else:
        raise ValueError(f'CELL_PARAMETERS must be in alat, bohr or angstrom! Yours was:\n{header}')
    cell_parameters.insert(0, header)
    return cell_parameters


def normalize_atomic_positions(positions) -> list:
    '''
    Takes a `positions` string or a list of strings with the atomic positions
    and possibly some additional rogue lines, and returns a list with the atomic positions,
    with the "ATOMIC_POSITIONS {alat|bohr|angstrom|crystal|crystal_sg}" on list[0],
    followed by the coordinates.
    '''
    if positions == None:
        return None
    if isinstance(positions, str):  # Convert to a list
        positions = positions.splitlines()
    if not isinstance(positions, list):
        raise ValueError(f'The provided atomic positions must be a list or a string! Yours was:\n{positions}')
    # Clean it
    pos_key = 'atomic_positions'
    stop_keys = ['atomic_species','cell_parameters']
    header = ''
    atomic_positions = []
    for line in positions:
        line = line.strip()
        if line == '' or line.startswith('!'):
            continue
        if any(key in line.lower() for key in stop_keys):
            break
        if pos_key in line.lower():
            header = line
            continue
        atom = extract.element(line)
        if not atom:
            raise ValueError(f'Atoms must be defined as the atom (H, He, Na...) or the isotope (H2, He4...)! Yours was:\n{line}')
        coords = extract.coords(line)
        if len(coords) < 3:
            raise ValueError(f'Each ATOMIC_POSITION must have at least three coordinates! Yours contained the line:\n{line}\nDetected coordinates were:\n{coords}')
        if len(coords) > 6:  # Including optional parameters
            coords = coords[:6]
        new_line = f"  {atom}   {coords[0]:.15f}   {coords[1]:.15f}   {coords[2]:.15f}"
        atomic_positions.append(new_line)
    # Check the header
    if 'bohr' in header.lower():
        header = 'ATOMIC_POSITIONS bohr'
    elif 'angstrom' in header.lower():
        header = 'ATOMIC_POSITIONS angstrom'
    elif 'alat' in header.lower():
        header = 'ATOMIC_POSITIONS alat'
    elif 'crystal_sg' in header.lower():
        header = 'ATOMIC_POSITIONS crystal_sg'
    elif 'crystal' in header.lower():
        header = 'ATOMIC_POSITIONS crystal'
    elif not header:
        header = 'ATOMIC_POSITIONS crystal'
    else:
        raise ValueError(f'ATOMIC_POSITIONS must be in alat, bohr, angstrom, crystal or crystal_sg. Yours was:\n{header}')
    atomic_positions.insert(0, header)
    return atomic_positions


def normalize_atomic_species(species) -> list:
    '''
    Takes a `species` string or a list of strings with the atomic species
    and possibly some additional rogue lines, and returns a list with the atomic species
    (without the ATOMIC_SPECIES header!).
    '''
    if species == None:
        return None
    if isinstance(species, str):  # Convert to a list
        species = species.splitlines()
    if not isinstance(species, list):
        raise ValueError(f'The provided atomic species must be a list or a string! Yours was:\n{species}')
    # Clean it
    stop_keys = ['atomic_positions','cell_parameters']
    atomic_species = []
    for line in species:
        line = line.strip()
        if line == '' or line.startswith('!'):
            continue
        if any(key in line.lower() for key in stop_keys):
            break
        if 'atomic_species' in line.lower():
            continue
        atom = None
        atom = extract.element(line)
        if not atom:
            raise ValueError(f'Atom must be specified at the beginning! Your line was:\n{line}')
        mass_list = extract.coords(line)
        if len(mass_list) == 1:
            mass = mass_list[0]
        else:  # Is the mass missing?
            raise ValueError(f'Mass is not properly specified: {line}')
        # Get the pseudo in the third position
        line_split = line.split()
        if len(line_split) < 3:
            raise ValueError(f'Does the ATOMIC_SPECIES contain the pseudopotential? Your line was:\n{line}')
        pseudopotential = line_split[2]
        full_line = f"  {atom}   {mass}   {pseudopotential}"
        atomic_species.append(full_line)
    return atomic_species


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
    folder_path = folder
    if not folder_path:
        folder_path = os.getcwd()
    relax_in = file.get(folder_path, relax_in)
    relax_out = file.get(folder_path, relax_out)
    data = read_dir(folder_path, relax_in, relax_out)
    # Create the scf.in from the previous relax.in
    scf_in = os.path.join(folder_path, 'scf.in')
    comment = f'! Automatic SCF input made with thotpy.qe {version}. https://github.com/pablogila/ThotPy'
    file.from_template(relax_in, scf_in, None, comment)
    scf_in = file.get(folder_path, scf_in)
    # Replace CELL_PARAMETERS, ATOMIC_POSITIONS, ATOMIC_SPECIES, alat, ibrav and calculation
    atomic_species = data['ATOMIC_SPECIES']
    cell_parameters = data['CELL_PARAMETERS_out']
    atomic_positions = data['ATOMIC_POSITIONS_out']
    alat = data['Alat']
    set_value(scf_in, 'ATOMIC_SPECIES', atomic_species)
    set_value(scf_in, 'CELL_PARAMETERS', cell_parameters)
    set_value(scf_in, 'ATOMIC_POSITIONS', atomic_positions)
    set_value(scf_in, 'celldm(1)', alat)
    set_value(scf_in, 'ibrav', 0)
    set_value(scf_in, 'calculation', "'scf'")
    # Terminal feedback
    print(f'Created input SCF file at:'
          f'{scf_in}\n')
    return None

