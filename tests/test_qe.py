import thotpy as th


folder = 'tests/samples/'


def test_read():
    ideal = {
        # relax.out
        'Energy'               : -1000.0,
        'Volume'               : 2.0,
        'Density'              : 1.0,
        'Alat'                 : 10,
        'BFGS converged'       : True,
        'BFGS failed'          : False,
        'Total force'          : 0.000001,
        'Total SCF correction' : 0.0,
        'ibrav'                : 1,
        'Runtime'              : '48m 8.71s',
        'CELL_PARAMETERS_out'  : [
            '1.000000000   0.000000000   0.000000000',
            '0.000000000   1.000000000   0.000000000',
            '0.000000000   0.000000000   1.000000000'],
        'ATOMIC_POSITIONS_out' : [
            'I                1.0000000000        0.0000000000        0.0000000000',
            'C                0.0000000000        1.0000000000        0.0000000000',
            'N                0.0000000000        0.0000000000        1.0000000000'],
        # relax.in
        'K_POINTS'             : '2 2 2 0 0 0',
        'etot_conv_thr'        : 1.0e-12,
        'max_seconds'          : 1000,
        'pseudo_dir'           : "'./pseudos/'",
        'CELL_PARAMETERS'      : [
            '2.000000000000000   0.000000000000000   0.000000000000000',
            '0.000000000000000   2.000000000000000   0.000000000000000',
            '0.000000000000000   0.000000000000000   2.000000000000000'],
        'ATOMIC_SPECIES'       : [
            'I  126.90400   I.upf',
            'N   14.00650   N.upf',
            'C   12.01060   C.upf'],
        'ATOMIC_POSITIONS'     : [
            'I   5.000000000000000   0.000000000000000   0.000000000000000',
            'C   0.000000000000000   5.000000000000000   0.000000000000000',
            'N   0.000000000000000   0.000000000000000   5.000000000000000'],
    }
    result = th.qe.read_dir(folder=folder, in_str='relax.in', out_str='relax.out')
    for key in ideal:
        assert ideal[key] == result[key]


def test_scf_from_relax():
    ideal = {
        'calculation'      : "'scf'",
        'celldm(1)'        : 10.0,
        'ibrav'            : 0,
        'occupations'      : "'fixed'",
        'conv_thr'         : 2.0e-12,
        'CELL_PARAMETERS'  : [
            '1.000000000   0.000000000   0.000000000',
            '0.000000000   1.000000000   0.000000000',
            '0.000000000   0.000000000   1.000000000'],
        'ATOMIC_SPECIES'   : [
            'I  126.90400   I.upf',
            'N   14.00650   N.upf',
            'C   12.01060   C.upf'],
        'ATOMIC_POSITIONS' : [
            'I                1.0000000000        0.0000000000        0.0000000000',
            'C                0.0000000000        1.0000000000        0.0000000000',
            'N                0.0000000000        0.0000000000        1.0000000000'],
    }
    th.qe.scf_from_relax(folder=folder)
    result = th.qe.read_in(folder + 'scf.in')
    for key in ideal:
        assert ideal[key] == result[key]
    assert 'A' not in result.keys()
    try:
        th.file.remove(folder + 'scf.in')
    except:
        pass


def test_add_atom():
    ideal_positions = [
        'I                5.0000000000        0.0000000000        0.0000000000',
        'C                0.0000000000        5.0000000000        0.0000000000',
        'N                0.0000000000        0.0000000000        5.0000000000',
        'O   0.0  0.0  0.0',
        'Cl  1.0  1.0  1.0']
    tempfile = folder + 'temp.in'
    th.file.copy(folder + 'relax.in', tempfile)
    position_1 = '  O   0.0   0.0   0.0'
    position_2 = ['Cl', 1.0, 1.0, 1.0]
    th.qe.add_atom(filename=tempfile, position=position_1)
    th.qe.add_atom(filename=tempfile, position=position_2)
    temp = th.qe.read_in(tempfile)
    nat = temp['nat']
    ntyp = temp['ntyp']
    atomic_positions = temp['ATOMIC_POSITIONS']
    assert nat == 5
    assert ntyp == 5
    for i, position in enumerate(atomic_positions):
        detected_atom = th.extract.element(position)
        ideal_atom = th.extract.element(ideal_positions[i])
        assert detected_atom == ideal_atom
        detected_coords = th.extract.coords(position)
        ideal_coords = th.extract.coords(ideal_positions[i])
        assert detected_coords == ideal_coords
    # Additional surrounding values, just in case
    assert temp['ibrav'] == 1
    assert temp['A'] == 10.0
    assert temp['ecutwfc'] == 60.0
    assert temp['input_dft'] == "'PBEsol'"
    th.file.remove(tempfile)


def test_normalize_cell_params():
    cell_params = 'CELL_PARAMETERS (alat= 10.0000)\n    1.00000000000   0.000000000   0.000000000\n   0.000000000   1.000000000   0.000000000 \n 0.000000000   0.000000000   1.0000000 '
    ideal_params = [
        'CELL_PARAMETERS alat= 10.0',
        '  1.000000000000000   0.000000000000000   0.000000000000000',
        '  0.000000000000000   1.000000000000000   0.000000000000000',
        '  0.000000000000000   0.000000000000000   1.000000000000000',]
    normalized_params = th.qe.normalize_cell_parameters(cell_params)
    assert normalized_params == ideal_params
    # Now check as a list
    cell_params = cell_params.splitlines()
    # With bohr values
    cell_params[0] = r' CELL_PARAMETERS {bohr}'
    ideal_params[0] = 'CELL_PARAMETERS bohr'
    normalized_params = th.qe.normalize_cell_parameters(cell_params)
    assert normalized_params == ideal_params
    # With armstrong values
    cell_params[0] = r' CELL_PARAMETERS {angstrom}'
    ideal_params[0] = 'CELL_PARAMETERS angstrom'
    normalized_params = th.qe.normalize_cell_parameters(cell_params)
    assert normalized_params == ideal_params


def test_normalize_atomic_positions():
    atomic_positions = " ATOMIC_POSITIONS {crystal} \n I   5.000000   0.0000000000000   0.000000000000000 \n C   0.000000000000000   5.000000000000000000   0.000000 "
    ideal_positions = [
        'ATOMIC_POSITIONS crystal',
        '  I   5.000000000000000   0.000000000000000   0.000000000000000',
        '  C   0.000000000000000   5.000000000000000   0.000000000000000']
    normalized_positions = th.qe.normalize_atomic_positions(atomic_positions)
    assert normalized_positions == ideal_positions

