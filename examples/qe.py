import thotpy as th

# Run in this folder
folder = th.call.here()

# Read the entire directory
dictionary = th.qe.read_dir(folder)

# Print some values. Full list in the documentation of the qe module.
print(f'\nK_POINTS = {dictionary['K_POINTS']}\n')
#print(f'\ncelldm(1) = {df['celldm(1)']}\n')
print(f'\nA = {dictionary['A']}\n')
print(f'\nAlat = {dictionary['Alat']}\n')
cell_in = dictionary['CELL_PARAMETERS']
print('\ncell in:\n')
for cell in cell_in:
    print(cell)
print('')
cell_out = dictionary['CELL_PARAMETERS_out']
print('\ncell out:\n')
for cell in cell_out:
    print(cell)
print('')
atomic_positions_in = dictionary['ATOMIC_POSITIONS']
print('\natomic positions:\n')
for position in atomic_positions_in:
    print(position)
print('')
atomic_positions_out = dictionary['ATOMIC_POSITIONS_out']
print('\natomic positions out:\n')
for position in atomic_positions_out:
    print(position)
print()

# You can also create a new SCF input calculation from the relaxed structure!
# th.qe.scf_from_relax()

