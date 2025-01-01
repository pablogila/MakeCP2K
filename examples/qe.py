import thotpy as th

folder = th.call.here()

dictionary = th.qe.read_dir(folder)

print(f'\nK_points = {dictionary['K_POINTS']}\n')

#print(f'\ncelldm(1) = {df['celldm(1)']}\n')

print(f'\nA = {dictionary['A']}\n')

print(f'\nalat = {dictionary['alat']}\n')

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
print('')