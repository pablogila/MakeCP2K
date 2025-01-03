import thotpy as th

# Run in this folder
th.call.here()

# Create the inputs
th.phonopy.make()

# Once ready, sbatch the calculations. Since testing=True, it will only print the files to sbatch.
th.phonopy.sbatch(testing=True)
