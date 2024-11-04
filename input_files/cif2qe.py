import os
from ase.io import read
from ase.data import atomic_masses, chemical_symbols

# Define directories
cif_dir = "cif"  # Directory containing CIF files
qe_output_dir = "qe_input"  # Output directory for Quantum ESPRESSO input files

# Create the output directory if it does not exist
if not os.path.exists(qe_output_dir):
    os.makedirs(qe_output_dir)

# Loop through each CIF file in the directory
for cif_file in os.listdir(cif_dir):
    if cif_file.endswith(".cif"):
        compound_name = os.path.splitext(cif_file)[0]  # Extract compound name without .cif extension

        # Load the CIF structure
        cif_path = os.path.join(cif_dir, cif_file)
        structure = read(cif_path)

        # Extract cell parameters and atomic positions
        cell = structure.get_cell()
        positions = structure.get_scaled_positions()
        symbols = structure.get_chemical_symbols()

        # Prepare Quantum ESPRESSO input filename
        qe_filename = f"{compound_name}.in"
        qe_filepath = os.path.join(qe_output_dir, qe_filename)

        # Create the ATOMIC_SPECIES block while avoiding duplicates
        atomic_species = {}
        species_list = []  # To maintain the order of first appearance
        for symbol in symbols:
            if symbol not in atomic_species:
                atomic_mass = atomic_masses[chemical_symbols.index(symbol)]  # Get atomic mass
                pseudo_file = f"{symbol}_pbe.UPF"  # Replace with actual pseudopotential filenames
                atomic_species[symbol] = f"{symbol}  {atomic_mass:.4f}  {pseudo_file}\n"
                species_list.append(symbol)  # Record the order of unique symbols

        # Write Quantum ESPRESSO input file
        with open(qe_filepath, "w") as f:
            f.write("CELL_PARAMETERS (angstrom)\n")
            for vector in cell:
                f.write(f"  {vector[0]:.6f} {vector[1]:.6f} {vector[2]:.6f}\n")

            f.write("\nATOMIC_SPECIES\n")
            # Write atomic species in the order of their first appearance
            for symbol in species_list:
                f.write(atomic_species[symbol])

            f.write("\nATOMIC_POSITIONS (crystal)\n")
            for symbol, pos in zip(symbols, positions):
                f.write(f"  {symbol} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")

        print(f"Generated Quantum ESPRESSO input file: {qe_filename}")

print("All CIF files have been processed.")
