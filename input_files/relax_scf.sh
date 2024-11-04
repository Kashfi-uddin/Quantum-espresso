#!/bin/bash

# Define directories
QE_INPUT_DIR="qe_input"  # Directory containing Quantum ESPRESSO input files
OUT_DIR="relax_scf"      # Output directory for calculation results
PSEUDO_DIR="./pseudos"   # Pseudopotentials directory path
NPP=2                    # Number of processors (-1 for serial)

# Create output directory if it does not exist
mkdir -p "$OUT_DIR"

# Loop through each input file in the QE_INPUT_DIR
for input_file in "$QE_INPUT_DIR"/*.in; do
    # Extract the compound name without extension
    compound_name=$(basename "$input_file" .in)

    # Create a subdirectory for the compound
    compound_dir="$OUT_DIR/$compound_name"
    mkdir -p "$compound_dir"

    # Prepare file paths for output
    vc_relax_out="$compound_dir/${compound_name}_vc-relax.out"
    scf_out="$compound_dir/${compound_name}_scf.out"

    # Read the input file and extract necessary information
    cell_parameters=()
    atomic_species=()
    atomic_positions=()
    recording_atomic_species=false
    recording_atomic_positions=false

    while IFS= read -r line; do
        if [[ "$line" == *"CELL_PARAMETERS"* ]]; then
            recording_atomic_species=false
            recording_atomic_positions=false
            continue
        fi

        if [[ "$line" == *"ATOMIC_SPECIES"* ]]; then
            recording_atomic_species=true
            continue
        fi

        if [[ "$line" == *"ATOMIC_POSITIONS"* ]]; then
            recording_atomic_positions=true
            recording_atomic_species=false
            continue
        fi

        # Only add non-empty and non-whitespace lines to atomic_species array
        if $recording_atomic_species && [[ -n "$line" && "$line" =~ [^[:space:]] ]]; then
            atomic_species+=("$line")
            continue
        fi

        if $recording_atomic_positions; then
            atomic_positions+=("$line")
            continue
        fi

        # Capture cell parameters (expecting 3 lines for 3 cell vectors)
        if [[ ${#cell_parameters[@]} -lt 3 ]]; then
            cell_parameters+=("$line")
        fi
    done < "$input_file"

    # Calculate the number of atoms
    nat=${#atomic_positions[@]}

    # Create vc-relax input file
    vc_relax_input="&CONTROL
calculation = 'vc-relax'
outdir = '$compound_dir'
prefix = '$compound_name'
pseudo_dir = '$PSEUDO_DIR'
restart_mode = 'from_scratch'
verbosity = 'high'
forc_conv_thr = 1.0d-5
/
&SYSTEM
ibrav = 0
nat = $nat
ntyp = ${#atomic_species[@]}
ecutwfc = 30
ecutrho = 240
input_dft = 'RPBE'
!lspinorb = .true.      ! Enable spin-orbit coupling
!noncolin = .true.      ! Enable non-collinear calculation
/
&ELECTRONS
conv_thr = 1.0D-8
mixing_beta = 0.7
/
&IONS
ion_dynamics = 'bfgs'
/
&CELL
cell_dynamics = 'bfgs'
press = 0.0
cell_factor = 2.0
/
ATOMIC_SPECIES
$(printf "%s\n" "${atomic_species[@]}")

ATOMIC_POSITIONS (alat)
$(printf "%s\n" "${atomic_positions[@]}")

CELL_PARAMETERS (angstrom)
$(printf "%s\n" "${cell_parameters[@]}")

K_POINTS automatic
4 4 4 0 0 0
"

    # Write the vc-relax input file
    vc_relax_input_file="$compound_dir/${compound_name}_vc-relax.in"
    echo "$vc_relax_input" > "$vc_relax_input_file"

    # Run vc-relax calculation
    if [[ $NPP -eq -1 ]]; then
        pw.x -in "$vc_relax_input_file" > "$vc_relax_out"
    else
        mpirun -np $NPP pw.x -in "$vc_relax_input_file" | tee "$vc_relax_out"
    fi

    echo "Variable-cell relaxation completed for $compound_name."

    # Extract relaxed structure from vc-relax output
    relaxed_positions=""
    relaxed_cell=""
    in_final_section=false

# Read the relax.out file
while IFS= read -r line; do
    # Check for the start of the final coordinates section
    if [[ "$line" == *"Begin final coordinates"* ]]; then
        in_final_section=true
    elif [[ "$line" == *"End final coordinates"* ]]; then
        in_final_section=false
    fi

   # If we're in the final section, extract CELL_PARAMETERS and ATOMIC_POSITIONS
    if $in_final_section; then
        if [[ "$line" == *"CELL_PARAMETERS (angstrom)"* ]]; then
            # Read the cell parameters
            read -r line  # Skip the line with CELL_PARAMETERS
            for i in {1..3}; do
                final_cell_parameters+="$line\n"
                read -r line
            done
        elif [[ "$line" == *"ATOMIC_POSITIONS (alat)"* ]]; then
            # Read the atomic positions
            read -r pos_line  # Skip the line with ATOMIC_POSITIONS
            for ((i=1; i<=nat; i++)); do
                final_positions+="$pos_line\n"
                read -r pos_line
            done
        fi
    fi
done < "$vc_relax_out"

# Remove the trailing newline characters for formatting
#final_positions="${final_positions//$'\n'/\\n}"
#final_cell_parameters="${final_cell_parameters//$'\n'/\\n}"

# Output the extracted final results for verification
#echo -e "Final CELL_PARAMETERS (angstrom):\n$final_cell_parameters"
#echo -e "Final ATOMIC_POSITIONS (alat):\n$final_positions"

    # Create SCF input file using relaxed positions and cell parameters
scf_input="&CONTROL
calculation = 'scf'
outdir = '$compound_dir'
prefix = '$compound_name'
pseudo_dir = '$PSEUDO_DIR'
restart_mode = 'from_scratch'
verbosity = 'high'
/
&SYSTEM
ibrav = 0
nat = $nat
ntyp = ${#atomic_species[@]}
ecutwfc = 30
ecutrho = 240
nbnd = 24
occupations = 'fixed'
input_dft = 'RPBE'
!lspinorb = .true.      ! Enable spin-orbit coupling
!noncolin = .true.      ! Enable non-collinear calculation
/
&ELECTRONS
conv_thr = 1.0D-8
mixing_beta = 0.7
/
ATOMIC_SPECIES
$(printf "%s\n" "${atomic_species[@]}")

CELL_PARAMETERS (angstrom)
$(printf "$final_cell_parameters")

ATOMIC_POSITIONS (alat)
$(printf "$final_positions")

K_POINTS automatic
5 5 5 0 0 0
"

    # Write the SCF input file
    scf_input_file="$compound_dir/${compound_name}_scf.in"
    echo "$scf_input" > "$scf_input_file"

    # Run SCF calculation
    if [[ $NPP -eq -1 ]]; then
        pw.x -in "$scf_input_file" > "$scf_out"
    else
        mpirun -np $NPP pw.x -in "$scf_input_file" | tee "$scf_out"
    fi

    echo "SCF calculation completed for $compound_name."
done

echo "All calculations completed."
