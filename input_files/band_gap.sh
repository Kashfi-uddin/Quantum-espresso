#!/bin/bash

# Define the directory containing compound subdirectories
OUT_DIR="relax_scf"
OUTPUT_CSV="band_gaps.csv"

# Create the CSV file and add the header
echo "Compound,Band Gap (eV)" > "$OUTPUT_CSV"

# Loop through each compound subdirectory in the output directory
for compound_dir in "$OUT_DIR"/*/; do
    # Get the compound name from the directory name
    compound_name=$(basename "$compound_dir")

    # Define the path to the scf.out file
    scf_out_file="${compound_dir}/${compound_name}_scf.out"

    # Check if the scf.out file exists
    if [[ -f "$scf_out_file" ]]; then
        # Extract the line containing "highest occupied, lowest unoccupied level"
        energy_line=$(grep "highest occupied, lowest unoccupied level (ev):" "$scf_out_file")

        # Check if the line was found
        if [[ -n "$energy_line" ]]; then
            # Extract the two energy levels (highest occupied and lowest unoccupied) using awk
            highest_occupied=$(echo "$energy_line" | awk '{print $7}')
            lowest_unoccupied=$(echo "$energy_line" | awk '{print $8}')

            # Calculate the band gap
            band_gap=$(echo "$lowest_unoccupied - $highest_occupied" | bc -l)

            # Append the result to the CSV file
            echo "$compound_name,$band_gap" >> "$OUTPUT_CSV"
            echo "Band gap for $compound_name calculated and recorded."
        else
            echo "Energy levels not found in $scf_out_file for $compound_name."
        fi
    else
        echo "No scf.out file found for $compound_name."
    fi
done

echo "Band gap extraction completed. Results saved in $OUTPUT_CSV."

