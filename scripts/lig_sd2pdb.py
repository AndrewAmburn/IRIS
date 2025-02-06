#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import sys

def convert_sd_to_pdb(base_directory):
    # Traverse the base directory and its subfolders
    for root, dirs, files in os.walk(base_directory):
        for file in files:
            # Check for files ending with "lig.sd" and generate the corresponding ".pdb" file
            if file.endswith("lig.sd"):
                input_file = os.path.join(root, file)
                output_file = os.path.join(root, file.replace("lig.sd", "lig.pdb"))

                print(f"Converting {input_file} to {output_file}...")

                # Build the command
                command = ["obabel", "-i", "sd", input_file, "-O", output_file]

                try:
                    # Run the command and capture output
                    result = subprocess.run(command, check=True, capture_output=True, text=True)
                    print(f"Conversion complete for {input_file}. Output:\n{result.stdout}")
                except subprocess.CalledProcessError as e:
                    print(f"Error converting {input_file}: {e.stderr}")
                except OSError as e:
                    print(f"Error processing file {input_file}: {e}")

def main():
    # Ensure the correct number of arguments are passed
    if len(sys.argv) != 2:
        print("Usage: python sd2pdb.py <directory_path>")
        sys.exit(1)

    # Get the directory from the command-line arguments
    base_directory = sys.argv[1]

    # Check if the input is a valid directory
    if os.path.isdir(base_directory):
        convert_sd_to_pdb(base_directory)
    else:
        print(f"Error: {base_directory} is not a valid directory.")
        sys.exit(1)

if __name__ == "__main__":
    main()
