"""
    This script reads a (output) text file from CRADLE, writing them to a new output file. 
    The output file will have the same name as the input file, but with "_short" appended 
    before the file extension.

    Columns in the input file:
        1st column: id of event (id)
        2nd column: half life of event (hl)
        3rd column: name of particle (name)
        4th column: i dont know (idk)
        5th column: kinetic energy of particle (p)
        6th column: total energy of particle (E=m+p) (E)
        7th column: momentum x of particle (px)
        8th column: momentum y of particle (py)
        9th column: momentum z of particle (pz)

    Usage: 
        python3 short.py input.txt
        python3 short.py input.txt particle p px py pz
"""

import sys
import os

# -----------------------------
# Column mapping
# -----------------------------
col_map = {
    "id": 0,
    "hl": 1,
    "particle": 2,
    "idk": 3,
    "p": 4,
    "E": 5,
    "px": 6,
    "py": 7,
    "pz": 8
}

# -----------------------------
# Input arguments
# -----------------------------
if len(sys.argv) < 2:
    print("Usage: python3 short.py input.txt [columns...]")
    sys.exit(1)

input_file = sys.argv[1]

# columns requested by user
requested_cols = sys.argv[2:]

# default behavior
if not requested_cols:
    requested_cols = ["particle", "p"]

# output file name
base_name = os.path.splitext(input_file)[0]
output_file = base_name + "_short.txt"

with open(input_file, "r") as f_in, open(output_file, "w") as f_out:
    for line in f_in:
        if not line.strip():
            continue

        cols = line.split()

        # safety check
        if len(cols) < 9:
            continue

        try:
            values = [cols[col_map[c]] for c in requested_cols]
        except KeyError as e:
            print(f"Unknown column: {e}")
            print(f"Available columns: {list(col_map.keys())}")
            sys.exit(1)

        f_out.write("\t".join(values) + "\n")

print(f"Written to: {output_file}")