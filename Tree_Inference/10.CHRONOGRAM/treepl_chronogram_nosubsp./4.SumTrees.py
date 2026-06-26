#!/usr/bin/env python3
import subprocess
import os

# Get the absolute path of the directory where THIS script is running
current_dir = os.path.dirname(os.path.abspath(__file__))

# Define paths relative to the local directory
# Based on your error log, sumtrees.py is inside DendroPy/src/dendropy/application/
sumtrees_script = os.path.join(current_dir, "DendroPy", "src", "dendropy", "application", "sumtrees.py")
output_file = os.path.join(current_dir, "SumTrees_chronogram.nex")
input_tree = os.path.join(current_dir, "chronograms", "combined_chronograms.tre")

# CRITICAL FIX: Tell Python where the 'dendropy' package live source is
dendropy_src_dir = os.path.join(current_dir, "DendroPy", "src")
custom_env = os.environ.copy()
custom_env["PYTHONPATH"] = dendropy_src_dir + os.pathsep + custom_env.get("PYTHONPATH", "")

# Reconstruct your exact command as a list of arguments
command = [
    "python3", sumtrees_script,
    "-i", "newick",
    "-F", "nexus",
    "-o", output_file,
    "--summarize-node-ages",
    "--root-target-at-outgroup=Talisia_nervosa",
    input_tree
]

print(f"Running SumTrees analysis in: {current_dir}...")

# Run the command with the custom environment paths included
result = subprocess.run(command, capture_output=True, text=True, env=custom_env)

# Check if it worked
if result.returncode == 0:
    print("Success! Your summarized chronogram is ready in the local directory.")
else:
    print("Something went wrong. Error output:")
    print(result.stderr)