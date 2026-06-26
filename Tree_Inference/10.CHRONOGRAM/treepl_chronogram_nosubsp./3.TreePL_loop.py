import os
import shutil
import subprocess

# Define the paths
tree_dir = "prunedTrees"
output_dir = "chronograms"
config_template = "config_template.txt"

# Ensure the "chronograms" directory exists
os.makedirs(output_dir, exist_ok=True)

# Get the list of tree files in the prunedTrees directory
tree_files = [f for f in os.listdir(tree_dir) if f.endswith(".tre")]

# Check if there are tree files
if not tree_files:
    print("No tree files found in the 'prunedTrees' directory.")
    exit()

# Loop through each tree file
for i, tree in enumerate(tree_files, start=1):
    # Create a unique config file name
    config_file = f"config_{i}.txt"
    
    # Copy the template config file to create a new config file
    if os.path.exists(config_template):
        shutil.copy(config_template, config_file)
    else:
        print("Error: config_template.txt not found.")
        continue
    
    # Read the content of the new config file
    with open(config_file, "r") as file:
        config_content = file.read()

    # Replace placeholders with actual file names
    config_content = config_content.replace("treefile =", f"treefile = {os.path.join(tree_dir, tree)}")
    config_content = config_content.replace("outfile =", f"outfile = {os.path.join(output_dir, f'output_{i}.tre')}")

    # Write the modified content back to the config file
    with open(config_file, "w") as file:
        file.write(config_content)
    
    # Confirm that the file was generated
    if os.path.exists(config_file):
        print(f"Generated {config_file} for tree {tree}")
    else:
        print(f"Failed to generate {config_file} for tree {tree}")
        continue
    
    # Run treePL with the unique config file
    print(f"Running treePL for {config_file}...")
    result = subprocess.run(["treePL", config_file], capture_output=True, text=True)
    
    # Print the output and errors (if any) from the treePL command
    if result.returncode == 0:
        print(f"treePL ran successfully for {config_file}")
    else:
        print(f"Error running treePL for {config_file}: {result.stderr}")
    
    # Optionally, clean up the config file if not needed
    os.remove(config_file)
