input_file = "pau_333s_351g_partitions.ufboot"

# Read the file and extract each line (tree)
with open(input_file, "r") as file:
    trees = file.readlines()

# Loop through each tree and save it as a separate file
for i, tree in enumerate(trees, start=1):
    # Create a unique file name for each tree
    output_file = f"tree_{i}.tre"
    
    # Write the tree to a new file
    with open(output_file, "w") as out_file:
        out_file.write(tree.strip())  # .strip() removes any leading/trailing whitespace
    
    print(f"Generated {output_file}")
