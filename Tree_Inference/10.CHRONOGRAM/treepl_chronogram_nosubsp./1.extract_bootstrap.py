import random

input_file = "pau_333s_351g_partitions.ufboot"
output_file = "pau_333s_351g_partitions.ufboot_100.tre"

# Read the file and extract each line (tree)
with open(input_file, "r") as file:
    # Read lines and keep only non-empty ones
    trees = [line.strip() for line in file if line.strip()]

# Check if there are enough trees to sample from
if len(trees) < 100:
    print(f"Warning: Only found {len(trees)} trees. Saving all of them instead.")
    sampled_trees = trees
else:
    # Set seed for reproducibility (optional, remove if you want it truly random every time)
    random.seed(42)
    # Pick 100 random trees from the list
    sampled_trees = random.sample(trees, 100)

# Write all 100 random trees into the single new file
with open(output_file, "w") as out_file:
    for tree in sampled_trees:
        out_file.write(tree + "\n")

print(f"Successfully generated {output_file} with {len(sampled_trees)} random trees!")
