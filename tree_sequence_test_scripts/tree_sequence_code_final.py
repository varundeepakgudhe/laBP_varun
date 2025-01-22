import tskit
import pandas as pd
from IPython.display import SVG

# Load the CSV file
data = pd.read_csv("data.csv")

# Initialize TableCollection with a large enough sequence length
tables = tskit.TableCollection(sequence_length=1000.011)

# Add a single population (population 0)
tables.populations.add_row(metadata=b"")

# Collect all unique node IDs from `nodeNumber` and `descendantNodeNumbers`
all_nodes = set(data["nodeNumber"])  # Start with parent nodes

for descendants in data["descendantNodeNumbers"]:
    if pd.notna(descendants):
        # Replace unexpected separators with a comma, then split the values
        cleaned_descendants = str(descendants).replace(";", ",")
        descendant_list = [int(x) for x in cleaned_descendants.split(",") if x.strip().isdigit()]
        all_nodes.update(descendant_list)

# Populate the individuals table (only sample nodes with time = 0)
for index, row in data.iterrows():
    time = float(row["time"])
    if time == 0:  # Add to individuals table only if the node is a sample
        tables.individuals.add_row(
            flags=0,  # Flags set to 0 for all individuals
            location=[],  # Empty location
            parents=[],  # Empty parents
            metadata=b""  # Empty metadata
        )

# Populate the nodes table
for node in sorted(all_nodes):  # Sorting ensures consistent ordering
    node_data = data[data["nodeNumber"] == node]
    if not node_data.empty:
        time = round(float(node_data.iloc[0]["time"]), 3)  # Ensure time has 3 decimal points
    else:
        time = 0.0  # Default to 0 if the node is not found
    tables.nodes.add_row(
        flags=1 if time == 0 else 0,  # Sample nodes have flags=1
        population=0,  # Population set to 0 for all nodes
        individual=node if time == 0 else tskit.NULL,  # Link sample nodes to individuals
        time=time,
        metadata=b"",  # Empty metadata
    )

# Populate the edges table
for index, row in data.iterrows():
    parent_node = row["nodeNumber"]
    descendants = row["descendantNodeNumbers"]

    if pd.notna(descendants):
        # Replace unexpected separators with a comma, then split the values
        cleaned_descendants = str(descendants).replace(";", ",")
        descendant_list = [int(x) for x in cleaned_descendants.split(",") if x.strip().isdigit()]
        for descendant in descendant_list:
            tables.edges.add_row(
                left=1000.01,
                right=1000.011,
                parent=parent_node,
                child=descendant
            )
    # Skip adding edges for rows without descendants

# Sort the tables (this fixes the TSK_ERR_EDGES_NOT_SORTED_CHILD error)
tables.sort()

# Generate the tree sequence
ts = tables.tree_sequence()

# Save the SVG to a file
svg_content = ts.draw_svg(
    size=(800, 800),  # Set the size of the SVG (width, height)
    y_axis=True  # Include the Y-axis without label rotation
)

# Save the SVG to a file
with open("tree_enlarged.svg", "w") as svg_file:
    svg_file.write(svg_content)

print("SVG tree saved as 'tree_enlarged.svg'. Open this file in a browser to view it.")