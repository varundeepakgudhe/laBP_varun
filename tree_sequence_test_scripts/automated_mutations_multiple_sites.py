import re
import tskit
import msprime
import pandas as pd

# 1. Create a TableCollection with a fixed sequence length.
sequence_length = 1000
tables = tskit.TableCollection(sequence_length=sequence_length)

# 2. Create populations (assuming context_pop is 0 or 1).
pop0 = tables.populations.add_row(metadata=b'pop0')
pop1 = tables.populations.add_row(metadata=b'pop1')

# 3. List of CSV files (each representing a tree at a given site position).
csv_files = [
    "1.genetree_site_400.csv",
    "1.genetree_site_600.csv",
    "1.genetree_site_800.csv"
    # Add more filenames here as needed.
]

# 4. Extract the site positions from filenames.
site_positions = {}
for filename in csv_files:
    # Expecting filenames like "1.genetree_site_400.csv"
    match = re.search(r"_site_(\d+)", filename)
    if match:
        site = int(match.group(1))
        site_positions[filename] = site
    else:
        raise ValueError(f"Could not extract site position from filename {filename}")

# 5. Compute the sorted site positions and then determine interval boundaries.
sorted_sites = sorted(site_positions.values())
boundaries = [0]  # L0 is always 0.
# For each adjacent pair of site positions, use the midpoint as a boundary.
for i in range(len(sorted_sites) - 1):
    boundaries.append((sorted_sites[i] + sorted_sites[i + 1]) / 2)
boundaries.append(sequence_length)  # Final boundary is the sequence length.

# 6. Map each site to its corresponding interval (left, right).
site_intervals = {}
for site in sorted_sites:
    idx = sorted_sites.index(site)
    left = boundaries[idx]
    right = boundaries[idx + 1]
    site_intervals[site] = (left, right)

# 7. Prepare containers to combine data across CSVs.
nodes_dict = {}         # Key: nodeNumber, Value: dict with node info.
individual_ids = {}     # Key: nodeNumber (for samples, time==0), Value: individual id.
edges_dict = {}         # Key: (parent, child), Value: list of intervals (tuples).

# Helper function to merge overlapping or contiguous intervals.
def merge_intervals(intervals):
    if not intervals:
        return []
    # Sort intervals by left coordinate.
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [intervals[0]]
    for current in intervals[1:]:
        prev = merged[-1]
        # Merge if the current interval touches or overlaps the previous one.
        if current[0] <= prev[1]:
            merged[-1] = (prev[0], max(prev[1], current[1]))
        else:
            merged.append(current)
    return merged

# 8. Process each CSV file.
for filename in csv_files:
    # Determine the site position and its effective interval.
    site = site_positions[filename]
    interval = site_intervals[site]  # (left, right)
    
    # Read CSV data.
    df = pd.read_csv(filename, skip_blank_lines=True)

    df = df.dropna(subset=["nodeNumber"])
    
    # Ensure columns are of the correct type.
    df["nodeNumber"] = df["nodeNumber"].astype(int)
    df["time"] = df["time"].astype(float)
    df["context_pop"] = df["context_pop"].astype(int)
    df["descendantNodeNumbers"] = df["descendantNodeNumbers"].fillna("").astype(str)
    
    # Convert DataFrame to a list of dictionary records.
    records = df[["nodeNumber", "time", "context_pop", "descendantNodeNumbers"]].to_dict(orient="records")
    
    # Process nodes and individuals.
    for row in records:
        node_id = row["nodeNumber"]
        time = row["time"]
        pop = row["context_pop"]
        if node_id not in nodes_dict:
            nodes_dict[node_id] = {"time": time, "pop": pop}
            # Create an individual if the node is a sample (time == 0).
            if time == 0:
                ind_id = tables.individuals.add_row()
                individual_ids[node_id] = ind_id
    
    # Process edges.
    for row in records:
        desc = row["descendantNodeNumbers"].strip()
        if desc != "":
            parent = row["nodeNumber"]
            # Convert descendant string into a list of integers.
            descendants = []
            for x in desc.split(";"):
                try:
                    descendants.append(int(x))
                except ValueError:
                    continue
            for child in descendants:
                key = (parent, child)
                if key not in edges_dict:
                    edges_dict[key] = []
                edges_dict[key].append(interval)

# 9. Add nodes to the table and build a mapping from CSV node number to table node ID.
csv_to_table = {}
sample_flag = tskit.NODE_IS_SAMPLE
for node_id, info in nodes_dict.items():
    time = info["time"]
    pop = info["pop"]
    if time == 0:
        flags = sample_flag
        individual = individual_ids[node_id]
    else:
        flags = 0
        individual = tskit.NULL
    # Add the node and capture the table's node id.
    table_node_id = tables.nodes.add_row(time=time, flags=flags, population=pop, individual=individual)
    csv_to_table[node_id] = table_node_id

# 10. Merge intervals for each edge and add them to the edges table.
for (parent_csv, child_csv), intervals in edges_dict.items():
    # Only add the edge if both parent and child are in the mapping.
    if parent_csv in csv_to_table and child_csv in csv_to_table:
        parent_table = csv_to_table[parent_csv]
        child_table = csv_to_table[child_csv]
        merged_intervals = merge_intervals(intervals)
        for (left, right) in merged_intervals:
            tables.edges.add_row(left=left, right=right, parent=parent_table, child=child_table)
    else:
        print(f"Skipping edge ({parent_csv} -> {child_csv}) because one or both nodes are missing in the mapping.")

# 11. Finalize and build the tree sequence.
tables.sort()
ts_manual = tables.tree_sequence()

# 12. Apply mutations.
ts_with_mut = msprime.sim_mutations(
    ts_manual,
    rate=1e-7,  # Adjust mutation rate as needed.
    model=msprime.InfiniteSites(),
    random_seed=42
)

# 13. Print summary information.
print("Number of nodes:", ts_with_mut.num_nodes)
print("Number of edges:", ts_with_mut.num_edges)
print("Number of sites:", ts_with_mut.num_sites)
print("Number of mutations:", ts_with_mut.num_mutations)

# 14. Optionally, save an SVG of the final tree.
svg_str = ts_with_mut.draw_svg(size=(800, 400))
with open("manual_tree_with_mutations.svg", "w") as f:
    f.write(svg_str)
print("Saved manual_tree_with_mutations.svg. Open it in a browser to view.")

#15. (Optional) Print required tables
import pandas as pd

# Ensure Pandas prints all rows and columns
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

tables = ts_with_mut.tables

# List of table names (modify as needed)
table_names = [
    "nodes",
    "edges",
    "sites",
    "mutations",
    "individuals",
    "populations",
    "migrations",
    "provenances"
]

for name in table_names:
    table = getattr(tables, name, None)
    if table is not None:
        print(f"----- {name.capitalize()} Table -----")
        try:
            # Try using the built-in method if available
            df = table.to_dataframe()
        except AttributeError:
            # If not, build the DataFrame row by row
            rows = [table[i].asdict() for i in range(table.num_rows)]
            df = pd.DataFrame(rows)
        # Print the full table as a string
        print(df.to_string())
        # Save the DataFrame to a CSV file
        csv_filename = f"{name}_table.csv"
        df.to_csv(csv_filename, index=False)
        print(f"Saved {name} table to {csv_filename}\n")