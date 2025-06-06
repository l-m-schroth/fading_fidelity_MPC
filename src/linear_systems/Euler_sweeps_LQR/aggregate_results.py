#!/usr/bin/env python3
import os
import glob
import pickle

# === USER-DEFINED SETTINGS ===
results_folder = "results"  # Folder where all the results are stored
results_prefix = "n50_nu1"  # Change this to match the prefix you want to aggregate
output_file = f"{results_folder}/{results_prefix}_aggregated.pkl"  # Output file

# === AGGREGATION LOGIC ===
def aggregate_results():
    all_data = []
    
    # Find all .pkl files that match the given prefix
    search_pattern = os.path.join(results_folder, f"{results_prefix}_lh*_results_*.pkl")
    pkl_files = glob.glob(search_pattern)

    if not pkl_files:
        print(f"No matching files found for prefix: {results_prefix}")
        return

    print(f"Found {len(pkl_files)} matching result files for prefix '{results_prefix}'.")

    # Load and combine all results
    for fpath in pkl_files:
        with open(fpath, "rb") as f:
            data = pickle.load(f)
        all_data.extend(data)

    # Dump the aggregated list directly into a pickle file
    with open(output_file, "wb") as f:
        pickle.dump(all_data, f)

    print(f"Aggregated {len(all_data)} results into: {output_file}")

# Run the aggregation
aggregate_results()
