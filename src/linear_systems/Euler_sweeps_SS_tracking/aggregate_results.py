#!/usr/bin/env python3
import os
import glob
import pickle

# === USER-DEFINED SETTINGS ===
RESULTS_FOLDER = "results"
RESULTS_PREFIX = "n50_nu1_nonoise_lh2.5"   # e.g. your base experiment name
OUTPUT_FILE = f"{RESULTS_FOLDER}/{RESULTS_PREFIX}_aggregated_total.pkl"

def aggregate_results():
    all_data = []
    
    # matches all pkl files that contain the prefix. 
    # you might want something like:
    search_pattern = os.path.join(
        RESULTS_FOLDER,
        f"{RESULTS_PREFIX}_sched*_*_results_*.pkl"
    )
    pkl_files = glob.glob(search_pattern)

    if not pkl_files:
        print(f"No matching files found for prefix: {RESULTS_PREFIX}")
        return

    print(f"Found {len(pkl_files)} matching result files for prefix '{RESULTS_PREFIX}'.")

    for fpath in pkl_files:
        with open(fpath, "rb") as f:
            data = pickle.load(f)
        all_data.append(data)

    with open(OUTPUT_FILE, "wb") as f:
        pickle.dump(all_data, f)

    print(f"Aggregated {len(all_data)} results into: {OUTPUT_FILE}")

if __name__ == "__main__":
    aggregate_results()
