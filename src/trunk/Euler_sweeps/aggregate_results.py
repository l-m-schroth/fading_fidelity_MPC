"""
This file can be used to aggregate the results of a Euler sweep
"""

#!/usr/bin/env python3
import os
import pickle
import argparse
import sys
from utils_shared import get_dir

parent_dir = os.path.dirname(os.path.dirname(__file__))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

def aggregate_results(results_dir, prefix, aggregated_filename):
    results = []
    if not os.path.exists(results_dir):
        print(f"Results directory '{results_dir}' does not exist.")
        return

    for filename in os.listdir(results_dir):
        if filename.startswith(prefix) and filename.endswith(".pkl"):
            full_path = os.path.join(results_dir, filename)
            with open(full_path, "rb") as f:
                try:
                    results.append(pickle.load(f))
                except Exception as e:
                    print(f"Error loading {full_path}: {e}")

    aggregated_dir = get_dir("data/trunk")
    os.makedirs(aggregated_dir, exist_ok=True)
    aggregated_path = os.path.join(aggregated_dir, aggregated_filename)
    with open(aggregated_path, "wb") as f:
        pickle.dump(results, f)

    print(f"Aggregated {len(results)} files into {aggregated_path}")

def main():
    default_euler = os.path.join(get_dir("src/trunk/Euler_sweeps"), "results_18_03")
    parser = argparse.ArgumentParser(description="Aggregate result files with a specific prefix.")
    parser.add_argument("--results_dir", type=str, default=default_euler)
    parser.add_argument("--prefix", type=str, default="results_16_link_ellipses_")
    parser.add_argument("--aggregated_filename", type=str, default="results_16_link_oval_aggregated.pkl")

    args = parser.parse_args()
    aggregate_results(args.results_dir, args.prefix, args.aggregated_filename)

if __name__ == "__main__":
    main()

