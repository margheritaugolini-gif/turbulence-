#!/usr/bin/env python3
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.io import loadmat, whosmat


COORDS = Path("schaefer1000_7N_coords_decolab.csv")
SC = Path("SC_schaefer1000_decolab.mat")
DATA = Path("timeseries_data_all.mat")


def cell_to_string(x):
    if isinstance(x, np.ndarray):
        if x.size == 0:
            return ""
        if x.dtype.kind in "US":
            return "".join(x.ravel().astype(str))
        return str(x.squeeze())
    return str(x)


def main():
    coords = pd.read_csv(COORDS)
    sc = loadmat(SC)["SC"]
    data_mat = loadmat(DATA, squeeze_me=False, struct_as_record=False)

    data = data_mat["Data"]
    filepaths = data_mat["filepaths"]
    subjects = [cell_to_string(x) for x in data_mat["subjects"][:, 0]]
    conditions = [cell_to_string(x) for x in data_mat["condition_order"][0]]

    print("Coordinate file")
    print(f"  rows: {len(coords)}")
    print(f"  columns: {list(coords.columns)}")
    print(f"  roi_label exactly 1..1000: {bool((coords['roi_label'].to_numpy() == np.arange(1, 1001)).all())}")
    print(f"  roi_name unique: {coords['roi_name'].nunique() == 1000}")
    print(f"  first parcel: {coords.loc[0, 'roi_name']}")
    print(f"  last parcel: {coords.loc[len(coords)-1, 'roi_name']}")

    print("\nSC file")
    print(f"  variables: {whosmat(SC)}")
    print(f"  shape: {sc.shape}")
    print(f"  symmetric: {bool(np.allclose(sc, sc.T))}")
    print(f"  diagonal all zero: {bool(np.allclose(np.diag(sc), 0))}")

    print("\nTime-series file")
    print(f"  Data shape: {data.shape}")
    print(f"  subjects: {len(subjects)}")
    print(f"  conditions: {conditions}")

    bad_cells = []
    bad_files = []
    used_files = []
    missing_cells = 0
    for i, subject in enumerate(subjects):
        for j, condition in enumerate(conditions):
            arr = data[i, j]
            if arr.size == 0:
                missing_cells += 1
            elif arr.ndim != 2 or arr.shape[1] != 1000:
                bad_cells.append((subject, condition, arr.shape))

            fp = cell_to_string(filepaths[i, j])
            if not fp:
                continue

            for path in fp.split("\n"):
                if not path:
                    continue
                used_files.append(path)
                name = Path(path).name
                if "schaefer1000-7networks_ts.txt" not in name:
                    bad_files.append((path, "wrong atlas filename"))
                match = re.search(r"_task-([^_]+)_", name)
                task = match.group(1) if match else None
                if task != condition:
                    bad_files.append((path, f"task {task} != condition {condition}"))

    print(f"  missing subject-condition cells: {missing_cells}")
    print(f"  selected source files: {len(used_files)}")
    print(f"  cells with wrong column count: {len(bad_cells)}")
    print(f"  source file/task issues: {len(bad_files)}")

    parts = coords["roi_name"].str.extract(r"7Networks_([LR]H)_([^_]+)", expand=True)
    hemi_network = (parts[0] + "_" + parts[1]).tolist()
    segments = []
    start = 1
    current = hemi_network[0]
    for idx, label in enumerate(hemi_network[1:], start=2):
        if label != current:
            segments.append((current, start, idx - 1, idx - start))
            current = label
            start = idx
    segments.append((current, start, len(hemi_network), len(hemi_network) - start + 1))

    print("\nCoordinate network order")
    for label, first, last, count in segments:
        print(f"  {label}: {first}-{last} ({count})")
    print(f"  each hemi-network is one contiguous block: {len(segments) == len(set(hemi_network))}")

    if bad_cells:
        print("\nCells with wrong shape")
        for item in bad_cells:
            print(f"  {item}")
    if bad_files:
        print("\nSource file issues")
        for item in bad_files:
            print(f"  {item}")


if __name__ == "__main__":
    main()
