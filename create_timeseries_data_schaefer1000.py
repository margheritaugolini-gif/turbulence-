#!/usr/bin/env python3
"""
Build the MATLAB input file used by the Schaefer-1000 analysis pipeline.

Purpose:
Collect subject-level time-series text files from the `time_series/` folder,
stack multiple runs from the same condition when present, and save the result
as `timeseries_data_all.mat` for the MATLAB preprocessing scripts.

Expected input layout:
time_series/
  sub-*/
    memory/*_task-memory_*schaefer1000-7networks_ts.txt
    counting/*_task-counting_*schaefer1000-7networks_ts.txt
    jhana/*_task-j1_*schaefer1000-7networks_ts.txt
    ...
    jhana/*_task-j8_*schaefer1000-7networks_ts.txt

Each text file must be a numeric matrix with shape:
  time points x 1000 Schaefer parcels

Main parameters to change:
ROOT:
    Folder containing the `sub-*` directories.
OUTPUT:
    Name of the generated MATLAB file.
BACKUP:
    One-time backup name used if OUTPUT already exists.
CONDITIONS:
    Condition/task order saved into the MATLAB file.
ATLAS_SUFFIX:
    Filename suffix used to select the intended atlas/parcellation files.

Outputs:
timeseries_data_all.mat containing:
  Data, subjects, condition_order, filepaths, run_counts, Tmap, source_root, atlas

Run:
From this project folder:
  python3 create_timeseries_data_schaefer1000.py
"""

from pathlib import Path
import re

import numpy as np
from scipy.io import savemat


ROOT = Path("time_series")
OUTPUT = Path("timeseries_data_all.mat")
BACKUP = Path("timeseries_data_all_schaefer200_backup.mat")
CONDITIONS = ["memory", "counting"] + [f"j{i}" for i in range(1, 9)]
CONDITION_LABELS = CONDITIONS
ATLAS_SUFFIX = "schaefer1000-7networks_ts.txt"


def task_name(path: Path) -> str:
    match = re.search(r"_task-([^_]+)_", path.name)
    return match.group(1) if match else ""


def load_condition_files(subject_dir: Path, condition: str):
    folder = "jhana" if condition.startswith("j") else condition
    condition_dir = subject_dir / folder
    if not condition_dir.exists():
        return [], None

    files = [
        path
        for path in sorted(condition_dir.glob(f"*{ATLAS_SUFFIX}"))
        if task_name(path) == condition
    ]

    if not files:
        return [], None

    runs = []
    for path in files:
        arr = np.loadtxt(path)
        if arr.ndim != 2 or arr.shape[1] != 1000:
            raise ValueError(f"{path} has shape {arr.shape}; expected time x 1000")
        runs.append(arr)

    return files, np.vstack(runs)


def make_cell_array(shape):
    cell = np.empty(shape, dtype=object)
    cell[:] = ""
    return cell


def main():
    if not ROOT.exists():
        raise FileNotFoundError(f"Missing folder: {ROOT}")

    subject_dirs = sorted(path for path in ROOT.glob("sub-*") if path.is_dir())
    if not subject_dirs:
        raise RuntimeError(f"No subject folders found under {ROOT}")

    n_sub = len(subject_dirs)
    n_cond = len(CONDITIONS)

    data = make_cell_array((n_sub, n_cond))
    filepaths = make_cell_array((n_sub, n_cond))
    run_counts = np.zeros((n_sub, n_cond), dtype=np.uint16)
    tmap = np.zeros((n_sub, n_cond), dtype=np.uint16)

    for i, subject_dir in enumerate(subject_dirs):
        for j, condition in enumerate(CONDITIONS):
            files, arr = load_condition_files(subject_dir, condition)
            if arr is None:
                data[i, j] = np.empty((0, 1000), dtype=np.float64)
                filepaths[i, j] = ""
                continue

            data[i, j] = arr.astype(np.float64, copy=False)
            filepaths[i, j] = "\n".join(str(path) for path in files)
            run_counts[i, j] = len(files)
            tmap[i, j] = arr.shape[0]

    subjects = np.array([[path.name] for path in subject_dirs], dtype=object)
    condition_order = np.array([CONDITION_LABELS], dtype=object)

    if OUTPUT.exists() and not BACKUP.exists():
        BACKUP.write_bytes(OUTPUT.read_bytes())

    savemat(
        OUTPUT,
        {
            "Data": data,
            "subjects": subjects,
            "condition_order": condition_order,
            "filepaths": filepaths,
            "run_counts": run_counts,
            "Tmap": tmap,
            "source_root": str(ROOT),
            "atlas": "schaefer1000-7networks",
        },
        do_compression=True,
    )

    complete = np.all(run_counts > 0, axis=1)
    print(f"Saved {OUTPUT}")
    if BACKUP.exists():
        print(f"Backup present: {BACKUP}")
    print(f"Subjects: {n_sub}")
    print(f"Conditions: {', '.join(CONDITION_LABELS)}")
    print(f"Complete subjects across all 10 selected conditions: {int(complete.sum())}")

    missing = []
    multiple = []
    for i, subject in enumerate(subjects[:, 0]):
        miss = [CONDITION_LABELS[j] for j in range(n_cond) if run_counts[i, j] == 0]
        mult = [
            f"{CONDITION_LABELS[j]}({run_counts[i, j]})"
            for j in range(n_cond)
            if run_counts[i, j] > 1
        ]
        if miss:
            missing.append(f"{subject}: {', '.join(miss)}")
        if mult:
            multiple.append(f"{subject}: {', '.join(mult)}")

    if missing:
        print("\nMissing selected conditions:")
        print("\n".join(missing))
    if multiple:
        print("\nConcatenated multiple runs:")
        print("\n".join(multiple))


if __name__ == "__main__":
    main()
