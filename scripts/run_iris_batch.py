import os
import subprocess

base_dir = "/Users/jasonamburn/Desktop/naffinity_apo"
iris_scripts = "/Users/jasonamburn/Desktop/IRIS/scripts"
pymol_path = "/Applications/PyMOL.app/Contents/bin/pymol"
plugin_path = "/Applications/PyMOL.app/Contents/lib/python3.9/site-packages/pmg_tk/startup"

print("\nStarting FULL IRIS pipeline...\n")

for folder in sorted(os.listdir(base_dir)):

    folder_path = os.path.join(base_dir, folder)

    if not os.path.isdir(folder_path):
        continue

    print(f"\n==============================")
    print(f"Processing: {folder}")
    print(f"==============================\n")

    # Step 1 — Generate features
    print("Running get_features.py ...")

    result = subprocess.run(
        [
            "python",
            os.path.join(iris_scripts, "get_features.py"),
            folder_path,
            pymol_path,
            plugin_path
        ]
    )

    if result.returncode != 0:
        print(f"Feature generation FAILED for {folder}")
        continue

    print("Feature generation complete.\n")

    # Step 2 — Run IRIS ranking
    print("Running IRIS.py ...")

    result = subprocess.run(
        [
            "python",
            os.path.join(iris_scripts, "IRIS.py"),
            folder_path
        ]
    )

    if result.returncode != 0:
        print(f"IRIS ranking FAILED for {folder}")
        continue

    print(f"IRIS completed successfully for {folder}")

print("\nPipeline complete.\n")