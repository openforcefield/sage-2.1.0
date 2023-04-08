import logging
import os
import time
from multiprocessing import Pool, cpu_count, set_start_method

import click
from openff.toolkit.topology import Molecule, Topology

logging.getLogger("openff").setLevel(logging.ERROR)


def report_charged_status(molecule, mol_dir):
    try:
        molecule.assign_partial_charges(partial_charge_method="am1bccelf10")
        return True
    except Exception:
        print("Failed for", molecule.to_smiles(mapped=True))
        print(mol_dir)
        return (mol_dir, molecule.to_smiles(mapped=True))


def check_molecule(inputs):
    mol_idx = inputs[0]
    molecule = inputs[1]
    mol_dir = inputs[2]

    # Prepare a title for this molecule
    if molecule.name == "":
        mol_name = f"molecule_{mol_idx + 1}"
    else:
        mol_name = molecule.name
    # Analyze missing parameters
    time_i = time.time()
    charged_or_not = report_charged_status(molecule, mol_dir)

    return (mol_name, charged_or_not)


set_start_method("fork")
num_threads = max(1, int(cpu_count() * 1.0))


@click.command()
@click.option(
    "-tdir",
    "--target_dir",
    "target_dir",
    type=click.STRING,
    default="./targets/",
)
def main(target_dir):
    subdirs = []
    for root, dirs, files in os.walk(target_dir):
        for file in files:
            if file.endswith(".sdf") or file.endswith(".mol2"):
                subdirs.append(os.path.join(root, file))

    molecules = [
        Molecule.from_file(subdir, file_format="sdf", allow_undefined_stereo=True)
        if subdir.endswith("sdf")
        else Molecule.from_file(subdir, file_format="mol2", allow_undefined_stereo=True)
        for subdir in subdirs
    ]

    start_time = time.time()
    p = Pool(num_threads)
    job_args = [(idx, molecule, subdirs[idx]) for idx, molecule in enumerate(molecules)]
    result_list = p.map(check_molecule, job_args)
    results = dict(result_list)

    print(
        f"Processing {len(molecules)} molecules took"
        f"{time.time() - start_time} seconds"
    )

    with open("failed_mol_dirs.txt", "w") as f:
        for item in results.values():
            if item != True:
                f.write(f"{item[0]}\n")
    with open("failed_mol_smiles.txt", "w") as f:
        for item in results.values():
            if item != True:
                f.write(f"{item[1]}\n")
    with open("failed_mol_opt_ids.txt", "w") as f:
        for item in results.values():
            if item != True:
                f.write(f"{os.path.basename(item[0][:-6])}\n")


if __name__ == "__main__":
    main()
