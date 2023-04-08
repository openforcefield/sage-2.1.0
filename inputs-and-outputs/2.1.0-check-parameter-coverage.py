import logging
import os
import time
from multiprocessing import Pool, cpu_count, set_start_method

import click
from openff.toolkit.topology import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField

logging.getLogger("openff").setLevel(logging.ERROR)


def report_assigned_parameters(molecule, forcefield):
    params = []
    try:
        topology = Topology.from_molecules([molecule])
        # Run the molecule labeling
        molecule_force_list = forcefield.label_molecules(topology)

        # Print out a formatted description of the parameters
        # applied to this molecule

        for mol_idx, mol_forces in enumerate(molecule_force_list):
            for force_tag, force_dict in mol_forces.items():
                if (
                    force_tag == "Angles"
                    or force_tag == "Bonds"
                    or force_tag == "ProperTorsions"
                    or force_tag == "ImproperTorsions"
                ):
                    for atom_indices, parameter in force_dict.items():
                        params.append(parameter.id)
        return set(params)
    except Exception:
        print("Failed for", molecule.to_smiles(mapped=True))
        return params


def check_molecule(inputs):
    mol_idx = inputs[0]
    molecule = inputs[1]
    forcefield = inputs[2]

    # Prepare a title for this molecule
    if molecule.name == "":
        mol_name = f"molecule_{mol_idx + 1}"
    else:
        mol_name = molecule.name
    print("\n" * 3)
    print("=" * 60)
    print("=" * 60)
    print(f'Processing "{mol_name}" with smiles {molecule.to_smiles()}')
    print("=" * 60)
    print("=" * 60)
    # Analyze missing parameters
    time_i = time.time()
    assigned_params = report_assigned_parameters(molecule, forcefield)
    print(f"Molecule analysis took {time.time() - time_i} seconds")

    return (mol_name, assigned_params)


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
@click.option(
    "-ff",
    "--ff",
    "ff_to_modify",
    type=click.STRING,
    default="openff_unconstrained-2.0.0.offxml",
)
@click.option(
    "-rl",
    "--restrict_linear",
    "restrict_linear_params",
    type=click.BOOL,
    default=True,
)
def main(target_dir, ff_to_modify, restrict_linear_params):
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
    forcefield = ForceField(ff_to_modify, allow_cosmetic_attributes=True)
    job_args = [(idx, molecule, forcefield) for idx, molecule in enumerate(molecules)]
    result_list = p.map(check_molecule, job_args)
    results = dict(result_list)

    all_params = []
    for key, val in results.items():
        all_params.extend(val)

    all_pars = sorted(list(set(all_params)))
    print(all_pars)
    # remove linear angles for fitting
    # a16, a17, a27, a35
    # remove linear torsions
    # t165, t166, t167 directly
    if restrict_linear_params:
        all_pars = [item for item in all_pars if item not in ["t165", "t166", "t167"]]

    forcefield = ForceField(ff_to_modify)
    for ptype in ["Angles", "Bonds", "ProperTorsions", "ImproperTorsions"]:
        for i in range(len(forcefield.get_parameter_handler(ptype).parameters)):
            item = forcefield.get_parameter_handler(ptype).parameters[i]
            dd = item.to_dict()
            if dd["id"] in all_pars:
                dd["parameterize"] = []
                for key, value in dd.items():
                    if (
                        key.startswith("k")
                        or key.startswith("angle")
                        or key.startswith("length")
                    ):
                        if key.startswith("angle") and dd["id"] not in [
                            "a16",
                            "a17",
                            "a27",
                            "a35",
                        ]:
                            dd["parameterize"].append(key)
                        elif not key.startswith("angle"):
                            dd["parameterize"].append(key)
                dd["parameterize"] = ",".join(dd["parameterize"])
                forcefield.get_parameter_handler(ptype).parameters[
                    i
                ].add_cosmetic_attribute("parameterize", dd["parameterize"])

    forcefield.to_file("cosmetic_attributes_added_to_forcefield.offxml")
    print(
        "FORCEFIELD with cosmetic attributes \
        written to cosmetic_attributes_added_to_forcefield.offxml"
    )

    print(
        f"Processing {len(molecules)} molecules took"
        f"{time.time() - start_time} seconds"
    )


if __name__ == "__main__":
    main()
