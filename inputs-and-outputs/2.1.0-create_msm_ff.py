# Modified code from https://github.com/jthorton/MSM_QCArchive/blob/master/Mod_sem.ipynb
# we need a function to calculate the modified seminario parameters
# type the molecule using openff and store the parameters under the smirks pattern
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.qcsubmit.results import BasicResultCollection, OptimizationResultCollection
from openff.qcsubmit.results.filters import LowestEnergyFilter
from qubekit.molecules import Ligand
from qubekit.bonded.mod_seminario import ModSeminario
from collections import defaultdict
import json
from simtk import unit
import numpy as np
import click
import os

# we need to remove the openeye wrapper to avoid stereochemistry issues
# as openeye gives nitrogen stereo flags and rdkit does not
GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper())
GLOBAL_TOOLKIT_REGISTRY
mod_sem = ModSeminario()


def calculate_parameters(qc_record, off_molecule, ff):
    """
    Calculate the modified seminario parameters for the given input molecule and store them by OFF SMIRKS.
    """
    # create the qube molecule, this should be in the same order as the off_mol
    qube_mol = Ligand.from_rdkit(off_molecule.to_rdkit())
    qube_mol.hessian = qc_record.return_result
    # calculate the modified seminario parameters and store in the molecule
    qube_mol = mod_sem.run(qube_mol)
    # label the openff molecule
    labels = ff.label_molecules(off_molecule.to_topology())[0]
    # loop over all bonds and angles and collect the results in nm/ kj/mol / radians(openMM units)
    bond_eq, bond_k, angle_eq, angle_k = (
        defaultdict(list),
        defaultdict(list),
        defaultdict(list),
        defaultdict(list),
    )
    for bond, parameter in labels["Bonds"].items():
        # bond is a tuple of the atom index the parameter is applied to
        qube_param = qube_mol.BondForce[bond]
        bond_eq[parameter.smirks].append(qube_param.length)
        bond_k[parameter.smirks].append(qube_param.k)
    for angle, parameter in labels["Angles"].items():
        qube_param = qube_mol.AngleForce[angle]
        angle_eq[parameter.smirks].append(qube_param.angle)
        angle_k[parameter.smirks].append(qube_param.k)

    return bond_eq, bond_k, angle_eq, angle_k


def update_parameters(master_params, new_params):
    """Update the new parameters into the master parameters."""
    for smirks, parameters in new_params.items():
        master_params[smirks].extend(parameters)


@click.command()
@click.option(
    "--initial_ff",
    "initial_ff",
    type=click.STRING,
    required=True,
)
@click.option(
    "--output_ff",
    "output_ff",
    type=click.STRING,
    required=True,
)
@click.option(
    "--opt_json",
    "opt_json",
    type=click.STRING,
    required=True,
)
@click.option(
    "--output_dir",
    "output_dir",
    type=click.STRING,
    required=True,
    default="./",
)
def main(initial_ff, output_ff, opt_json, output_dir):
    default_filters = [LowestEnergyFilter()]
    optimization_set = OptimizationResultCollection.parse_file(opt_json)
    optimization_set = optimization_set.filter(*default_filters)
    # now we want only those entries we calculated hessians for
    hessian_set = optimization_set.to_basic_result_collection(driver="hessian")
    # save the hessian dataset
    with open(
        output_dir + "hessian-set-used-in-creating-msm-starting-point.json", "w"
    ) as file:
        file.write(hessian_set.json())
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # when restarting pull from here
    # hessian_set = BasicResultCollection.parse_file(
    #    "hessian_set.json")
    # we should have one result per molecule now
    print(hessian_set.n_molecules)
    print(hessian_set.n_results)

    # set up the class
    ff = ForceField(initial_ff, allow_cosmetic_attributes=True)

    # pull down each record and molecule from the dataset
    # this gives a list of tuples (record, molecule)
    records_and_molecules = hessian_set.to_records()

    # loop over all of the records and calculate the bond and angle terms
    all_bond_eq, all_bond_k, all_angle_eq, all_angle_k = (
        defaultdict(list),
        defaultdict(list),
        defaultdict(list),
        defaultdict(list),
    )
    for record, molecule in records_and_molecules:
        try:
            bonds_eq, bonds_k, angles_eq, angles_k = calculate_parameters(
                record, molecule, ff
            )
        except KeyError:
            print(record)
            print(molecule)
            raise KeyError
        update_parameters(all_bond_eq, bonds_eq)
        update_parameters(all_bond_k, bonds_k)
        update_parameters(all_angle_eq, angles_eq)
        update_parameters(all_angle_k, angles_k)

    # make one dict with all the data and store in json
    all_parameters = {
        "bonds_eq": all_bond_eq,
        "bonds_k": all_bond_k,
        "angles_eq": all_angle_eq,
        "angles_k": all_angle_k,
    }
    with open(output_dir + "/seminario_parameters.json", "w") as output:
        output.write(json.dumps(all_parameters))

    # now lets edit the offxml and save to file
    edit_ff = ForceField(initial_ff, allow_cosmetic_attributes=True)
    # first do bonds
    # lengths in nanometer, k in kj/mol nm**2
    # converting to angstroms and kcal/mol/ang^2
    bond_handler = edit_ff.get_parameter_handler("Bonds")
    for smirks in all_bond_eq.keys():
        off_bond_param = bond_handler.parameters[smirks]
        bond_length = unit.Quantity(np.mean(all_bond_eq[smirks]), unit=unit.nanometer)
        off_bond_param.length = bond_length.in_units_of(unit.angstrom)
        bond_k = unit.Quantity(
            np.mean(all_bond_k[smirks]),
            unit=unit.kilojoule_per_mole * unit.nanometer**-2,
        )
        off_bond_param.k = bond_k.in_units_of(
            unit.kilocalorie_per_mole * unit.angstrom**-2
        )

    # now angles
    # radians and kj/mol radian**2
    angle_handler = edit_ff.get_parameter_handler("Angles")
    for smirks in all_angle_eq.keys():
        off_angle_param = angle_handler.parameters[smirks]
        angle = unit.Quantity(np.mean(all_angle_eq[smirks]), unit=unit.radian)
        off_angle_param.angle = angle.in_units_of(unit.degree)
        angle_k = unit.Quantity(
            np.mean(all_angle_k[smirks]),
            unit=unit.kilojoules_per_mole * unit.radian**-2,
        )
        off_angle_param.k = angle_k.in_units_of(
            unit.kilocalorie_per_mole * unit.radian**-2
        )

    edit_ff.to_file(output_dir + "/" + output_ff)


if __name__ == "__main__":
    main()
