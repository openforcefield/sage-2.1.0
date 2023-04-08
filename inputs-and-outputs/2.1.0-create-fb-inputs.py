import json
import os.path
from pathlib import Path

from openff.bespokefit.optimizers.forcebalance import ForceBalanceInputFactory
from openff.bespokefit.schema.fitting import OptimizationSchema, OptimizationStageSchema
from openff.bespokefit.schema.optimizers import ForceBalanceSchema
from openff.bespokefit.schema.smirnoff import (
    AngleHyperparameters,
    AngleSMIRKS,
    BondHyperparameters,
    ImproperTorsionHyperparameters,
    ProperTorsionHyperparameters,
    BondSMIRKS,
    ProperTorsionSMIRKS,
)
from openff.bespokefit.schema.targets import (
    OptGeoTargetSchema,
    TorsionProfileTargetSchema,
)
from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
from openff.qcsubmit.results.filters import (
    SMARTSFilter,
    SMILESFilter,
)
from openff.toolkit.typing.engines.smirnoff import ForceField


def main():
    Path("./schemas/optimizations/").mkdir(parents=True, exist_ok=True)
    tag = "fb-fit"
    port_number = 55387
    Path("./" + tag).mkdir(parents=True, exist_ok=True)
    torsion_training_set = TorsionDriveResultCollection.parse_file(
        tag + "/data-sets/" + "td-set-for-fitting-2.1.0.json"
    )
    # remove charging-failures
    torsion_training_set.entries["https://api.qcarchive.molssi.org:443/"] = [
        entry
        for entry in torsion_training_set.entries[
            "https://api.qcarchive.molssi.org:443/"
        ]
        if entry.record_id
        not in [
            "6098580",
            "2703504",
            "2703505",
            "18045478",
            "2703253",
            "2703343",
            "2703386",
            "2703439",
            "2703449",
            "2703545",
            "2703546",
            "2703616",
            "35045000",
            "18433638",
            "2703106",
            "2703093",
            "2703125",
            "2703122",
            "2703104",
            "2703087",
            "2703086",
            "2703103",
            "3504499",
            "2703124",
            "2703107",
            "2703090",
            "2703100",
            "2703098",
            "2703084",
            "2703097",
            "35044997",
        ]
    ]

    # Filtering out unusual chemistries
    smarts_to_exclude = [
        "[#8+1:1]=[#7:2]",
        "[#15:1]=[#6:2]",
        "[#16+1:1]~[*:2]",
        "[*:1]=[#15:2]-[#7:3]~[*:4]",
        "[#17:1]~[#1:2]",
        "[#16-1:1]~[#15:2]",  # No BCC
        "[#7:1]=[#7:2]#[#7:3]",
        "[#16-1:1]~[#16:2]",  # ELF10 failure
    ]
    smiles_to_exclude = ["[S](=[N-])(=O)", "[H]C([H])(C#N)N=N#N", "O[O-]"]
    torsion_training_set = torsion_training_set.filter(
        SMARTSFilter(smarts_to_exclude=smarts_to_exclude),
        SMILESFilter(smiles_to_exclude=smiles_to_exclude),
    )

    optimization_training_set = OptimizationResultCollection.parse_file(
        tag + "/data-sets/" + "opt-set-for-fitting-2.1.0.json"
    )

    optimization_training_set = optimization_training_set.filter(
        SMARTSFilter(smarts_to_exclude=smarts_to_exclude),
        SMILESFilter(smiles_to_exclude=smiles_to_exclude),
    )
    optimization_training_set.entries["https://api.qcarchive.molssi.org:443/"] = [
        entry
        for entry in optimization_training_set.entries[
            "https://api.qcarchive.molssi.org:443/"
        ]
        if entry.record_id
        not in [
            "2002949",
            "2002950",
            "18433638",
            "18433906",
            "2002933",
            "2002934",
            "2002937",
            "2003047",
            "2003043",
            "95602295",
            "95602250",
            "18433502",
            "18434090",
            "2002949",
            "2002950",
            "18433638",
            "18433906",
            "2002933",
            "2002934",
            "2002937",
            "2003047",
            "2003043",
            "95602295",
            "95602250",
            "18433502",
            "18434090",
            "18433675",
            "18433675",
            "2003404",
            "2002930",
            "2002929",
            "2002979",
        ]  # Some more charging failures
    ]
    # to pick initial values and parameters to optimize
    # enter
    custom_force_field = (
        "../modified_initial_force_field/msm_output_opt_set_10"
        "/msm_output_opt_set_10.offxml"
    )
    initial_force_field = ForceField(custom_force_field, allow_cosmetic_attributes=True)

    # Define the parameters to train
    with open(tag + "/data-sets/opt-set-angle-smirks.json") as file:
        angle_smirks = json.load(file)
    with open(tag + "/data-sets/opt-set-bond-smirks.json") as file:
        bond_smirks = json.load(file)
    with open(tag + "/data-sets/td-set-torsion-smirks.json") as file:
        torsion_smirks = json.load(file)

    # a16, a17, a27, a35
    linear_angle_smirks = [
        "[*:1]~[#6X2:2]~[*:3]",  # a16
        "[*:1]~[#7X2:2]~[*:3]",  # a17
        "[*:1]~[#7X2:2]~[#7X1:3]",  # a27
        "[*:1]=[#16X2:2]=[*:3]",
    ]  # a35, this one anyways doesn't have a training target for ages

    target_parameters = [
        *[
            AngleSMIRKS(smirks=smirks, attributes={"k", "angle"})
            if smirks not in linear_angle_smirks
            else AngleSMIRKS(smirks=smirks, attributes={"k"})
            for smirks in angle_smirks["Angles"]
        ],
        *[
            BondSMIRKS(smirks=smirks, attributes={"k", "length"})
            for smirks in bond_smirks["Bonds"]
        ],
        *[
            ProperTorsionSMIRKS(
                smirks=smirks,
                attributes={
                    f"k{i + 1}"
                    for i in range(
                        len(
                            initial_force_field.get_parameter_handler("ProperTorsions")
                            .parameters[smirks]
                            .k
                        )
                    )
                },
            )
            for smirks in torsion_smirks["ProperTorsions"]
        ],
    ]

    # Define the full schema for the optimization.
    optimization_schema = OptimizationSchema(
        id=tag,
        initial_force_field=os.path.abspath(custom_force_field),
        # Define the optimizer / ForceBalance specific settings.
        stages=[
            OptimizationStageSchema(
                optimizer=ForceBalanceSchema(
                    max_iterations=50,
                    step_convergence_threshold=0.01,
                    objective_convergence_threshold=0.1,
                    gradient_convergence_threshold=0.1,
                    n_criteria=2,
                    initial_trust_radius=-1.0,
                    finite_difference_h=0.01,
                    extras={
                        "wq_port": str(port_number),
                        "asynchronous": "True",
                        "search_tolerance": "0.1",
                        "backup": "0",
                        "retain_micro_outputs": "0",
                    },
                ),
                # Define the torsion profile targets to fit against.
                targets=[
                    TorsionProfileTargetSchema(
                        reference_data=torsion_training_set,
                        energy_denominator=1.0,
                        energy_cutoff=8.0,
                        extras={"remote": "1"},
                    ),
                    OptGeoTargetSchema(
                        reference_data=optimization_training_set,
                        weight=0.01,
                        extras={"batch_size": 30, "remote": "1"},
                        bond_denominator=0.05,
                        angle_denominator=5.0,
                        dihedral_denominator=10.0,
                        improper_denominator=10.0,
                    ),
                ],
                # Define the parameters to refit and the priors to place on them.
                parameters=target_parameters,
                parameter_hyperparameters=[
                    AngleHyperparameters(priors={"k": 100, "angle": 5}),
                    BondHyperparameters(priors={"k": 100, "length": 0.1}),
                    ProperTorsionHyperparameters(priors={"k": 5}),
                    ImproperTorsionHyperparameters(priors={"k": 5}),
                ],
            )
        ],
    )

    with open(
        os.path.join("./schemas", "optimizations", f"{optimization_schema.id}.json"),
        "w",
    ) as file:
        file.write(optimization_schema.json())

    # Generate the ForceBalance inputs
    ForceBalanceInputFactory.generate(
        os.path.join(optimization_schema.id),
        optimization_schema.stages[0],
        ForceField(optimization_schema.initial_force_field),
    )


if __name__ == "__main__":
    main()
