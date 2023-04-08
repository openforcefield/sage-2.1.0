import click
import numpy as np
from openff.toolkit.typing.engines.smirnoff.forcefield import ForceField
from simtk import unit


@click.command()
@click.option(
    "-ff_ref",
    "--ff_ref",
    "ff1",
    type=click.STRING,
    default="forcefield/force-field.offxml",
)
@click.option(
    "-ff",
    "--ff",
    "ff2",
    type=click.STRING,
    default="result/optimize/force-field.offxml",
)
def main(ff1, ff2):
    ff_1 = ForceField(ff1, allow_cosmetic_attributes=True)
    ff_2 = ForceField(ff2, allow_cosmetic_attributes=True)
    angle_diff = []
    kval_diff = []
    param_angle = []
    smirks_a = []
    # Angles
    j = 0
    for i, item in enumerate(ff_1.get_parameter_handler("Angles").parameters):
        if (
            ff_1.get_parameter_handler("Angles")
            .parameters[item.smirks]
            .attribute_is_cosmetic("parameterize")
        ):
            param_angle.append(item.id)
            smirks_a.append(item.smirks)
            angle_diff.append(
                (
                    ff_2.get_parameter_handler("Angles").parameters[item.smirks].angle
                    - item.angle
                )
            )
            kval_diff.append(
                (
                    ff_2.get_parameter_handler("Angles").parameters[item.smirks].k
                    - item.k
                )
            )
            print(
                item.id,
                "{: 6.2f}".format(angle_diff[i - j].value_in_unit(unit.degrees)),
                "degrees, change in value  ",
                "{: .0%}".format(
                    angle_diff[i - j].value_in_unit(unit.degrees)
                    / item.angle.value_in_unit(unit.degrees)
                ),
                "    ",
                "{: 6.2f}".format(
                    kval_diff[i - j].value_in_unit(kval_diff[i - j].unit)
                ),
                kval_diff[i - j].unit,
                ", change is  ",
                "{: .0%}".format(
                    kval_diff[i - j].value_in_unit(kval_diff[i - j].unit)
                    / item.k.value_in_unit(kval_diff[i - j].unit)
                ),
            )
        else:
            j += 1
    if param_angle:
        print(
            param_angle[np.argmax(np.abs(angle_diff))],
            smirks_a[np.argmax(np.abs(angle_diff))],
            "    ",
            max(angle_diff, key=abs),
        )
        print(
            param_angle[np.argmax(np.abs(kval_diff))],
            smirks_a[np.argmax(np.abs(kval_diff))],
            "    ",
            max(kval_diff, key=abs),
        )

    # Bonds
    bond_diff = []
    kval_b_diff = []
    param_bond = []
    smirks_b = []
    j = 0
    for i, item in enumerate(ff_1.get_parameter_handler("Bonds").parameters):
        if (
            ff_1.get_parameter_handler("Bonds")
            .parameters[item.smirks]
            .attribute_is_cosmetic("parameterize")
        ):
            param_bond.append(item.id)
            smirks_b.append(item.smirks)
            bond_diff.append(
                (
                    ff_2.get_parameter_handler("Bonds").parameters[item.smirks].length
                    - item.length
                ).__abs__()
            )
            kval_b_diff.append(
                (
                    ff_2.get_parameter_handler("Bonds").parameters[item.smirks].k
                    - item.k
                ).__abs__()
            )
            print(
                item.id,
                "{: 6.2f}".format(bond_diff[i - j].value_in_unit(unit.angstrom)),
                "angstroms, change in value is ",
                "{: .0%}".format(
                    bond_diff[i - j].value_in_unit(unit.angstrom)
                    / item.length.value_in_unit(unit.angstrom)
                ),
                "{: 6.2f}".format(
                    kval_b_diff[i - j].value_in_unit(kval_b_diff[i - j].unit)
                ),
                kval_b_diff[i - j].unit,
                "{: .0%}".format(
                    kval_b_diff[i - j].value_in_unit(kval_b_diff[i - j].unit)
                    / item.k.value_in_unit(kval_b_diff[i - j].unit)
                ),
            )
        else:
            j += 1
    if param_bond:
        print(
            param_bond[np.argmax(np.abs(bond_diff))],
            smirks_b[np.argmax(np.abs(bond_diff))],
            "   ",
            max(bond_diff, key=abs),
        )
        print(
            param_bond[np.argmax(np.abs(kval_b_diff))],
            smirks_b[np.argmax(np.abs(kval_b_diff))],
            "   ",
            max(kval_b_diff, key=abs),
        )

    # Proper Torsions
    tor_diff = []
    param_tor = []
    smirks_t = []
    largediff = []
    print("kvalues in kcal/mol")
    j = 0
    for i, item in enumerate(ff_1.get_parameter_handler("ProperTorsions").parameters):
        if (
            ff_1.get_parameter_handler("ProperTorsions")
            .parameters[item.smirks]
            .attribute_is_cosmetic("parameterize")
        ):
            param_tor.append(item.id)
            smirks_t.append(item.smirks)
            k_list1 = item.k
            k_list2 = (
                ff_2.get_parameter_handler("ProperTorsions").parameters[item.smirks].k
            )
            klist_diff = [
                "%.2f"
                % (k_list2[it] - k_list1[it]).value_in_unit(unit.kilocalorie_per_mole)
                for it in range(len(k_list1))
            ]
            if all(
                iem.value_in_unit(unit.kilocalorie_per_mole) != 0.0 for iem in k_list1
            ):
                klist_diff_percent = [
                    "{: .0%}".format(((k_list2[it] - k_list1[it]) / k_list1[it]))
                    for it in range(len(k_list1))
                ]
            else:
                klist_diff_percent = "N/A"
            tor_diff.append(klist_diff)
            print(item.id, tor_diff[i - j])
            print("    ", klist_diff_percent)
            if klist_diff_percent != "N/A" and any(
                abs(float(itm.strip("%"))) >= 50 for itm in klist_diff_percent
            ):
                largediff.append(
                    [
                        item.id,
                        item.smirks,
                        klist_diff_percent,
                        "ff_1 params: ",
                        [
                            "%.2f"
                            % k_list1[it].value_in_unit(unit.kilocalorie_per_mole)
                            for it in range(len(k_list1))
                        ],
                        ", ff_2 params:",
                        [
                            "%.2f"
                            % k_list2[it].value_in_unit(unit.kilocalorie_per_mole)
                            for it in range(len(k_list1))
                        ],
                    ]
                )
        else:
            j += 1
    print("\n")
    print("torsion params list with >= 50 % change in any k val")
    for item in largediff:
        print(item)
    # Improper Torsions
    tor_diff = []
    param_tor = []
    smirks_t = []
    largediff = []
    print("kvalues in kcal/mol")
    j = 0
    for i, item in enumerate(ff_1.get_parameter_handler("ImproperTorsions").parameters):
        if (
            ff_1.get_parameter_handler("ImproperTorsions")
            .parameters[item.smirks]
            .attribute_is_cosmetic("parameterize")
        ):
            param_tor.append(item.id)
            smirks_t.append(item.smirks)
            k_list1 = item.k
            k_list2 = (
                ff_2.get_parameter_handler("ImproperTorsions").parameters[item.smirks].k
            )
            klist_diff = [
                "%.2f"
                % (k_list2[it] - k_list1[it]).value_in_unit(unit.kilocalorie_per_mole)
                for it in range(len(k_list1))
            ]
            if all(
                iem.value_in_unit(unit.kilocalorie_per_mole) != 0.0 for iem in k_list1
            ):
                klist_diff_percent = [
                    "{: .0%}".format(((k_list2[it] - k_list1[it]) / k_list1[it]))
                    for it in range(len(k_list1))
                ]
            else:
                klist_diff_percent = "N/A"
            tor_diff.append(klist_diff)
            print(item.id, tor_diff[i - j])
            print("    ", klist_diff_percent)
        else:
            j += 1


if __name__ == "__main__":
    main()
