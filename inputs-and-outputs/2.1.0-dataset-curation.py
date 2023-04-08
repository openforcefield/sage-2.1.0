import copy
import functools
import logging
import random
from collections import defaultdict
from multiprocessing import Pool
from tempfile import NamedTemporaryFile

from openff.qcsubmit.results import (
    TorsionDriveResultCollection,
    OptimizationResultCollection,
)
from openff.qcsubmit.results.filters import (
    ConnectivityFilter,
    RecordStatusFilter,
    UnperceivableStereoFilter,
    HydrogenBondFilter,
    ElementFilter,
    ConformerRMSDFilter,
    ResultRecordFilter,
)
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils.exceptions import UnassignedMoleculeChargeException
from qcportal import FractalClient
from qcportal.models import TorsionDriveRecord
from qcportal.models.records import RecordStatusEnum
from tqdm import tqdm

explicit_ring_torsions = [
    "t15",
    "t44",
    "t49",
    "t80",
    "t84",
    "t85",
    "t87",
    "t87a",
    "t111",
    "t117",
    "t128",
    "t129",
    "t136",
    "t137",
    "t140",
]


class ChargeCheckFilter(ResultRecordFilter):
    def _filter_function(self, result, record, molecule) -> bool:
        # Some of the molecules fail charging with am1bccelf10 either
        # because of no bccs or failed conformer generation, sometimes it
        # cannot be captured with just the cmiles present in the record
        # metadata, so reading from file and checking it
        can_be_charged = True

        molecule = copy.deepcopy(molecule)
        molecule._conformers = [molecule.conformers[0]]

        try:
            with NamedTemporaryFile(suffix=".sdf") as file:
                molecule.to_file(file.name, "SDF")
                molecule.from_file(file.name)
                molecule.assign_partial_charges(partial_charge_method="am1bccelf10")

        except UnassignedMoleculeChargeException:
            can_be_charged = False

        return can_be_charged


def label_and_tag_torsion_ids(record_and_molecule, force_field, parameter_types):
    record, molecule = record_and_molecule
    full_labels = force_field.label_molecules(molecule.to_topology())[0]

    parameter_ids = set()

    for parameter_type in parameter_types:
        parameter_labels = full_labels[parameter_type]

        for indices, parameter in parameter_labels.items():
            if isinstance(record, TorsionDriveRecord) and {*indices[1:3]} != {
                *record.keywords.dihedrals[0][1:3]
            }:
                continue

            if parameter.id not in explicit_ring_torsions:
                # some general parameters overlap with in-ring torsions and
                # there are many torsion scans from Gen1 sets that have
                # in-ring torsions and we want to exclude them in training
                # as they result in higher k values unless the parameters
                # have smirks explicitly for an in-ring torsion. It is to be
                # noted that training on in-ring torsions is needed to
                # properly model puckering in rings with hetero atoms
                if (
                    molecule.get_bond_between(indices[0], indices[1]).is_in_ring()
                    and molecule.get_bond_between(indices[1], indices[2]).is_in_ring()
                    and molecule.get_bond_between(indices[2], indices[3]).is_in_ring()
                ):
                    continue
            rdmol = molecule.to_rdkit()
            num_heavy_atoms = rdmol.GetNumHeavyAtoms()
            parameter_ids.add((parameter.id, record.id, num_heavy_atoms))

    return [*parameter_ids]


def get_parameter_distribution(training_set, parameter_types, force_field):
    coverage = defaultdict(int)
    parameter_records = defaultdict(list)
    heavy_atom_count = defaultdict(list)

    with Pool(8) as pool:
        for parameter_ids in tqdm(
            pool.imap(
                functools.partial(
                    label_and_tag_torsion_ids,
                    force_field=force_field,
                    parameter_types=parameter_types,
                ),
                training_set.to_records(),
            ),
            total=training_set.n_results,
        ):
            for parameter_id, record, n_heavy in parameter_ids:
                coverage[parameter_id] += 1
                parameter_records[parameter_id].append(record)
                heavy_atom_count[parameter_id].append(n_heavy)

    # Sorting the record ids in descending order of heavy atom count
    # so that when the torsions are capped picking up the torsion
    # scans with a higher number of heavy atoms (that was also
    # another objective of Gen2 sets, to generate larger molecules
    # which would have more non-bonded interactions at play. So,
    # instead of randomly picking up the records we can pick up the
    # first few with higher heavy atom count.
    for parameter_id, value in parameter_records.items():
        heavy_atom_count[parameter_id], parameter_records[parameter_id] = zip(
            *sorted(
                zip(heavy_atom_count[parameter_id], parameter_records[parameter_id]),
                reverse=True,
            )
        )
        parameter_records[parameter_id] = list(parameter_records[parameter_id])
        heavy_atom_count[parameter_id] = list(heavy_atom_count[parameter_id])

    return coverage, parameter_records


def cap_torsions_per_parameter(
    force_field, torsion_set_to_filter, cap_size, method="pick_random"
):
    coverage, parameter_records = get_parameter_distribution(
        torsion_set_to_filter, ["ProperTorsions"], force_field
    )
    tor_full = [
        "t1",
        "t2",
        "t3",
        "t4",
        "t5",
        "t6",
        "t7",
        "t8",
        "t9",
        "t10",
        "t11",
        "t12",
        "t13",
        "t14",
        "t15",
        "t16",
        "t17",
        "t18",
        "t18a",
        "t18b",
        "t19",
        "t19a",
        "t20",
        "t21",
        "t22",
        "t23",
        "t24",
        "t25",
        "t26",
        "t27",
        "t28",
        "t29",
        "t30",
        "t31",
        "t31a",
        "t32",
        "t33",
        "t34",
        "t35",
        "t36",
        "t37",
        "t38",
        "t39",
        "t40",
        "t41",
        "t42",
        "t42a",
        "t43",
        "t44",
        "t45",
        "t46",
        "t47",
        "t48",
        "t48a",
        "t49",
        "t50",
        "t51",
        "t52",
        "t53",
        "t54",
        "t55",
        "t56",
        "t57",
        "t58",
        "t59",
        "t60",
        "t61",
        "t62",
        "t63",
        "t64",
        "t65",
        "t66",
        "t67",
        "t68",
        "t69",
        "t70",
        "t71",
        "t72",
        "t73",
        "t74",
        "t75",
        "t76",
        "t77",
        "t78",
        "t79",
        "t80",
        "t81",
        "t82",
        "t82a",
        "t83",
        "t83a",
        "t84",
        "t85",
        "t86",
        "t87",
        "t87a",
        "t88",
        "t89",
        "t90",
        "t91",
        "t92",
        "t93",
        "t94",
        "t95",
        "t96",
        "t97",
        "t98",
        "t99",
        "t100",
        "t101",
        "t102",
        "t103",
        "t104",
        "t105",
        "t106",
        "t107",
        "t108",
        "t109",
        "t110",
        "t111",
        "t112",
        "t113",
        "t114",
        "t115",
        "t116",
        "t117",
        "t118",
        "t119",
        "t120",
        "t121",
        "t122",
        "t123",
        "t123a",
        "t124",
        "t125",
        "t126",
        "t127",
        "t128",
        "t129",
        "t130",
        "t131",
        "t132",
        "t133",
        "t134",
        "t135",
        "t136",
        "t137",
        "t138",
        "t138a",
        "t139",
        "t140",
        "t141",
        "t141a",
        "t141b",
        "t141c",
        "t142",
        "t143",
        "t144",
        "t145",
        "t146",
        "t147",
        "t148",
        "t149",
        "t150",
        "t151",
        "t152",
        "t153",
        "t154",
        "t155",
        "t156",
        "t157",
        "t158",
        "t159",
        "t160",
        "t161",
        "t162",
        "t163",
        "t164",
        "t165",
        "t166",
        "t167",
    ]
    tor_keys = list(coverage.keys())
    records_to_keep = []

    for key in tor_keys:
        if len(parameter_records[key]) <= cap_size:
            print(key, len(parameter_records[key]), parameter_records[key])
            records_to_keep.extend(parameter_records[key])
        else:
            print(
                key,
                min(cap_size, len(parameter_records[key])),
                parameter_records[key][:cap_size],
            )
            if method == "pick_heavy":
                records_to_keep.extend(parameter_records[key][:cap_size])
            elif method == "pick_random":
                records_to_keep.extend(random.sample(parameter_records[key], cap_size))
    print("Length of records to keep: ", len(records_to_keep))
    indices_to_delete = []
    for i, entry in enumerate(
        torsion_set_to_filter.entries["https://api.qcarchive.molssi.org:443/"]
    ):
        if entry.record_id not in records_to_keep:
            indices_to_delete.append(i)
    indices_to_delete = list(set(indices_to_delete))

    for item in sorted(indices_to_delete, reverse=True):
        torsion_set_to_filter.entries["https://api.qcarchive.molssi.org:443/"].pop(item)

    print(
        f"After capping to {cap_size} torsion scans per parameter, the "
        f"distribution of torsion scans per parameter is as follows:"
    )
    for key in tor_full:
        if key in tor_keys:
            print(key, coverage[key])
    print(tor_full)
    print(coverage)

    return torsion_set_to_filter


def main():
    logging.getLogger("openff").setLevel(logging.ERROR)
    from pathlib import Path

    Path("./data-sets").mkdir(parents=True, exist_ok=True)

    default_filters = [
        RecordStatusFilter(status=RecordStatusEnum.complete),
        ConnectivityFilter(tolerance=1.2),
        UnperceivableStereoFilter(),
    ]
    ff = ForceField(
        "../modified_initial_force_field/initial-force-field" ".offxml",
        allow_cosmetic_attributes=True,
    )

    # Pull down the main torsion drive and optimization sets and filter out any records
    # which have not completed or which inadvertently contain intra-molecular h-bonds.
    client = FractalClient()

    # SMIRNOFF Coverage torsions set inconsistent IDs, ELF failures,
    # Gen3 no-param-exotic-mol and other known errors in FB fits
    tdrecs_to_remove = [
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
    ]
    optrecs_to_remove = [
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
    ]
    # Following the strategy to keep everything from Gen2 sets and
    # augmenting them with upto 5 additional records from Gen1 and
    # other sets

    #######
    # torsion_subset_1
    #######
    torsion_subset_1 = TorsionDriveResultCollection.from_server(
        client=client,
        datasets=[
            "OpenFF Gen 2 Torsion Set 1 Roche 2",
            "OpenFF Gen 2 Torsion Set 2 Coverage 2",
            "OpenFF Gen 2 Torsion Set 3 Pfizer Discrepancy 2",
            "OpenFF Gen 2 Torsion Set 4 eMolecules Discrepancy 2",
            "OpenFF Gen 2 Torsion Set 5 Bayer 2",
            "OpenFF Gen 2 Torsion Set 6 supplemental 2",
        ],
        spec_name="default",
    )

    # Drop record ids with inconsistent optimization histories or which cause failures
    # in ForceBalance.
    torsion_subset_1.entries[client.address] = [
        entry
        for entry in torsion_subset_1.entries[client.address]
        if entry.record_id not in tdrecs_to_remove
    ]

    torsion_subset_1 = torsion_subset_1.filter(
        HydrogenBondFilter(method="baker-hubbard"),
        *default_filters,
        ElementFilter(
            # The elements supported by SMIRNOFF
            # Excluding Iodine here since we don't have Iodine torsions and any record with iodine is tainted on the
            # datasets listed above because of the auxiliary basis set issue
            allowed_elements=["H", "C", "N", "O", "S", "P", "F", "Cl", "Br"]
        ),
    )
    with open("data-sets/td-subset-1-datasets.json", "w") as file:
        file.write(torsion_subset_1.json())

    #######
    # torsion_subset_2
    #######
    torsion_subset_2 = TorsionDriveResultCollection.from_server(
        client=client,
        datasets=[
            "SMIRNOFF Coverage Torsion Set 1",
            "OpenFF Group1 Torsions",
            "OpenFF Group1 Torsions 2",
            "OpenFF Group1 Torsions 3",
            "Pfizer discrepancy torsion dataset 1",
            "OpenFF Gen3 Torsion Set v1.0",
            "OpenFF Amide Torsion Set v1.0",
            "OpenFF WBO Conjugated Series v1.0",
            "OpenFF DANCE 1 eMolecules t142 v1.0",
        ],
        spec_name="default",
    )
    torsion_subset_2.entries[client.address] = [
        entry
        for entry in torsion_subset_2.entries[client.address]
        if entry.record_id not in tdrecs_to_remove
    ]

    torsion_subset_2 = torsion_subset_2.filter(
        HydrogenBondFilter(method="baker-hubbard"),
        *default_filters,
        ElementFilter(
            # The elements supported by SMIRNOFF
            # Excluding Iodine here since we don't have Iodine torsions and any record with iodine is tainted on the
            # datasets listed above because of the auxiliary basis set issue
            allowed_elements=["H", "C", "N", "O", "S", "P", "F", "Cl", "Br"]
        ),
    )
    with open("data-sets/td-subset-2-datasets-before-capping.json", "w") as file:
        file.write(torsion_subset_2.json())

    torsion_subset_2 = cap_torsions_per_parameter(
        force_field=ff, torsion_set_to_filter=torsion_subset_2, cap_size=5
    )

    with open("data-sets/td-subset-2-datasets-after-capping.json", "w") as file:
        file.write(torsion_subset_2.json())

    torsion_subset_1.entries["https://api.qcarchive.molssi.org:443/"].extend(
        torsion_subset_2.entries["https://api.qcarchive.molssi.org:443/"]
    )
    torsion_set = torsion_subset_1
    # t126 don't have any matches even after adding all these datasets,
    # so adding the following record manually for t126,
    # which got filtered due to H-bond,
    # {"type": "torsion", "record_id": "2703402", "cmiles": "[H:6][C:2]([H:7])([H:8])[C:1](=[O:3])[O:5][O:4][H:9]",
    # "inchi_key": "KFSLWBXXFJQRDL-UHFFFAOYNA-N"}

    # removing duplicate entries, there were 18, 1318 to 1300
    torsion_set.entries["https://api.qcarchive.molssi.org:443/"].append(
        {
            "type": "torsion",
            "record_id": "2703402",
            "cmiles": "[H:6][C:2]([H:7])([H:8])[C:1](=[O:3])[O:5][O:4][H:9]",
            "inchi_key": "KFSLWBXXFJQRDL-UHFFFAOYNA-N",
        }
    )
    unique = {
        each["record_id"]: each
        for each in torsion_set.entries["https://api.qcarchive.molssi.org:443/"]
    }.values()
    torsion_set.entries["https://api.qcarchive.molssi.org:443/"] = list(unique)
    torsion_set = torsion_set.filter(ChargeCheckFilter())
    with open("data-sets/td-set-for-fitting-charge-check-2.1.0.json", "w") as file:
        file.write(torsion_set.json())
    # torsion_set = TorsionDriveResultCollection.parse_file("data-sets/td-set-for-fitting-2.1.0.json")

    #######
    # opt_subset_1: Gen2 sets without iodine containing mols
    #######
    opt_subset_1 = OptimizationResultCollection.from_server(
        client=FractalClient(),
        datasets=[
            "OpenFF Gen 2 Opt Set 1 Roche",
            "OpenFF Gen 2 Opt Set 2 Coverage",
            "OpenFF Gen 2 Opt Set 3 Pfizer Discrepancy",
            "OpenFF Gen 2 Opt Set 4 eMolecules Discrepancy",
            "OpenFF Gen 2 Opt Set 5 Bayer",
        ],
        spec_name="default",
    )
    opt_subset_1 = opt_subset_1.filter(
        ElementFilter(
            # The elements supported by SMIRNOFF
            # Excluding Iodine here since we don't have Iodine torsions and any record with iodine is tainted on the
            # datasets listed above because of the auxiliary basis set issue
            # New sets added below in opt_subset_2 has Iodine containing molecules that are safe
            allowed_elements=["H", "C", "N", "O", "S", "P", "F", "Cl", "Br"]
        ),
    )

    opt_subset_1.entries[client.address] = [
        entry
        for entry in opt_subset_1.entries[client.address]
        if entry.record_id not in optrecs_to_remove
    ]

    #######
    # opt_subset_2: Gen2 sets with iodine containing mols and extra protomers
    #######
    opt_subset_2 = OptimizationResultCollection.from_server(
        client=FractalClient(),
        datasets=[
            "OpenFF Gen2 Optimization Dataset Protomers v1.0",
            "OpenFF Iodine Chemistry Optimization Dataset v1.0",
        ],
        spec_name="default",
    )

    opt_subset_2 = opt_subset_2.filter(
        ElementFilter(
            # The elements supported by SMIRNOFF
            allowed_elements=["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]
        ),
    )

    ## Filtering Gen2 sets to cap conformers at 10 based on greedy approach
    # to pick 10 confs that differ by RMSD
    opt_subset_1.entries["https://api.qcarchive.molssi.org:443/"].extend(
        opt_subset_2.entries["https://api.qcarchive.molssi.org:443/"]
    )
    optimization_set = opt_subset_1
    optimization_set = optimization_set.filter(
        RecordStatusFilter(status=RecordStatusEnum.complete),
        ConnectivityFilter(tolerance=1.2),
        UnperceivableStereoFilter(),
        ConformerRMSDFilter(max_conformers=12),
        ChargeCheckFilter(),
    )
    with open("data-sets/opt-subsets-1-and-2.json", "w") as file:
        file.write(optimization_set.json())
    #######
    # opt_subset_3: Gen 1 and Aniline para sets for more molecules
    #######
    opt_subset_3 = OptimizationResultCollection.from_server(
        client=FractalClient(),
        datasets=[
            "OpenFF Optimization Set 1",
            "SMIRNOFF Coverage Set 1",
            "OpenFF Aniline Para Opt v1.0",
        ],
        spec_name="default",
    )
    opt_subset_3 = opt_subset_3.filter(
        ElementFilter(
            # The elements supported by SMIRNOFF
            # Excluding Iodine here since we don't have Iodine torsions and any record with iodine is tainted on the
            # datasets listed above because of the auxiliary basis set issue
            # New sets added below in opt_subset_2 has Iodine containing molecules that are safe
            allowed_elements=["H", "C", "N", "O", "S", "P", "F", "Cl", "Br"]
        ),
    )

    opt_subset_3.entries[client.address] = [
        entry
        for entry in opt_subset_3.entries[client.address]
        if entry.record_id not in optrecs_to_remove
    ]
    opt_subset_3 = opt_subset_3.filter(
        RecordStatusFilter(status=RecordStatusEnum.complete),
        ConnectivityFilter(tolerance=1.2),
        UnperceivableStereoFilter(),
        ConformerRMSDFilter(max_conformers=12),
        ChargeCheckFilter(),
    )
    with open("data-sets/opt-subset-3.json", "w") as file:
        file.write(opt_subset_3.json())

    optimization_set.entries["https://api.qcarchive.molssi.org:443/"].extend(
        opt_subset_3.entries["https://api.qcarchive.molssi.org:443/"]
    )

    with open("data-sets/opt-set-for-fitting-charge-check-2.1.0.json", "w") as file:
        file.write(optimization_set.json())

    print("done!")


if __name__ == "__main__":
    main()
