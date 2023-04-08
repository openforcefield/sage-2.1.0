# Input files to create the forcebalance inputs and the optimization run output
    -  2.1.0-check-elf10-charging.py: checking whether the targets generated can charge with AM1BCC-ELF10 (included the same in dataset-curation)
    -  2.1.0-check-parameter-coverage.py: checking which valence parameters match to the target molecules and tag them with parameterize (excludes some linear angles/torsions)
    -  2.1.0-create-fb-inputs.py: script that creates forcebalance inputs by reading the record information in data-sets directory
    -  2.1.0-create_msm_ff.py: script that would create a starting forcefield based on the hessians of target optimization records using modified-seminario method
    -  2.1.0-dataset-curation.py: script that is used to curate the training datasets, Gen2 + Gen1 datasets were used in the training for a broader coverage
    -  2.1.0-forcefield-diff.py: utility script to check the difference between any two forcefield files
    -  2.1.0-remove_cosmetic_attributes.py: utility script to remove cosmetic attributes from a forcefield file
    -  data-sets/ : directory that contains the opt-geo and torsion profile targets information (the smirks files are generic files used to generate inputs, parameters to optimize were tagged using check-parameter-coverage script)
    -  fb-fit/ : forcebalance inputs created and the final output
    -  msm_starting_point/ : output of the create_msm_ff script, which is used as starting point for the forcebalance run

