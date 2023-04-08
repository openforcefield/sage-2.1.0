# Sage-2.1.0

Updates to openff-2.0.0 which include improved chemical perception for sulfonamides and phosphates, extended training set coverage which now contains training data from Gen2 as well as Gen1 datasets, and improved fitting procedures with the use of a physically intuitive starting point from modified seminario and including dihedral deviations in optimized geometry targets, and optimizing impropers as well. Contributions for this release include changes from @pavankum, @trevorgokey, @chapincavender, @jthorton and valuable feedback from team @openforcefield. 

Changes include:
 - A better starting point for the angles and bonds from the modified Seminario Method (https://doi.org/10.1021/acs.jctc.7b00785) as implemented in QUBEKit (https://github.com/qubekit/QUBEKit) (https://pubs.acs.org/doi/10.1021/acs.jcim.8b00767). Initial idea proposed by @jthorton and supporting analysis here, https://github.com/jthorton/MSM_QCArchive.
 - Included lot of new parameters for delocalized chemistries from @chapincavender's work which would assign the same parameters for atoms that are in a delocalized configuration.
 - Improved typing for sulfonamides and phosphates from @pavankum and @trevorgokey with changes in angle parameters and additional periodicities for certain torsion parameters.
 - Introduced new parameters for bridgehead nitrogens and 
 - 
