# Sage-2.1.0

Updates to openff-2.0.0 which include improved chemical perception for sulfonamides and phosphates, extended training set coverage which now contains training data from Gen2 as well as Gen1 datasets, and improved fitting procedures with the use of a physically intuitive starting point from modified seminario and including dihedral deviations in optimized geometry targets, and optimizing impropers as well. Contributions for this release include changes from @pavankum, @trevorgokey, @chapincavender, @jthorton and valuable feedback from team @openforcefield. 

Changes can be broadly classified as: 
 - Chemical typing related
 - sulfonamides 
 - phosphates 
 - bridgehead nitrogens 
 - bridgehead carbons 
 - groups with delocalized charges
     - Fitting procedure related
     - use of physically intuitive bonds and angles from modified-seminario as a starting point
     - data-driven values for improper torsions
     - including dihedral deviations from optimized geometries to better resolve torsion parameters
     - broader coverage of parameters with extended training targets

Parameters with modifications:
```
        <!-- modification: CC b13a added guanadinium specific bond-->
        <!-- modification: PB b43 changed from [#8X2:1]-[#8X2:2] to [#8X2:1]-[#8X2,#8X1-1:2] -->
        <!-- modification: PB b53 changed from [#16X2:1]-[#7:2] to [#16X2,#16X1-1:1]-[#7:2] -->
        <!-- modification: PB b57a child parameter to separate out double bondish/double bonds from b57  -->
        <!-- modification: PB-TG a18 changed from [*:1]-[#7X4,#7X3,#7X2-1:2]-[*:3] to [*:1]~[#7X4,#7X3,#7X2-1:2]~[*:3] -->
        <!-- modification: PB-TG a18a child parameter to separate out some ring matches from a18  -->
        <!-- modification: PB-TG a22a child parameter to separate out some ring matches from a22  -->
        <!-- modification: PB-TG a32 smirks definition changed to drain [*]-[S]=[*] matches from a31  -->
        <!-- modification: CC added carboxylate torsions t18a to address delocalized charges -->
        <!-- modification: CC added amidinium torsion t18b to address delocalized charges -->
        <!-- modification: CC added carboxylate torsions t19a to address delocalized charges -->
        <!-- modification: CC added carboxylate torsions t31a to address delocalized charges -->
        <!-- modification: CC added carboxylate torsions t42a to address delocalized charges -->
        <!-- modification: CC added carboxylate torsions t48a to address delocalized charges -->
        <!-- modification: PB t51 changed from [*:1]-[#6X4:2]-[#7X3:3]-[*:4] to [*:1]~[#6X4:2]-[#7X3:3]~[*:4] -->
        <!-- modification: CC added nitro torsion t82a to address delocalized charges -->
        <!-- modification: CC added nitro torsion t83a to address delocalized charges -->
        <!-- modification: CC added guanidinium torsion t87a to address delocalized charges -->
        <!-- modification: PB-TG t123a added with a periodicity of 3 to better match that chemistry -->
        <!-- modification: PB-TG t124, changed periodicity from 1 with phase_0, to periodicity of 2 with phase_0, and 3 with phase_0 -->
        <!-- modification: PB t130 changed from [*:1]-[#7X4,#7X3:2]-[#7X4,#7X3:3]-[*:4] to [*:1]~[#7X4,#7X3:2]-[#7X4,#7X3:3]~[*:4] -->
        <!-- modification: PB t138a as a child parameter to include [#7X2]-[#7X4] chemistry, other general force fields can parameterize this -->
        <!-- modification: PB new parameter t141a for bridgehead Nitrogen chemistry based on t134 and tweaking the central Nitrogen to be 7x3 -->
        <!-- modification: PB new parameter t141b for bridgehead Nitrogen chemistry based on t138 and tweaking the central Nitrogen to be 7x3 -->
        <!-- modification: PB t141c, bridgehead carbons with a heteroatom neighbor has non-planar geometry, so a specific torsion term for that -->
        <!-- modification: PB-TG t159, additional periodicities of 1 with phase_0, and 2 with phase_0 -->
        <!-- modification: PB-TG t160, additional periodicity of 1 with phase_0 -->
        <!-- modification: PB t161 changed from "[*:1]~[#7X3:2]-[#15:3]~[*:4]" to "[*:1]~[#7:2]-[#15:3]~[*:4]" to make it more general -->
        <!-- modification: PB additional periodicities for t143 and t157 based on QM profiles for sulfonamides -->
```

Details of changes (still updating this section...):
 - A better starting point for the angles and bonds from the modified Seminario Method (https://doi.org/10.1021/acs.jctc.7b00785) as implemented in QUBEKit (https://github.com/qubekit/QUBEKit) (https://pubs.acs.org/doi/10.1021/acs.jcim.8b00767). Initial work by @jthorton and supporting analysis here, https://github.com/jthorton/MSM_QCArchive. This is helpful in avoiding manual corrections to parameters and brings the force constants of bonds and angles in physically intuitive domain. Since this is multi-dimensional optimization we get solutions on the pareto-optimal surface which resulted in some double bonds having lower force constants than single bonds. Modified seminario uses Hessian data to evaluate the force constants and the mean of the force constants (and lengths) for the parameters was taken as the starting point for forcebalance optimization.
 - Included lot of new parameters for chemistries that involve delocalized charges from @chapincavender's work which would assign the same parameters for atoms that are in a delocalized configuration (https://github.com/chapincavender/protein-param-fit/blob/attenuation/test-delocalized-charge-assignments.py). Functional groups targeted were:
     - Amidinium
     - Carboxylate
     - Nitro
     - Phosphate
     - Sulfate
 - Improved typing for sulfonamides and phosphates from @pavankum and @trevorgokey with changes in angle parameters and additional periodicities for certain torsion parameters.
 - Introduced new parameters for bridgehead nitrogens and bridgehead carbons where pyramidal geometry is seen in QM.
 - 
