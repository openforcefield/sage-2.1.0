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

# Details of changes (still updating this section...)
 - A better starting point for the angles and bonds from the modified Seminario Method (https://doi.org/10.1021/acs.jctc.7b00785) as implemented in QUBEKit (https://github.com/qubekit/QUBEKit) (https://pubs.acs.org/doi/10.1021/acs.jcim.8b00767). Initial work by @jthorton and supporting analysis here, https://github.com/jthorton/MSM_QCArchive. This is helpful in avoiding manual corrections to parameters and brings the force constants of bonds and angles to physically intuitive domain. Since this is multi-dimensional optimization we get solutions on the pareto-optimal surface which resulted in some double bonds having lower force constants than single bonds. Modified seminario uses projection of Hessian data to evaluate the force constants and the mean of the force constants (and lengths) for the parameters was taken as the starting point for forcebalance optimization. Here is an example of bond parameter force constant changes from b34 to b40, which includes bonds between nitrogens (similar large changes can be observed in other bond and angle parameters):
 
 |     Bond smarts     |      Sage 2.0.0     |      Sage 2.1.0     |
|:-------------------:|:-------------------:|:-------------------:|
|     (b34 - b40)     | k (kcal/mol/ang**2) | k (kcal/mol/ang**2) |
|   "[#7:1]-[#7:2]"   |         845         |         578         |
| "[#7X3:1]-[#7X2:2]" |         830         |         620         |
| "[#7X2:1]-[#7X2:2]" |         675         |         473         |
|   "[#7:1]:[#7:2]"   |         732         |         662         |
|   "[#7:1]=[#7:2]"   |         698         |         1089        |
| "[#7+1:1]=[#7-1:2]" |         766         |         2440        |
|   "[#7:1]#[#7:2]"   |         760         |         3237        |

 - Included lot of new parameters for chemistries that involve delocalized charges from @chapincavender's work which would assign the same parameters for atom substructures that are in a delocalized configuration (https://github.com/chapincavender/protein-param-fit/blob/attenuation/test-delocalized-charge-assignments.py). Functional groups targeted were:
     - Amidinium (added t18b)
     - Guanidinium (added b13a, t87a)
     
     ![image](https://user-images.githubusercontent.com/16142894/231004432-ff59cd33-9b4a-44fa-888c-0f156504cd28.png)
     
     - Carboxylate (added t18a, t19a, t31a, t42a, t48a)
     
     ![image](https://user-images.githubusercontent.com/16142894/231004590-dbb74503-70c9-4820-a7c3-d08c2280545f.png)

     - Nitro (added t82a, t83a)
     
     ![image](https://user-images.githubusercontent.com/16142894/231004529-1e1dabe1-ae23-4ad5-8e95-df3a98723c06.png)

     - Phosphate (no change needed)
     - Sulfate (no change needed)
 - Improved typing for sulfonamides and phosphates from @pavankum and @trevorgokey with changes in angle parameters and additional periodicities for certain torsion parameters. The changed parameters are b57a, change in a32 smirks definition to match [*]-[S]=[*] substructure, which has a mean near 100°, from a31, which has a mean near 120°.
 ![image](https://user-images.githubusercontent.com/16142894/231005046-097894a6-ff2c-4da4-a1da-d390784d245b.png)

 - Introduced new parameters for bridgehead nitrogens and bridgehead carbons where pyramidal geometry is observed in QM but the MM models it as planar.
 - Minor modifications to generalize parameters to include chemistries of interest (commonly encountered on PubChem and ChemBL) that don't have any parameters assigned and fail:
      - b43 changed from `[#8X2:1]-[#8X2:2]` to `[#8X2:1]-[#8X2,#8X1-1:2]`, to accommodate the missing bond for `( O), ( O)` in a molecule such as (`[N:1](=[O:2])[O:3][O-:4]`)
      
      ![image](https://user-images.githubusercontent.com/16142894/231008845-06584b78-0a8b-486b-bb18-91083014fde7.png)

      - b53 changed from `[#16X2:1]-[#7:2]` to `[#16X2,#16X1-1:1]-[#7:2]`, to accommodate the missing bond for `( N), ( S)` in a molecule such as 
      
      ![image](https://user-images.githubusercontent.com/16142894/231009017-ca8d2e71-33b1-4021-a9e5-eec549422966.png)

      - a18 changed from `[*:1]-[#7X4,#7X3,#7X2-1:2]-[*:3]` to `[*:1]~[#7X4,#7X3,#7X2-1:2]~[*:3]` to accommodate angle assignment for  `( P), ( N), ( O)` in a molecule such as (`[H:34][c:13]1[c:12]([c:11]([c:10]([c:15]([c:14]1[H:35])[H:36])[P:8](=[O:9])([N+:16](=[O:17])[O-:18])[O:7][C:6]([H:30])([H:31])[C:5]([H:28])([H:29])[N+:2]([C:1]([H:19])([H:20])[H:21])([C:3]([H:22])([H:23])[H:24])[C:4]([H:25])([H:26])[H:27])[H:32])[H:33]`)
      
       ![image](https://user-images.githubusercontent.com/16142894/231007723-54c9d086-9c3f-46e2-8725-8c5add572b62.png)
      
      - t51 changed from `[*:1]-[#6X4:2]-[#7X3:3]-[*:4]` to `[*:1]~[#6X4:2]-[#7X3:3]~[*:4]` to accommodate `( N), ( N), ( C), ( H)` torsions in molecules such as (`[H:7][C:1]([H:8])([H:9])[N+:2](=[N:3][C:4]([H:10])([H:11])[O:5][H:12])[O-:6]`) 
      
      ![image](https://user-images.githubusercontent.com/16142894/231007298-02ebe095-4364-4082-a990-291524646bab.png)
 
      - t130 changed from `[*:1]-[#7X4,#7X3:2]-[#7X4,#7X3:3]-[*:4]` to `[*:1]~[#7X4,#7X3:2]-[#7X4,#7X3:3]~[*:4]`, to accommodate `( C), ( N), ( N), ( O)` and `( H), ( N), ( N), ( O)` torsions in a molecule such as (`[H:11][C:1]([H:12])([H:13])[N:2]([C:3](=[N:4][H:14])[N:5]([H:15])[N+:6](=[O:7])[O-:8])[N:9]=[O:10]`) 
      
      ![image](https://user-images.githubusercontent.com/16142894/231009518-4e99ae4d-a8f7-4f50-963c-39d273d158d4.png)

      - t138a with the smirks `"[*:1]~[#7X2:2]-[#7X4:3]~[*:4]"` to accommodate torsion assignment for `( C), ( N), ( N), ( C)` in a molecule such as (`[H:26][C:10](=[C:9](/[C:8](=[N:7]/[N+:4]([C:5]([H:20])([H:21])[H:22])([C:6]([H:23])([H:24])[H:25])[C:3]([H:18])([H:19])[C:2]([H:17])([C:1]([H:14])([H:15])[H:16])[O:13][H:31])/[O-:12])[C:11]([H:28])([H:29])[H:30])[H:27]`)
      
      ![image](https://user-images.githubusercontent.com/16142894/231008143-be95151d-9a15-49bc-90cd-8d927c67af5e.png)

      - t161 changed from `"[*:1]~[#7X3:2]-[#15:3]~[*:4]"` to `"[*:1]~[#7:2]-[#15:3]~[*:4]"` to make it more general, to accommodate `( C), ( N), ( P), ( O)` and `( C), ( N), ( P), ( N)` torsions in a molecule such as (`[H:19][C:7]1=[N:6][P:4](=[O:5])([O:3][C:2]([C:1]1([H:15])[H:16])([H:17])[H:18])[N:8]([C:9]([H:20])([H:21])[C:10]([H:22])([H:23])[Cl:11])[C:12]([H:24])([H:25])[C:13]([H:26])([H:27])[Cl:14]`) 
      
      ![image](https://user-images.githubusercontent.com/16142894/231010478-3d339032-875d-4946-9f8b-f7af1d847025.png)

