# cath-funsite-predictor

This is the code base for generating protein functional site predictions using the CATH Funsite predictor.

This git repository contains the following:

* `datasets`
  This contains all the datasets mentioned in the Funsite manuscript along with jupyter notebooks showing the steps for generating them.
  Some dataset files are very large, these can be downloaded from [**ftp://orengoftp.biochem.ucl.ac.uk/cath/supplementary-materials/2020-cath-funsite-predictor/**](ftp://orengoftp.biochem.ucl.ac.uk/cath/supplementary-materials/2020-cath-funsite-predictor/)

* `funsite_models`
  This contains the Funsite model generation, prediction scripts and benchmark results using the dataset files that are reported in the Funsite manuscript.

* `scripts`
  This contains all scripts used for generating Funsite models and predictions.


### Python environment for Funsite predictor
```
virtualenv -p python3.6 FunsiteEnv
source FunsiteEnv/bin/activate

# install dependencies
pip install -r scripts/FunsiteEnv_requirements.txt
```

### Getting CATH FunFam assignments for proteins

To get [CATH](http://www.cathdb.info/wiki) FunFam assignments for proteins, one can use the [online search tool](http://www.cathdb.info/search/by_sequence) for few sequences or the [cath-genomescan tool](https://github.com/UCLOrengoGroup/cath-tools-genomescan) for getting assignments for large sequence datasets.

### External softwares

This repository uses external bioinformatics tools that are not written and maintained by the authors of this project. If you use the results of these tools, please reference the relevant papers.

#### [GroupSim](https://compbio.cs.princeton.edu/specificity/)
Characterization and Prediction of Residues Determining Protein Functional Specificity. Capra JA and Singh M (2008). Bioinformatics, 24(13): 1473-1480, 2008.

#### [Scorecons](https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/valdar/scorecons_server.pl)
Scoring residue conservation. Valdar WSJ (2002)
Proteins: Structure, Function, and Genetics. 43(2): 227-241, 2002.

#### [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html)
Kabsch, W, and C Sander. 1983. ‘DSSP: Definition of Secondary Structure of Proteins given a Set of 3D Coordinates’. Biopolymers 22: 2577–2637.

#### [NACCESS](http://www.bioinf.manchester.ac.uk/naccess/)
Hubbard, S J. 1992. ‘NACCESS: Program for Calculating Accessibilities’. Department of Biochemistry and Molecular Biology, University College of London.

#### [Speedfill](https://www.ebi.ac.uk/thornton-srv/software/SURFNET/)
Laskowski, R A. 1995. ‘SURFNET: A Program for Visualizing Molecular Surfaces, Cavities, and Intermolecular Interactions’. J. Mol. Graph. 13 (5): 307-308,323-330.

#### [FOLDX ALASCAN](http://foldxsuite.crg.eu/command/AlaScan)
Schymkowitz, Joost, Jesper Borg, Francois Stricher, Robby Nys, Frederic Rousseau, and Luis Serrano. 2005. ‘The FoldX Web Server: An Online Force Field’. Nucleic Acids Res. 33 (Web Server issue): W382-8.

#### [bio3d](https://cran.r-project.org/web/packages/bio3d/index.html)
Skjærven, Lars, Shashank Jariwala, Xin-Qiu Yao, Julien Idé, and Barry J Grant. 2016. ‘The Bio3D Project: Interactive Tools for Structural Bioinformatics’. Biophys. J. 110 (3): 379a.

#### [PSAIA](http://bioinfo.zesoi.fer.hr/index.php/en/10-category-en-gb/tools-en/19-psaia-en)
Mihel, Josip, Mile Sikić, Sanja Tomić, Branko Jeren, and Kristian Vlahovicek. 2008. ‘PSAIA - Protein Structure and Interaction Analyzer’. BMC Struct. Biol. 8 (April): 21.

#### [PSI-BLAST](http://www.biology.wustl.edu/gcg/psiblast.html)
Altschul, S.F., Madden, T.L., Schäffer, A.A., Zhang, J., Zhang, Z., Miller, W. and Lipman, D.J., 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic acids research, 25(17), pp.3389-3402.

#### [Centre of mass](https://github.com/rasbt/protein-science/tree/master/scripts-and-tools/center_of_mass)

### References

The most recent papers describing the CATH protein structure database and CATH FunFams:

1. [CATH: expanding the horizons of structure-based functional annotations for genome sequences.](https://doi.org/10.1093/nar/gky1097)

2. [Functional classification of CATH superfamilies: a domain-based approach for protein function annotation](https://doi.org/10.1093/bioinformatics/btv398)

3. [CATH FunFHMMer web server: protein functional annotations using functional family assignments](https://doi.org/10.1093/nar/gkv488)
