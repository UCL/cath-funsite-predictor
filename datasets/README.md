
Generation of features for domains for Funsite predictor
========================================================

### 1) For each domain PDB file (e.g. 12asA00.pdb), the following files need to be generated:

a) generate [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/DSSP_3.html) file (12asA00.dssp)
```
dsspcmbi 12asA00.pdb 12asA00.dssp
```

b) generate [NACCESS](http://www.bioinf.manchester.ac.uk/naccess/) file (12asA00.rsa)
```
naccess 12asA00.pdb
```

c) Generate [centre of mass](https://github.com/rasbt/protein-science/tree/master/scripts-and-tools/center_of_mass) file (12asA00.center_of_mass)
```
python center_of_mass.py 12asA00.pdb -i ATOM > 12asA00.center_of_mass
```

d) generate [Speedfill](https://www.ebi.ac.uk/thornton-srv/software/SURFNET/) file (12asA00.clefts.depth)
```
speedfill -f 12asA00.pdb -d -ntop 5 -outdir $OUTDIR
```

e) generate [FOLDX ALASCAN](http://foldxsuite.crg.eu/command/AlaScan) file (12asA00_AS.fxout)
```
foldx --command=AlaScan --vdwDesign=0 --pdb=12asA00.pdb --pdb-dir=$PDB_DIR --output-dir=$FOLDXOUTDIR
```

f) generate centrality measures using [bio3d](https://cran.r-project.org/web/packages/bio3d/index.html) (12asA00)
```
Rscript bio3d.R $PDB_DIR 12asA00.pdb $OUTDIR
```

g) generate protein geometric parameters using [PSAIA](http://bioinfo.zesoi.fer.hr/index.php/en/10-category-en-gb/tools-en/19-psaia-en) (12ASa00.psaia)
```
psa.sh psa.cfg $PSAIA_PDB_PATH
```

h) generate PSSM feature file (12asA00.pssmfeatures)

I.  get the ATOM fasta sequence for domain (12asA00.faa)
II. generate PSSM for domain using [PSI-BLAST](http://www.biology.wustl.edu/gcg/psiblast.html)
```
psiblast -query 12asA00.faa -db nr.db -num_iterations=3 -evalue=0.001 -out_ascii_pssm 12asA00.pssm -inclusion_ethresh=0.001 -num_threads 4 -use_sw_tback
```
III. generate PSSM features
```
python generate_pssm_features.py 12asA00.pssm > 12asA00.pssmfeatures
```

i) generate DISTMAP(all-against-all pair-wise inter-residue distances) file (12asA00.distmap)
```
perl generate_distmap_and_other_features_for_dompdb.pl 12asA00.pdb $WORKDIR $OUTDIR
```

### 2) Generate features for domains

```
perl process_funsite-features_for_domain_list.pl 4.2 $DOMAINLIST $RESULTDIR $FUNSITE_DATASET_OUTFILE
```
