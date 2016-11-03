The required files are `waterbox.psf` and `init.pdb`.

## Initial PDB (a.pdb)

Use a single water as the seed to create a.pdb.
```
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1
ATOM   1232  OH2 TIP3W   1      14.108  17.970  10.538  1.00  0.00      WT1  O
ATOM   1233  H1  TIP3W   1      13.664  18.339   9.735  1.00  0.00      WT1  H
ATOM   1234  H2  TIP3W   1      13.532  18.267  11.261  1.00  0.00      WT1  H
END
```
The box size on the first line does not matter.


## Initial PSF (water.pdb and water.psf)

Write `a.pgn`
```
package require psfgen
topology top_all27_prot_lipid.inp
segment U {pdb a.pdb}
coordpdb a.pdb U	 
writepdb water.pdb 
writepsf water.psf
```

Run
```
vmd -dispdev text -e a.pgn
```

Type `exit` to quit VMD.

This will create `water.pdb` and `water.psf`.

See 
http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node6.html


## Solvate the system

Need `water.pdb` and `water.psf`.

To place the solute in a 40x40x40 box
In VMD, Tk Console
```
package require solvate
solvate water.psf water.pdb -minmax {{0 0 0} {40 40 40}} -o waterbox
```
This creates `waterbox.pdb` and `waterbox.psf`

More Options of the solvate module
https://github.com/thatchristoph/vmd-cvs-github/blob/master/plugins/solvate/solvate.tcl

http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/node7.html#SECTION00043200000000000000

### Tweaking waterbox.pdb

1. Change the chain name (at the end of line) from "U  " to "WT1" for atoms 1 to 3.

2. Change `TIP3U` to `TIP3W` for the first three lines.

### Tweaking waterbox.psf

1. Change `5 !NTITLE` to `4 !NTITLE`

2. Delete `REMARKS segment U {....`

3. Change atoms 1 to 3 from `U    1...` to `WT1  1...`
   for the first three lines after `!NATOM`.

### Cleaning up

`rm waterbox.log`

## Equilibration


### Editing `equil.conf`

Change the cell size.  For a 40x40x40 box, it should be
```
cellBasisVector1    40.0    0.   0.0
cellBasisVector2     0.0  40.0   0.0
cellBasisVector3     0.0    0   40.0
cellOrigin          20.0  20.0  20.0
```

Change the number of minimizations and productions and run `equil.conf`

```
namd2 +p2 equil.conf
```

This creates `wb.pdb`


### Creating a PDB from the final coordinates

1.
```
vmd waterbox.pdb 
```

2. Load `wb.coor`.
Right clicking the PDB, and selecting `Load Data Into Molecule`.

3. Right clicking PDB and selecting `Save Coordinates` 
First `1`, Last `1`, save name `init.pdb`.

### Cleaning up

```
rm -f wb*
```
