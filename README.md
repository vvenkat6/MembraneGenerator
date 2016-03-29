
This is a pipeline that creates a membrane file and inserts a protein structure into it if given by the user.
This membrane can then be used for molecular dynamics simulation in Amber package.
This script need python 2.7 to work
If you have different versions of python installed, the following can be run to use python2.7:
scl enable python27 bash

Dependencies:   1. Requires vmd, chimera and R installed in bin
                2. hacked amberlipid2charmm.py script that converts to amber and recognizes POP residue (instead of POPC)
                3. rpy2 package installed on python
                4. bio3d package installed on R

