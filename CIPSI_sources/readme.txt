Bugfixes/correction (RV 01/2016) :

Added mrec $vpol namelist (cval) needed for allocate memory
new : all blocks in unit=8 (ijkl,cval and scf) are read and stored in memory at the beginning

ERROR not fixed in pshf : dname, fname on atom.f  are uninitialized !!!!!
set to blanc as temporary fix

pshf/cale.f
"alpha" replace by "amix"

pshf/dafile.f
ante initialized to "  pshf"

pshf/atom.f
loop on dname changed (from i=1,7 to i=1,6)
dfg2 initialazed to blanc

pshf/ dafile.f;daclose.f;ptgrp.f;cale.f
corrected some declarations :
logical ymono
character*8 aname,direct,title

increase dimension of iodb to 5000 (allow more iteration)

pshf/scf.f
replaced hdiag subroutine by a lapack call (diagonaliser.f90)