RG 22 fev. 2008
---------------

Also here, most of integer*2 arrays have been replaced by integer*4 ones.
Except in reijkl.f and ai.f where it is the interface with the ijkl part of
the program.

Also, the call of the subroutine morue.f have been disabled. A mysterious
error :

ld: morue.o r_value (0xb71f18c9) field of relocation entry 115 in section (__TEXT,__text) out of range

was generated during the compilation because of this subroutine. It seems that
the bd (diagonalization) calculation seems to be fine with that…

Dimensions have, of course, been increased. Check the pshf.prm file :)

