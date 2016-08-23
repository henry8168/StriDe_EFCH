# StriDe_EFCH
EFCH replaced old FM-extension .

Before installing this version of StriDe, TBB library is needed.
( https://www.threadingbuildingblocks.org )

step:
  1. ./autogen.sh 
  2. ./configure
  3. Add "-ltbb" into LIBS on "StriDe_EFCH/stride/MAKEFILE".
  4. make
