# SIRENA
Sub-package of files developed @IFCA(CSIC-UC) for the onboard event processing of Athena/XIFU data
Release versions integrated in the simulator SIXTE/SIRENA http://www.sternwarte.uni-erlangen.de/research/sixte/index.php

To install the last version of this package:

1. Install SIXTE following instructions at http://www.sternwarte.uni-erlangen.de/research/sixte and define SIXTE environ variable
 according to instructions.
2. Go to SIRENA installation dir:
   > cd $SIXTE
   > cd ../sixt
   > git clone https://github.com/bcobo/SIRENA.git
   or
   > unzip SIRENA-master.zip (if you downloaded the ZIP file)

3. Recompile SIXTE following SIXTE instructions:
   > make install
