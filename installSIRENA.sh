#/bin/bash

# script to install SIRENA under SIXTE distribution

# Last stable functionalities of SIRENA are always inside SIXTE package.
# However if you want to install a new/experimental feature (only present
# in this repository) follow these steps:

# check SIXTE environ variable is set

if [ -z ${SIXTE+x} ]; then
    echo "SIXTE variable is unset"
    echo "Please define SIXTE var according to the SIXTE installation"
else
    echo "SIXTE is set to '$SIXTE'"
    rsync -av libsixt/  $SIXTE/../sixt/libsixt/
    rsync -av tools/ $SIXTE/../sixt/tools/
    cd $SIXTE/../sixt/
    make install
    echo "SIRENA installation is complete"
fi
