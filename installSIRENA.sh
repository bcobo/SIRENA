#/bin/bash

# script to install SIRENA under SIXTE distribution

# check SIXTE environ variable is set

if [ -z ${SIXTE+x} ]; then
    echo "SIXTE variable is unset"
    echo "Please define SIXTE var according to the SIXTE installation"
else
    echo "SIXTE is set to '$SIXTE'"
    rsync -av libsixt/  $SIXTE/../sixt/libsixt/
    rsync -av tools/ $SIXTE/../sixt/tools/
    cd $SIXTE/../sixt/libsixt
    make install
    echo "SIRENA installation is complete"
fi
