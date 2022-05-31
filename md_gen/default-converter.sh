#!/bin/bash

mdconvert -o frame0.xtc frame0.dcd
# check if frame0.xtc exists
if [[ -f frame0.xtc ]]; then
    # make sure it has
    size_dcd=$(stat -c %s frame0.dcd)
    size_xtc=$(stat -c %s frame0.xtc)
    safety_factor=$((size_xtc*5))
    if [[ size_dcd -gt safety_factor ]]; then
        echo Warning: trajectory file size was not large enough to plausibly imply successful conversion
        echo sizes were frame0.dcd: $size_dcd, frame0.xtc: $size_xtc
        echo Doing NOTHING!
    else
        rm frame0.dcd
    fi
else
    echo WARNING: frame0.xtc does not exist. Conversion failed.
fi