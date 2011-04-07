#!/bin/bash
F=$1; shift;
if [[ "$F" == "" ]]; then echo "Usage: merge.sh file[.root]  (will merge file_*.root into file.root)"; exit 1; fi;
F="${F/.root/}";
FILES=$(ls ${F}_[A-Z]*.root);
if [[ "X$1" != "X" ]]; then
    FILES=$(ls ${F}_[A-Z]*.root | grep $*);
fi;
echo "Will merge $FILES into $F.root"
hadd -f $F.root $FILES > /dev/null 
