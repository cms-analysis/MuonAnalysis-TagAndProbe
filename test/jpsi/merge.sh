#!/bin/bash
F=$1;
if [[ "$F" == "" ]]; then echo "Usage: merge.sh file[.root]  (will merge file_*.root into file.root)"; exit 1; fi;
F="${F/.root/}";
echo "Will merge $(ls ${F}_[A-Z]*.root) into $F.root"
hadd -f $F.root ${F}_[A-Z]*.root > /dev/null 
