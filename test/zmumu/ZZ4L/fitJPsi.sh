for X in data2012 mc2012; do 
    for B in pt_abseta pt_abseta2; do
        for I in PF Tight2012noIP Tight2012withIP ZZLoose Glb TMA; do
            echo ~/sh/cmsBSub  -8nh fitMuonID_JPsi.py $X $I $B; 
        done
    done
done
