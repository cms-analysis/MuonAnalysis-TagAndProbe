scenario="$1"; if [[ "$scenario" == "" ]]; then scenario="data2011"; fi;

IDS=""
IDS="$IDS PF SIP4_from_PF"
IDS="$IDS Glb TMA"
IDS="$IDS ZZLoose"
IDS="$IDS Tight2012noIP Tight2012withIP"
#IDS="$IDS SIP4only_from_PF SIP4_from_PF_and_SIP4only"
IDS="$IDS PFIso40_from_PF_and_SIP4"
IDS="$IDS PFIso20DB_from_PF_and_SIP4"
IDS="$IDS PFIso12DB_from_PF_and_SIP4"
IDS="$IDS PFIso40_from_Tight2012withIP"
IDS="$IDS PFIso20DB_from_Tight2012withIP"
IDS="$IDS PFIso12DB_from_Tight2012withIP"
if echo $scenario | grep -q 2011; then
    IDS="$IDS HLT_Mu17_from_PF_and_SIP4_and_PFIso40"
    IDS="$IDS HLT_Mu17_from_PF_and_SIP4_and_PFIso40_and_2011A"
    IDS="$IDS HLT_Mu17_from_PF_and_SIP4_and_PFIso40_and_2011B"
    IDS="$IDS HLT_Mu8_from_PF_and_SIP4_and_PFIso40"
    IDS="$IDS HLT_Mu8_from_PF_and_SIP4_and_PFIso40_and_2011A"
    IDS="$IDS HLT_Mu8_from_PF_and_SIP4_and_PFIso40_and_2011B"
else
    IDS="$IDS HLT_Mu17_from_PF_and_SIP4_and_PFIso40 "
    IDS="$IDS HLT_Mu17_from_PF_and_SIP4_and_PFIso40_and_PostTS "
    IDS="$IDS HLT_Mu8_from_PF_and_SIP4_and_PFIso40 "
    IDS="$IDS HLT_Mu8_from_PF_and_SIP4_and_PFIso40_and_PostTS "
    IDS="$IDS HLT_OrMu17_from_PF_and_SIP4_and_PFIso40"
    IDS="$IDS HLT_OrMu17_from_PF_and_SIP4_and_PFIso40_and_PostTS"
    IDS="$IDS HLT_OrMu8_from_PF_and_SIP4_and_PFIso40"
    IDS="$IDS HLT_OrMu8_from_PF_and_SIP4_and_PFIso40_and_PostTS"
    IDS="$IDS HLT_TkMu17_from_PF_and_SIP4_and_PFIso40 "
    IDS="$IDS HLT_TkMu17_from_PF_and_SIP4_and_PFIso40_and_PostTS "
    IDS="$IDS HLT_TkMu8_from_PF_and_SIP4_and_PFIso40 "
    IDS="$IDS HLT_TkMu8_from_PF_and_SIP4_and_PFIso40_and_PostTS "
fi;
if echo $scenario | grep -q echo; then  echo $IDS | perl -npe 's/\s+/\n/g'; exit 0; fi;
for X in $IDS; do
    for B in pt_abseta pt_abseta2 eta vtx; do
        echo ~/sh/cmsBSub -8nh  fitMuonID.py $scenario $X $B
    done
done



