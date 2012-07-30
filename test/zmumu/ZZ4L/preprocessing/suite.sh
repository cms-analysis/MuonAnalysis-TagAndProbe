if test \! -f $1; then echo "Usage: $0 src/file.root"; exit 1; fi
DIR=$(basename $1 .root).d
mkdir $DIR || exit 1;
cd $DIR && \
root.exe -b -l -q ../$1 '../slimTree.cxx(0)'              && \
root.exe -b -l -q tnpZ_slim_step0.root ../subTree.C       && \
root.exe -b -l -q subTree_IsoMu24.root ../addEAIso.cxx+   && \
rm subTree_IsoMu24.root                                   && \
root.exe -b -l -q tnpZ_withEAIso.root  ../addTriggerEmu.cxx         && \
root.exe -b -l -q tnpZ_withTriggerEmu.root  '../slimTree.cxx(1)'      && \
rm tnpZ_withTriggerEmu.root                                   && \
mv -v tnpZ_slim_step1.root ../$(basename $1 .root)_forID.root -v

#root.exe -b -l -q tnpZ_withEAIso.root ../loadMVAlibs.cxx ../addMVAIso.cxx+ && \
#rm tnpZ_withEAIso.root                                               && \
