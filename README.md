Installation:

```

cmsrel CMSSW_11_2_0
cd CMSSW_11_2_0/src/
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline l1t-integration-CMSSW_11_2_0
git cms-merge-topic -u cms-l1t-offline:l1t-integration-v105.5
git cms-addpkg L1Trigger/L1TMuon
git clone https://github.com/cms-l1t-offline/L1Trigger-L1TMuon.git L1Trigger/L1TMuon/data
git cms-addpkg L1Trigger/L1TCalorimeter
git clone https://github.com/cms-l1t-offline/L1Trigger-L1TCalorimeter.git L1Trigger/L1TCalorimeter/data
git cms-checkdeps -A -a
git remote add panoskatsoulis https://github.com/panoskatsoulis/cmssw.git
git cms-merge-topic panoskatsoulis:fix_empty_bmtf_ntuples
git remote add pallabidas https://github.com/pallabidas/cmssw.git
git cms-merge-topic pallabidas:l1t-integration-test-Run3NTuplizer
cd L1Trigger
git clone -b L1boosted_L1TOffline git@github.com:pallabidas/Run3Ntuplizer.git
cd ..
USER_CXXFLAGS="-Wno-error=reorder -Wno-unused-variable" scram b -j 8

cd L1Trigger/Run3Ntuplizer/test

## to run zero bias data
cmsRun testL1TCaloSummary-ZeroBias.py

## to run boosted ggHbb events
cmsRun testL1TCaloSummary-ggHbb.py


```
