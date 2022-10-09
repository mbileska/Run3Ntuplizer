Installation:

```

cmsrel CMSSW_12_4_0
cd CMSSW_12_4_0/src/
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline l1t-integration-CMSSW_12_4_0
git cms-merge-topic -u cms-l1t-offline:l1t-integration-v134
git cms-addpkg L1Trigger/L1TCalorimeter
git clone https://github.com/cms-l1t-offline/L1Trigger-L1TCalorimeter.git L1Trigger/L1TCalorimeter/data
git cms-checkdeps -A -a
git remote add pallabidas https://github.com/pallabidas/cmssw.git
git cms-merge-topic pallabidas:l1t-integration-test-09Oct
cd L1Trigger
git clone -b L1boosted_L1TOffline git@github.com:pallabidas/Run3Ntuplizer.git
cd ..
scram b -j 12

cd L1Trigger/Run3Ntuplizer/test

## to run zero bias data
cmsRun testL1TCaloSummary-ZeroBias.py

## to run boosted ggHbb events
cmsRun testL1TCaloSummary-ggHbb.py


```
