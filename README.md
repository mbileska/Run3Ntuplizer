Installation:

```

cmsrel CMSSW_13_3_0
cd CMSSW_13_3_0/src/
cmsenv
git cms-init
git cms-merge-topic -u pallabidas:l1boosted-133X
cd L1Trigger
git clone -b L1boosted_133X git@github.com:pallabidas/Run3Ntuplizer.git
cd ..
scram b -j 12

cd L1Trigger/Run3Ntuplizer/test

## to run zero bias data
cmsRun testL1TCaloSummary-ZeroBias.py

## to run boosted ggHbb events
cmsRun testL1TCaloSummary-ggHBB.py


```
