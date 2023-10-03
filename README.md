Installation:

```

cmsrel CMSSW_12_4_3
cd CMSSW_12_4_3/src/
cmsenv
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git remote add pallabidas https://github.com/pallabidas/cmssw.git
git cms-merge-topic pallabidas:l1boosted-regioninfo
cd L1Trigger
git clone -b L1boosted_124X_MLtagger git@github.com:pallabidas/Run3Ntuplizer.git
cd ..
scram b -j 12

cd L1Trigger/Run3Ntuplizer/test

## to run boosted ggHbb events
cmsRun testL1TCaloSummary-ggHbb.py


```
