# Ponia-OniaPhoton

This package is mean to be run using UL MINIAOD (for now 2018 is tested OK)

* Setup: (it has being tested on 10_6_16 should run in any of the recent cmssw releases)

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_10_6_16
cd CMSSW_10_6_16/src/
cmsenv
git clone https://github.com/xuliyan/chi-analysis-miniaod.git Ponia/OniaPhoton
scram b -j8

```

* Run: (use your favorite input sample)

```
cmsRun Ponia/OniaPhoton/test/run-chic-miniaod.py (for chic reconstruction)
```

In test directory you can find other examples to run over data and mc.

