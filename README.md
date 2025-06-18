# Framework to generate private B physics samples


## Installation

Set up the environment

```
cmsrel CMSSW_10_6_28
cd CMSSW_10_6_28/src
cmsenv
git cms-init
```

Install the tool

```
git clone git@github.com:CPVAnalysis/CPVGen.git
```

## Private production on CRAB

Activate the proxy

```
voms-proxy-init --voms cms --valid 186:00
```

Then, 

```
cd CPVGen/crab
```

Run the tool as

```
python CRABLauncher.py --pl <prod_label> --nevents <nevents> <--dogenonly> <--dosubmit>
```

with 
   * nevents the requested number of events at analysis level.
   * --dogenonly the option to run the GEN step only instead of the full GEN->miniAOD chain.

