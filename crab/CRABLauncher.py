import os
import sys
from os import path
import random


def getOptions():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='CRAB Launcher', add_help=True)

  parser.add_argument('--pl'        , type=str, dest='pl'         , help='production label'                 ,                      default='V00_v00')
  parser.add_argument('--nevents'   , type=str, dest='nevents'    , help='requested events (analysis level)',                      default='20000')
  parser.add_argument('--year'      , type=str, dest='year'       , help='year (2018, 2022, or 2024)'       ,                      default='2018')
  parser.add_argument('--dogenonly' ,           dest='dogenonly'  , help='run GEN step only'                , action='store_true', default=False)
  parser.add_argument('--dosubmit'  ,           dest='dosubmit'   , help='submit to CRAB'                   , action='store_true', default=False)
  parser.add_argument('--doresubmit',           dest='doresubmit' , help='resubmit to CRAB'                 , action='store_true', default=False)

  return parser.parse_args()


class CMSDriver(object):
    def __init__(self, era, global_tag, beamspot_conditions):
        self.era = era
        self.global_tag = global_tag
        self.beamspot_conditions = beamspot_conditions


class CRABLauncher(object):
  def __init__(self, options):
    self.opt = options
    for k,v in sorted(vars(opt).items()):
      setattr(self,k,v)

    # some fixed parameters
    self.nevents_perminiaod = 5

    # filter efficiency obtained from private GEN production test_trgmu_v0 with trigger muon filter and without acceptance cuts
    #eff = 1.86e-4 

    # obtained with test_fragment_v2 (no cond on trgmu mother + kaon pt > 0.5 GeV)
    eff = 0.00012451599

    # acceptance of kaons (pt > 0.5 GeV |eta| < 2.5) otained with test_trgmu_v0
    #acc = 0.54 

    # acceptance of kaons (pt > 0.5 GeV |eta| < 2.5) otained with test_fragment_v2
    acc = 0.97 

    self.eff_filter = eff * acc
    self.eff_nanoaod = 0.3 #FIXME too low?

    #self.GT = '106X_upgrade2018_realistic_v11_L1v1'
    #self.GT = '102X_upgrade2018_realistic_v11'
    #self.CMSSW = 'CMSSW_10_2_28_patch1' 

    #TODO add sqrt(s) as parameter, as it is different between Run2 and Run3
    self.pl += '_{}'.format(self.year)

    self.year = int(self.year)
    if self.year not in [2018, 2022, 2024]:
        raise RuntimeError('Incorrect year. Please choose among [2018, 2022, 2024]')

    if self.year != 2018 and not self.dogenonly:
        raise RuntimeError('The DIGI->miniAOD steps are not yet supported for years other than 2018')

    if self.year == 2018:
        self.driver = CMSDriver(
                era = 'Run2_2018',
                global_tag = '106X_upgrade2018_realistic_v11_L1v1',
                beamspot_conditions = 'Realistic25ns13TeVEarly2018Collision'
                )
    elif self.year == 2022:
        self.driver = CMSDriver(
                era = 'Run3',
                global_tag = '124X_mcRun3_2022_realistic_v12',
                beamspot_conditions = 'Realistic25ns13p6TeVEarly2022Collision'
                )
    elif self.year == 2024:
        self.driver = CMSDriver(
                era = 'Run3_2024',
                global_tag = '140X_mcRun3_2024_realistic_v26',
                beamspot_conditions = 'Realistic'
                )

    if self.year == 2018:
      self.fragment_name = 'BsToPhiPhiTo4K_cfi.py'
      self.CMSSW = 'CMSSW_10_6_28'
    else:
      self.fragment_name = 'BsToPhiPhiTo4K_Run3_cfi.py'
      self.CMSSW = 'CMSSW_14_0_21'



  def createOuputDir(self):
    dirname = './production/{}'.format(self.pl)
    if not path.exists(dirname):
      os.system('mkdir -p {}'.format(dirname))


  def createCRABConfig(self):
    if not self.dogenonly:
      self.nevents_togenerate = float(self.nevents) / (float(self.eff_filter) * float(self.eff_nanoaod))
    else:
      self.nevents_togenerate = float(self.nevents) / float(self.eff_filter)

    self.nevents_perjob = int(self.nevents_perminiaod / float(self.eff_filter))

    if self.nevents_togenerate > self.nevents_perjob:
      self.njobs = int(self.nevents_togenerate / self.nevents_perjob)
    else:
      self.nevents_perjob = int(self.nevents_togenerate)
      self.njobs = 1

    config = [
      'from CRABClient.UserUtilities import config',
      'config = config()',
      '',
      'import datetime, time',
      'ts = time.time()',
      'st = datetime.datetime.fromtimestamp(ts).strftime("%Y%m%d_%H%M%S")',
      'config.General.requestName = "{pl}_BsToPhiPhiTo4K_"+ st',
      'config.General.transferOutputs = True',
      'config.General.transferLogs = True',
      'config.General.workArea = "crab_workdir"',
      '',
      'config.JobType.pluginName = "PrivateMC"',
      'config.JobType.psetName = "step1.py"',
      '{inputfiles}',
      'config.JobType.outputFiles = {outputfiles}',
      'config.JobType.scriptExe = "submitter.sh"',
      'config.JobType.disableAutomaticOutputCollection = True',
      'ncores = 1',
      'mempecore = 3500',
      'config.JobType.maxMemoryMB  = ncores*mempecore',                                                                                                                                  
      'config.JobType.allowUndistributedCMSSW = True',
      '',
      'config.Data.outputPrimaryDataset = "BsToPhiPhiTo4K"',
      'config.Data.outLFNDirBase = "/store/user/anlyon/CPVGen/{pl}"',
      'NUMEVENTSPERJOB={nevtsjob}',
      'NUMJOBS={njobs}',
      'config.Data.totalUnits = NUMEVENTSPERJOB*NUMJOBS',
      'config.Data.publication = False',
      'config.Data.unitsPerJob = NUMEVENTSPERJOB',
      'config.Data.splitting = "EventBased"',
      'config.Data.ignoreLocality = False',
      '',
      'config.Site.storageSite = "T3_CH_CERNBOX"',
      ]

    config = '\n'.join(config)
    config = config.format(
        inputfiles = 'config.JobType.inputFiles = ["../../data/FrameworkJobReport.xml", "step1.py", "../../cmsDrivers/step_DIGI.py", "../../data/list_pu_files.py", "../../cmsDrivers/step_HLT.py", "../../cmsDrivers/step_RECO.py", "../../cmsDrivers/step_MINI.py"]' if not self.dogenonly else 'config.JobType.inputFiles = ["../../data/FrameworkJobReport.xml", "step1.py"]',
        outputfiles = '["step1.root", "miniaod.root"]' if not self.dogenonly else '["step1.root"]',
        pl = self.pl,
        time = 30000, # in minutes
        nevtsjob = self.nevents_perjob,
        njobs = self.njobs,
        )

    config_filename = './production/{pl}/crab_config.py'.format(pl=self.pl)
    config_file = open(config_filename, 'w+')
    config_file.write(config)
    config_file.close()


  def createSubmitter(self):
    if not self.dogenonly:
      submitter = [
        '#!/bin/bash',
        'echo " "',
        'echo "Source new environment"',
        'source /cvmfs/cms.cern.ch/cmsset_default.sh',
        'if [ -r {cmssw}/src ] ; then',
        '  echo release {cmssw} already exists',
        'else',
        '  scram p CMSSW {cmssw}',
        'fi',
        'cd {cmssw}/src',
        'eval `scram runtime -sh`',
        'scram b',
        'cd ../..',
        'echo " "',
        'echo "content of /srv/{cmssw}/src"',
        'ls -al /srv/{cmssw}/src',
        'echo " "',
        'echo "will run step1"',
        'cmsRun -j step1.log step1.py',
        'echo "end run step1"',
        'echo "content of dir"',
        'ls -al',
        'echo " "',
        'echo "will run step2 (DIGI)"',
        'cmsRun -j step2.log step_DIGI.py',
        'echo "end run step2 (DIGI)"',
        'echo "content of dir"',
        'ls -al',
        'echo " "',
        'echo "Source new environment"',
        'source /cvmfs/cms.cern.ch/cmsset_default.sh',
        'if [ -r CMSSW_10_2_16_UL/src ] ; then',
        '  echo release CMSSW_10_2_16_UL already exists',
        'else',
        '  scram p CMSSW CMSSW_10_2_16_UL',
        'fi',
        'cd CMSSW_10_2_16_UL/src',
        'eval `scram runtime -sh`',
        'scram b',
        'cd ../..',
        'echo " "',
        'echo "will run step3 (HLT)"',
        'cmsRun -j step3.log step_HLT.py',
        'echo "end run step3 (HLT)"',
        'echo "content of dir"',
        'ls -al',
        'echo " "',
        'echo "Source new environment"',
        'source /cvmfs/cms.cern.ch/cmsset_default.sh',
        'if [ -r {cmssw}/src ] ; then',
        '  echo release {cmssw} already exists',
        'else',
        '  scram p CMSSW {cmssw}',
        'fi',
        'cd {cmssw}/src',
        'eval `scram runtime -sh`',
        'scram b',
        'cd ../..',
        'echo " "',
        'echo "will run step4 (RECO)"',
        'cmsRun -j step4.log step_RECO.py',
        'echo "end run step4 (RECO)"',
        'echo "content of dir"',
        'ls -al',
        'echo " "',
        'echo "will run step5 (MINI)"',
        'cmsRun -e -j FrameworkJobReport.xml step_MINI.py',
        'echo "end run step5 (MINI)"',
        'echo "Done"',
        ]
    else:
      submitter = [
        '#!/bin/bash',
        'echo " "',
        'echo "Source new environment"',
        'source /cvmfs/cms.cern.ch/cmsset_default.sh',
        'if [ -r {cmssw}/src ] ; then',
        '  echo release {cmssw} already exists',
        'else',
        '  scram p CMSSW {cmssw}',
        'fi',
        'cd {cmssw}/src',
        'eval `scram runtime -sh`',
        'scram b',
        'cd ../..',
        'echo "content of /srv/{cmssw}/src"',
        'ls -al /srv/{cmssw}/src',
        'echo " "',
        'echo "will run step1"',
        'cmsRun -j step1.log step1.py',
        'echo "end run step1"',
        'echo "content of dir"',
        'ls -al',
        'echo " "',
        'echo "Done"',
        ]

    submitter = '\n'.join(submitter)
    submitter = submitter.format(
        cmssw = self.CMSSW,
        )

    submitter_filename = './production/{pl}/submitter.sh'.format(pl=self.pl)
    submitter_file = open(submitter_filename, 'w+')
    submitter_file.write(submitter)
    submitter_file.close()


  def createDriver(self):

    #TODO ok for both 2018 and 2024?
    if not self.dogenonly:
      command = 'cmsDriver.py Configuration/GenProduction/python/fragment.py --fileout file:step1.root --mc --eventcontent FEVTDEBUG --datatier GEN-SIM --conditions {gt} --beamspot {bs} --step GEN,SIM --geometry DB:Extended --era {era} --python_filename step1.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n {nevts} --mc --customise_commands "from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService); randSvc.populate()"'
    else:
      command = 'cmsDriver.py Configuration/GenProduction/python/fragment.py --fileout file:step1.root --mc --eventcontent FEVTDEBUG --datatier GEN --conditions {gt} --beamspot {bs} --step GEN --geometry DB:Extended --era {era} --python_filename step1.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n {nevts} --mc --customise_commands "from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService); randSvc.populate()"'

    command = command.format(
          gt = self.driver.global_tag,
          bs = self.driver.beamspot_conditions,
          era = self.driver.era,
          nevts = self.nevents_perjob,
          )
    command_replace = ''

    if self.doresubmit or not path.exists('./production/{pl}/step1.py'.format(pl=self.pl)):
      submitter_tmp = [
        '#!/bin/bash',
        'STARTDIR=$PWD',
        'PRODDIR=$PWD/production/{pl}',
        'WORKDIR=$CMSSW_BASE/src',
        'cd $WORKDIR',
        'echo "Source new environment"',
        'source /cvmfs/cms.cern.ch/cmsset_default.sh',
        'if [ -r {cmssw}/src ] ; then',
        '  echo release {cmssw} already exists',
        'else',
        '  scram p CMSSW {cmssw}',
        'fi',
        'cd {cmssw}/src',
        'eval `scram runtime -sh`',
        'scram b',
        'cd ../..',
        'mkdir -p Configuration/GenProduction/python/',
        'cp $STARTDIR/../fragments/{frgmt} Configuration/GenProduction/python/fragment.py',
        'scram b -j 8',
        '{cmd}',
        '{cmd_rpl}',
        'cp step1.py $PRODDIR/step1.py',
        'cd $STARTDIR',
        ]

      submitter_tmp = '\n'.join(submitter_tmp)
      submitter_tmp = submitter_tmp.format(
          pl = self.pl,
          frgmt = self.fragment_name,
          cmd = command,
          cmd_rpl = command_replace,
          cmssw = self.CMSSW,
          )

      submitter_tmp_file = open('submitter_tmp.sh', 'w+')
      submitter_tmp_file.write(submitter_tmp)
      submitter_tmp_file.close()
      
      os.system('sh submitter_tmp.sh > log.txt')

      #os.system('rm submitter_tmp.sh') 
      os.system('rm log.txt') 


  def createStep1(self):
    part1 = [ 
       "import FWCore.ParameterSet.Config as cms",
       " ",
       "from Configuration.Eras.Era_{era}_cff import {era}",
       " ",
       "process = cms.Process('GEN', {era})",
       " ",
       "process.load('Configuration.StandardSequences.Services_cff')",
       "process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')",
       "process.load('FWCore.MessageService.MessageLogger_cfi')",
       "process.load('Configuration.EventContent.EventContent_cff')",
       "process.load('SimGeneral.MixingModule.mixNoPU_cfi')",
       "process.load('Configuration.StandardSequences.GeometryRecoDB_cff')",
       "process.load('Configuration.StandardSequences.MagneticField_cff')",
       "process.load('Configuration.StandardSequences.Generator_cff')",
       "process.load('IOMC.EventVertexGenerators.VtxSmeared{bs}_cfi')",
       "process.load('GeneratorInterface.Core.genFilterSummary_cff')",
       "process.load('Configuration.StandardSequences.EndOfProcess_cff')",
       "process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')",
       " ",
       "process.maxEvents = cms.untracked.PSet(",
       "    input = cms.untracked.int32({nevents})",
       ")",
       " ",
       "process.source = cms.Source('EmptySource')",
       " ",
       "process.options = cms.untracked.PSet(",
       " ",
       ")",
       " ",
       "process.configurationMetadata = cms.untracked.PSet(",
       "    annotation = cms.untracked.string('Configuration/GenProduction/python/fragment.py nevts:{nevents}'),",
       "    name = cms.untracked.string('Applications'),",
       "    version = cms.untracked.string('$Revision: 1.19 $')",
       ")",
       " ",
       " ",
       "process.FEVTDEBUGoutput = cms.OutputModule('PoolOutputModule',",
       "    SelectEvents = cms.untracked.PSet(",
       "        SelectEvents = cms.vstring('generation_step')",
       "    ),",
       "    dataset = cms.untracked.PSet(",
       "        dataTier = cms.untracked.string('{step}'),",
       "        filterName = cms.untracked.string('')",
       "    ),",
       "    fileName = cms.untracked.string('file:step1.root'),",
       "    outputCommands = process.FEVTDEBUGEventContent.outputCommands,",
       "    splitLevel = cms.untracked.int32(0)",
       ")",
       " ",
       "process.genstepfilter.triggerConditions=cms.vstring('generation_step')",
       "from Configuration.AlCa.GlobalTag import GlobalTag",
       "process.GlobalTag = GlobalTag(process.GlobalTag, '{gt}', '')",
       " ",
       " ",
       ]
      
    part1 = '\n'.join(part1)
    part1 = part1.format(
          gt = self.driver.global_tag,
          bs = self.driver.beamspot_conditions,
          era = self.driver.era,
          nevents = self.nevents_perjob,
          step = 'GEN-SIM' if not self.dogenonly else 'GEN',
        )

    fragment = [
       "process.generator = cms.EDFilter('Pythia8GeneratorFilter',",
       "    ExternalDecays = cms.PSet(",
       "        EvtGen130 = cms.untracked.PSet(",
       "            convertPythiaCodes = cms.untracked.bool(False),",
       "            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),",
       "            list_forced_decays = cms.vstring(",
       "                'MyB_s0',",
       "                'Myanti-B_s0',", 
       "                'MyPhi'",
       "            ),",
       "            operates_on_particles = cms.vint32(531, -531, 333),",
       "            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2020.pdl'),",
       "            user_decay_embedded = cms.vstring(",
       "                'Define Hp 0.558',", 
       "                'Define Hz 0.616',", 
       "                'Define Hm 0.557',", 
       "                'Define pHp 2.469',", 
       "                'Define pHz 0.0',", 
       "                'Define pHm 2.75',", 
       "                '',", 
       "                'Alias      MyB_s0   B_s0',", 
       "                'Alias      Myanti-B_s0   anti-B_s0',", 
       "                'ChargeConj Myanti-B_s0   MyB_s0',", 
       "                'Alias      MyPhi phi',", 
       "                'ChargeConj MyPhi MyPhi',", 
       "                '',", 
       "                'Decay MyB_s0',", 
       "                '1.000       MyPhi      MyPhi    PVV_CPLH 0.02 1 Hp pHp Hz pHz Hm pHm;',",
       "                'Enddecay',", 
       "                'Decay Myanti-B_s0',",
       "                '1.000       MyPhi      MyPhi    PVV_CPLH 0.02 1 Hp pHp Hz pHz Hm pHm;',",
       "                'Enddecay',", 
       "                'Decay MyPhi',",
       "                '1.000      K+         K-       VSS;',",
       "                'Enddecay',", 
       "                'End'",
       "            )",
       "        ),",
       "        parameterSets = cms.vstring('EvtGen130')",
       "    ),",
       "    PythiaParameters = cms.PSet(",
       "        parameterSets = cms.vstring(",
       "            'pythia8CommonSettings',", 
       "            'pythia8CP5Settings',", 
       "            'processParameters'",
       "        ),",
       "        processParameters = cms.vstring(",
       "            'SoftQCD:nonDiffractive = on',", 
       "            'PTFilter:filter = on',", 
       "            'PTFilter:quarkToFilter = 5',", 
       "            'PTFilter:scaleToFilter = 5.0'",
       "        ),",
       "        pythia8CP5Settings = cms.vstring(",
       "            'Tune:pp 14',", 
       "            'Tune:ee 7',", 
       "            'MultipartonInteractions:ecmPow=0.03344',",
       "            'MultipartonInteractions:bProfile=2',", 
       "            'MultipartonInteractions:pT0Ref=1.41',", 
       "            'MultipartonInteractions:coreRadius=0.7634',", 
       "            'MultipartonInteractions:coreFraction=0.63',", 
       "            'ColourReconnection:range=5.176',", 
       "            'SigmaTotal:zeroAXB=off',", 
       "            'SpaceShower:alphaSorder=2',", 
       "            'SpaceShower:alphaSvalue=0.118',", 
       "            'SigmaProcess:alphaSvalue=0.118',", 
       "            'SigmaProcess:alphaSorder=2',", 
       "            'MultipartonInteractions:alphaSvalue=0.118',",
       "            'MultipartonInteractions:alphaSorder=2',", 
       "            'TimeShower:alphaSorder=2',", 
       "            'TimeShower:alphaSvalue=0.118',", 
       "            'SigmaTotal:mode = 0',", 
       "            'SigmaTotal:sigmaEl = 21.89',",
       "            'SigmaTotal:sigmaTot = 100.309',",
       "            'PDF:pSet=LHAPDF6:NNPDF31_nnlo_as_0118'",
       "        ),",
       "        pythia8CommonSettings = cms.vstring(",
       "            'Tune:preferLHAPDF = 2',", 
       "            'Main:timesAllowErrors = 10000',",
       "            'Check:epTolErr = 0.01',", 
       "            'Beams:setProductionScalesFromLHEF = off',",
       "            #'SLHA:keepSM = on',", 
       "            'SLHA:minMassSM = 1000.',", 
       "            'ParticleDecays:limitTau0 = on',", 
       "            'ParticleDecays:tau0Max = 10',", 
       "            'ParticleDecays:allowPhotonRadiation = on'",
       "        )",
       "    ),",
       "    comEnergy = cms.double({com}),",
       "    maxEventsToPrint = cms.untracked.int32(0),",
       "    pythiaHepMCVerbosity = cms.untracked.bool(False),",
       "    pythiaPylistVerbosity = cms.untracked.int32(0)",
       ")",
       " ",
       "process.BFilter = cms.EDFilter('PythiaFilter',",
       "    MaxEta = cms.untracked.double(10.0),",
       "    MinEta = cms.untracked.double(-10.0),",
       "    ParticleID = cms.untracked.int32(531)",
       ")",
       " ",
       " ",
       "process.BsDecayFilter = cms.EDFilter('PythiaDauVFilter',",
       "    DaughterIDs = cms.untracked.vint32(333, 333),",
       "    MaxEta = cms.untracked.vdouble(10.0, 10.0),",
       "    MinEta = cms.untracked.vdouble(-10.0, -10.0),",
       "    MinPt = cms.untracked.vdouble(0.0, 0.0),",
       "    NumberDaughters = cms.untracked.int32(2),",
       "    ParticleID = cms.untracked.int32(531),",
       "    verbose = cms.untracked.int32(1)",
       ")",
       " ",
       "process.TriggerMuonFilter = cms.EDFilter('PythiaFilter',",
       "    MaxEta = cms.untracked.double(1.55),",
       "    MinEta = cms.untracked.double(-1.55),",
       "    MinPt = cms.untracked.double(6.8),",
       "    ParticleID = cms.untracked.int32(13)",
       ")",
       " ",
       " ",
       "process.PhiDecayFilter = cms.EDFilter('MCMultiParticleFilter',",
       "    AcceptMore = cms.bool(True),",
       "    EtaMax = cms.vdouble(2.5),",
       "    MotherID = cms.untracked.vint32(333),",
       "    NumRequired = cms.int32(4),",
       "    ParticleID = cms.vint32(321),",
       "    PtMin = cms.vdouble(0.3),",
       "    Status = cms.vint32(0)",
       ")",
       " ",
       " ",
       "process.ProductionFilterSequence = cms.Sequence(process.generator+process.BFilter+process.BsDecayFilter+process.PhiDecayFilter+process.TriggerMuonFilter)",
       ]

    fragment = '\n'.join(fragment)
    fragment = fragment.format(
          com = '13000.0' if self.year == 2018 else '13600.0',
          gt = self.driver.global_tag,
          bs = self.driver.beamspot_conditions,
          era = self.driver.era,
          nevents = self.nevents_perjob,
          step = 'GEN-SIM' if not self.dogenonly else 'GEN',
        )

    part3 = [ 
       " ",
       " ",
       "process.generation_step = cms.Path(process.pgen)",
       "process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)",
       "process.endjob_step = cms.EndPath(process.endOfProcess)",
       "process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)",
       " ",
       "process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.FEVTDEBUGoutput_step)",
       "from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask",
       "associatePatAlgosToolsTask(process)",
       "for path in process.paths:",
       "	getattr(process,path).insert(0, process.ProductionFilterSequence)",
       " ",
       "from Configuration.DataProcessing.Utils import addMonitoring",
       " ",
       "process = addMonitoring(process)",
       " ",
       "from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randSvc = RandomNumberServiceHelper(process.RandomNumberGeneratorService); randSvc.populate()",
       "from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete",
       "process = customiseEarlyDelete(process)",
       " ",
       ]

    part3 = '\n'.join(part3)

    step1_file = open('./production/{pl}/step1.py'.format(pl=self.pl), 'w+')
    step1_file.write(part1)
    step1_file.write('\n')
    #with open('../fragments/{}'.format(self.fragment_name), 'r') as fragment:
    #  for line in fragment:
    #    step1_file.write(line)
    step1_file.write(fragment)
    step1_file.write(part3)
    step1_file.close()


  def submit(self):
    dirname = './production/{pl}/'.format(pl=self.pl)
    os.chdir(dirname)
    command_submit = 'crab submit -c crab_config.py > job_submit.log'
    os.system(command_submit)
    os.chdir('../../')


  def process(self):
    print('------------------------------------------')
    print('              CRAB Launcher               ')
    print('------------------------------------------')
    print(' ')
    print('Will submit production with:')
    print('           - prod label: {}'.format(self.pl)) 
    print('------------------------------------------')
    print(' ') 

    print(' -> Creating output directories')
    self.createOuputDir()

    print('\n -> Creating CRAB configuration files')
    self.createCRABConfig()

    print('\n -> Creating submitter files')
    self.createSubmitter()

    print('\n -> Creating cmsDrivers')
    #self.createDriver()
    self.createStep1()

    if self.dosubmit or self.doresubmit:
      print('\n --> Submitting...')
      self.submit()

    print('\nDone')



if __name__ == "__main__":
  
  opt = getOptions()

  launcher = CRABLauncher(options=opt)
  launcher.process()



