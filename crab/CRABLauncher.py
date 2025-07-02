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
    def __init__(self, era, global_tag, beamspot_conditions, fragment_name, CMSSW):
        self.era = era
        self.global_tag = global_tag
        self.beamspot_conditions = beamspot_conditions
        self.fragment_name = fragment_name
        self.CMSSW = CMSSW


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
    acc = 0.86 #0.97 

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
                beamspot_conditions = 'Realistic25ns13TeVEarly2018Collision',
                fragment_name = 'BsToPhiPhiTo4K_cfi.py',
                CMSSW = 'CMSSW_10_6_28',
                )
    elif self.year == 2022:
        self.driver = CMSDriver(
                era = 'Run3',
                global_tag = '124X_mcRun3_2022_realistic_v12',
                beamspot_conditions = 'Realistic25ns13p6TeVEarly2022Collision',
                fragment_name = 'BsToPhiPhiTo4K_Run3_cfi.py',
                CMSSW = 'CMSSW_12_4_24',
                )
    elif self.year == 2024:
        self.driver = CMSDriver(
                era = 'Run3_2024',
                global_tag = '140X_mcRun3_2024_realistic_v26',
                beamspot_conditions = 'Realistic',
                fragment_name = 'BsToPhiPhiTo4K_Run3_cfi.py',
                CMSSW = 'CMSSW_14_0_21',
                )


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
        cmssw = self.driver.CMSSW,
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
        #'echo "Source new environment"',
        #'source /cvmfs/cms.cern.ch/cmsset_default.sh',
        #'if [ -r {cmssw}/src ] ; then',
        #'  echo release {cmssw} already exists',
        #'else',
        #'  scram p CMSSW {cmssw}',
        #'fi',
        #'cd {cmssw}/src',
        #'eval `scram runtime -sh`',
        #'scram b',
        #'cd ../..',
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
          frgmt = self.driver.fragment_name,
          cmd = command,
          cmd_rpl = command_replace,
          cmssw = self.driver.CMSSW,
          )

      submitter_tmp_file = open('submitter_tmp.sh', 'w+')
      submitter_tmp_file.write(submitter_tmp)
      submitter_tmp_file.close()
      
      os.system('sh submitter_tmp.sh > log.txt')

      #os.system('rm submitter_tmp.sh') 
      os.system('rm log.txt') 


  def createStep1(self):

    step1_template_name = './cmsDrivers/step_GEN_{}.py'.format(self.year)
    step1_template = open(step1_template_name, 'r')
    lines = step1_template.readlines()

    step1_file = open('./production/{pl}/step1.py'.format(pl=self.pl), 'w+')

    for line in lines:
      if 'input = cms.untracked.int32(-1)' not in line:
        step1_file.write(line)
      else:
        step1_file.write('    input = cms.untracked.int32({}),\n'.format(self.nevents_perjob))

    step1_file.close()
    print('--> ./production/{pl}/step1.py created'.format(pl=self.pl))


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
    if not self.dogenonly:
      self.createDriver()
    else:
      self.createStep1()

    if self.dosubmit or self.doresubmit:
      print('\n --> Submitting...')
      self.submit()

    print('\nDone')



if __name__ == "__main__":
  
  opt = getOptions()

  launcher = CRABLauncher(options=opt)
  launcher.process()



