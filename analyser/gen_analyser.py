# http://home.thep.lu.se/~torbjorn/talks/fnal04lha.pdf

import os
import sys
import glob
import ROOT
import numpy as np
#import pandas as pd
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi
from scipy.constants import c as speed_of_light
from os import path


#def isAncestor(a, p):
#    if a == p :
#        return True
#    for i in xrange(0,p.numberOfMothers()):
#        if isAncestor(a,p.mother(i)):
#            return True
#    return False

#def weight_to_new_ctau(old_ctau, old_v2, new_v2, ct):
#    '''
#    Returns an event weight based on the ratio of the normalised lifetime distributions.
#    old_ctau: reference ctau
#    old_v2  : reference coupling squared
#    new_v2  : target coupling squared
#    ct      : heavy neutrino lifetime in the specific event
#    '''
#    new_ctau = old_ctau * old_v2 / new_v2
#    weight = old_ctau/new_ctau * np.exp( (1./old_ctau - 1./new_ctau) * ct )
#    return weight, new_ctau

def get_dxy(input):
  # taken from https://github.com/cms-sw/cmssw/blob/29f5fc15b34591745c5cd3c2c6eb9793aa6f371b/DataFormats/TrackReco/interface/TrackBase.h#L608
  # taking (0, 0, 0) as reference
  return (-input.vx() * input.py() + input.vy() * input.px()) / input.pt()


def get_list_daughters(mother, pdgid):
    list_daughters = [mother.daughter(jj) for jj in range(mother.numberOfDaughters()) if abs(mother.daughter(jj).pdgId()) == pdgid]
    sorted_list_daughters = sorted([ii for ii in list_daughters], key = lambda x : x.pt(), reverse=True)

    return sorted_list_daughters
    #return [mother.daughter(jj) for jj in range(mother.numberOfDaughters()) if abs(mother.daughter(jj).pdgId()) == pdgid]


def get_decay_chain(part, pdgid_list=[]):
  ancestors = sorted([part.mother(jj) for jj in range(part.numberOfMothers())], key = lambda x : x.pt(), reverse=True)
  if len(ancestors) != 0:
    ancestor = ancestors[0]
    ancestor_pdgid = ancestor.pdgId()
    if abs(ancestor_pdgid) in [1, 2, 3, 4, 5, 21, 1103, 2101, 2103, 2203, 3101, 3103, 3201, 3203, 3303, 4101, 4103, 4201, 4203, 4301, 4303, 4403, 5101, 5201, 5203, 5301, 5303, 5401, 5403, 5503]:
      return pdgid_list
    else:
      pdgid_list.append(ancestor_pdgid)
      return get_decay_chain(ancestor, pdgid_list)
  else:
      return []


branches = [
    'run',  
    'lumi', 
    'event',
     
    'Bs_pt',
    'Bs_eta',
    'Bs_phi',
    'Bs_mass',
    'Bs_q',
    'Bs_pdgid',

    'phi1_pt',
    'phi1_eta',
    'phi1_phi',
    'phi1_mass',
    'phi1_q',
    'phi1_pdgid',

    'phi2_pt',
    'phi2_eta',
    'phi2_phi',
    'phi2_mass',
    'phi2_q',
    'phi2_pdgid',

    'k1_pt',
    'k1_eta',
    'k1_phi',
    'k1_mass',
    'k1_q',
    'k1_pdgid',

    'k2_pt',
    'k2_eta',
    'k2_phi',
    'k2_mass',
    'k2_q',
    'k2_pdgid',

    'k3_pt',
    'k3_eta',
    'k3_phi',
    'k3_mass',
    'k3_q',
    'k3_pdgid',

    'k4_pt',
    'k4_eta',
    'k4_phi',
    'k4_mass',
    'k4_q',
    'k4_pdgid',

    'invmass_phi1phi2',
    'invmass_k1k2k3k4',
    'invmass_k1k2',
    'invmass_k3k4',

    'charge_k1k2',
    'charge_k3k4',
    'charge_phi1phi2',

    # unmixed Bs decays
    'chain_Bs_to_phiphi',
    'chain_Bsstar_to_Bs_to_phiphi',
    'chain_Bc_to_Bs_to_phiphi',
    'chain_Bcstar_to_Bc_to_Bs_to_phiphi',
    'chain_Bc_to_Bsstar_to_Bs_to_phiphi',
    'chain_Bcstar_to_Bc_to_Bsstar_to_Bs_to_phiphi',

    # mixed Bs decays
    'chain_Bsbar_to_Bs_to_phiphi',
    'chain_Bsstarbar_to_Bsbar_to_Bs_to_phiphi',
    'chain_Bcbar_to_Bsbar_to_Bs_to_phiphi',
    'chain_Bcbar_to_Bsstarbar_to_Bsbar_to_Bs_to_phiphi',
    'chain_Bcstarbar_to_Bcbar_to_Bsbar_to_Bs_to_phiphi',
    'chain_Bcstarbar_to_Bcbar_to_Bsstarbar_to_Bsbar_to_Bs_to_phiphi',

    # unmixed Bsbar decays
    'chain_Bsbar_to_phiphi',
    'chain_Bsstarbar_to_Bsbar_to_phiphi',
    'chain_Bcbar_to_Bsbar_to_phiphi',
    'chain_Bcstarbar_to_Bcbar_to_Bsbar_to_phiphi',
    'chain_Bcbar_to_Bsstarbar_to_Bsbar_to_phiphi',
    'chain_Bcstarbar_to_Bcbar_to_Bsstarbar_to_Bsbar_to_phiphi',

    # mixed Bsbar decays
    'chain_Bs_to_Bsbar_to_phiphi',
    'chain_Bsstar_to_Bs_to_Bsbar_to_phiphi',
    'chain_Bc_to_Bs_to_Bsbar_to_phiphi',
    'chain_Bc_to_Bsstar_to_Bs_to_Bsbar_to_phiphi',
    'chain_Bcstar_to_Bc_to_Bs_to_Bsbar_to_phiphi',
    'chain_Bcstar_to_Bc_to_Bsstar_to_Bs_to_Bsbar_to_phiphi',


    #TODO add deltaR

]

## couplings to be tested, for which the reweight is run
#new_v2s = [
#    1e-10, 
#    5e-10, 
#    1e-9, 
#    5e-9, 
#    1e-8, 
#    5e-8, 
#    1e-7, 
#    5e-7, 
#    1e-6, 
#    5e-6, 
#    6e-06, 
#    8e-06, 
#    1e-5, 
#    2e-5, 
#    3e-5, 
#    4e-5, 
#    5e-5, 
#    7e-05, 
#    0.0001, 
#    0.0002, 
#    0.00025, 
#    0.0003, 
#    0.0005, 
#    0.0012,
#]
#
## ctau weights, each coupling gets its own weight
#weights = OrderedDict(zip(new_v2s, np.ones(len(new_v2s))))

# add weight branches
#for vv in new_v2s:
#    branches.append('weight_%s'      %(str(vv).replace('-', 'm')))
#    branches.append('ctau_%s'        %(str(vv).replace('-', 'm')))
#    branches.append('xs_scale_to_%s' %(str(vv).replace('-', 'm')))


#handles = OrderedDict()
#handles['genP'] = ('genParticles' , Handle('std::vector<reco::GenParticle>'))
#
## output file and tree gymnastics
##outfile = ROOT.TFile.Open('outputfiles/test_modfilter_v4_m3_ctau1000.root', 'recreate')
#outfile = ROOT.TFile.Open('outputfiles/test.root', 'recreate')
#ntuple  = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
#tofill = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       
#
## get to the real thing
#print 'loading the file ...'
#infiles = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/BHNLsGen/test_modfilter_v3_n10000000_njt500/mass4.5_ctau1.0/step1*.root'
##infiles = '/pnfs/psi.ch/cms/trivcat/store/user/mratti/BHNLsGen/V33_stats_Lxy1300_tkPt500MeV_lepPt400MeV/mass4.5_ctau1.0/step1*root'


#def treeProducer(infiles, outdir, outfilename, lepton):
def treeProducer(infiles, outdir, outfilename):
  # get the files
  print 'loading the file ...'
  files = glob.glob(infiles)
  events = Events(files)
  print 'nevts: ',events.size()
  print '... done!'

  # get the handles
  handles = OrderedDict()
  handles['genP'] = ('genParticles' , Handle('std::vector<reco::GenParticle>'))

  # output file and tree gymnastics
  outfile = ROOT.TFile.Open('{}/{}'.format(outdir, outfilename), 'recreate')
  ntuple  = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
  tofill = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99       

  # get reconstruction weights
  #bins_pt  = [(0,7), (7,10), (10,15),(15,30), (30,1000)] #rows
  #bins_lxy = [(0,10),(10,30),(30,50),(50,100),(100,150),(150,300),(300,500),(500,1e10)] #columns, everything in mm
  #reco_weights = pd.read_csv('reco_weights_updated.csv', sep=',', comment='#')

  # file ct weight
  #file_ct = ROOT.TFile.Open('weight_ct_m3_ctau1.root', 'READ')
  #hist_weight_ct = file_ct.Get('hist_weight')

  count_tot = 0
  count_mu_trg = 0
  count_B_meson = 0
  count_B_baryon = 0
  count_other = 0
  list_mu_mother_pdgid = []

  count_acceptance = 0

  do_trgmu_study = False

  for i, event in enumerate(events):
    #if float(i)>1000: continue
    #print '\n\n Event {}'.format(i)

    # access the handles
    for k, v in handles.iteritems():
        event.getByLabel(v[0], v[1])
        setattr(event, k, v[1].product())

    if i%1000==0:
        percentage = float(i)/events.size()*100.
        print '\t===> processing %d / %d event \t completed %.1f%s' %(i, events.size(), percentage, '%')

    count_tot += 1

    # searching for trigger muon
    if do_trgmu_study:
      muons = [ip for ip in event.genP if (abs(ip.pdgId())==13 and ip.pt()>6.8 and abs(ip.eta())<1.55)]
      the_muons = sorted([ii for ii in muons], key = lambda x : x.pt(), reverse=True)
      if len(the_muons) != 0:
        count_mu_trg += 1
        #print 'number of triggering muon: {}'.format(len(the_muons))
        for the_muon in the_muons:
         mu_mothers = [the_muon.mother(jj) for jj in range(the_muon.numberOfMothers())]
         the_mu_mothers = sorted([ii for ii in mu_mothers], key = lambda x : x.pt(), reverse=True)
         #print 'the mothers are:'
         for the_mu_mother in the_mu_mothers:
           #print the_mu_mother.pdgId()
           if abs(the_mu_mother.pdgId()) in [511, 521, 531]:
             count_B_meson += 1
           elif abs(the_mu_mother.pdgId()) in [5122, 5132, 5232]:
             count_B_baryon += 1
           else:
             count_other += 1 
             print 'other mother pdgid: {}'.format(abs(the_mu_mother.pdgId()))

           if abs(the_mu_mother.pdgId()) not in list_mu_mother_pdgid:
             list_mu_mother_pdgid.append(abs(the_mu_mother.pdgId()))

           #mu_grandmothers = [the_mu_mother.mother(jj) for jj in range(the_mu_mother.numberOfMothers())]
           #the_mu_grandmothers = sorted([ii for ii in mu_grandmothers], key = lambda x : x.pt(), reverse=True)
           #for the_mu_grandmother in the_mu_grandmothers:
           #  print 'grandmother:'
           #  print the_mu_grandmother.pdgId()

      #else:
      #  continue #TODO remove eventually

    # list of mother pdgid:
    # [521, -521, -431, -511, 511, -421, -5122, -411, 411, -531, 5122, 443, 13, 15, 421, 531, -5132, 5232, 431]
    # 511: B0
    # 521: B+
    # 531: Bs
    # 421: D0
    # 431: Ds
    # 443: J/psi
    # 5122: Lambda b
    # 5132: chsi b
    # 5232: chsi b 0
    # 13: mu
    # 15: tau

    # get the Bs mesons
    the_Bs_mesons = [ip for ip in event.genP if abs(ip.pdgId())==531 and ip.isLastCopy()]
    #print 'number of Bs mesons: {}'.format(len(the_Bs_mesons))

    # find Bs meson that has two phi daughters
    idx_Bs = -99
    count_Bs_has_two_phi_daughters = 0
    for iBs, the_Bs_meson in enumerate(the_Bs_mesons):
      #Bs_daughters = [the_Bs_meson.daughter(jj) for jj in range(the_Bs_meson.numberOfDaughters()) if the_Bs_meson.daughter(jj).pdgId() == 333]
      number_phi_daughters = len(get_list_daughters(mother=the_Bs_meson, pdgid=333))
      if number_phi_daughters == 2:
        #print 'B meson {} has 2 phi daughters'.format(iBs)
        idx_Bs = iBs
        count_Bs_has_two_phi_daughters += 1
      #elif number_phi_daughters > 0 and number_phi_daughters != 2:
      #  print 'WARNING: there are not exactly 2 phi daughters'
     # elif number_phi_daughters > 0 and number_phi_daughters < 2:
     #   print 'WARNING: there are less than 2 phi daughtersi for Bs meson {}'.format(iBs)
     #   print 'number of Bs mesons: {}'.format(len(the_Bs_mesons))

    if count_Bs_has_two_phi_daughters == 0:
      print 'WARNING - no Bs meson decaying to two phi mesons was found'
      print '--> skipping event'
      continue
    if count_Bs_has_two_phi_daughters != 1:
      print 'number of Bs mesons with 2 phi daughters: {}'.format(count_Bs_has_two_phi_daughters)
        
    event.the_Bs = the_Bs_mesons[idx_Bs]

    # searching for mother of Bs meson
    event.the_Bs.mothers = [event.the_Bs.mother(jj) for jj in range(event.the_Bs.numberOfMothers())]
    the_Bs_mothers = sorted([ii for ii in event.the_Bs.mothers], key = lambda x : x.pt(), reverse=True)
    event.the_Bs_mother = the_Bs_mothers[0]

    decay_chain = get_decay_chain(event.the_Bs, [])

    # unmixed Bs decays
    chain_Bs_to_phiphi = 0
    chain_Bsstar_to_Bs_to_phiphi = 0
    chain_Bc_to_Bs_to_phiphi = 0
    chain_Bcstar_to_Bc_to_Bs_to_phiphi = 0
    chain_Bc_to_Bsstar_to_Bs_to_phiphi = 0
    chain_Bcstar_to_Bc_to_Bsstar_to_Bs_to_phiphi = 0

    # mixed Bs decays
    chain_Bsbar_to_Bs_to_phiphi = 0
    chain_Bsstarbar_to_Bsbar_to_Bs_to_phiphi = 0
    chain_Bcbar_to_Bsbar_to_Bs_to_phiphi = 0
    chain_Bcbar_to_Bsstarbar_to_Bsbar_to_Bs_to_phiphi = 0
    chain_Bcstarbar_to_Bcbar_to_Bsbar_to_Bs_to_phiphi = 0
    chain_Bcstarbar_to_Bcbar_to_Bsstarbar_to_Bsbar_to_Bs_to_phiphi = 0

    # unmixed Bsbar decays
    chain_Bsbar_to_phiphi = 0
    chain_Bsstarbar_to_Bsbar_to_phiphi = 0
    chain_Bcbar_to_Bsbar_to_phiphi = 0
    chain_Bcstarbar_to_Bcbar_to_Bsbar_to_phiphi = 0
    chain_Bcbar_to_Bsstarbar_to_Bsbar_to_phiphi = 0
    chain_Bcstarbar_to_Bcbar_to_Bsstarbar_to_Bsbar_to_phiphi = 0

    # mixed Bsbar decays
    chain_Bs_to_Bsbar_to_phiphi = 0
    chain_Bsstar_to_Bs_to_Bsbar_to_phiphi = 0
    chain_Bc_to_Bs_to_Bsbar_to_phiphi = 0
    chain_Bc_to_Bsstar_to_Bs_to_Bsbar_to_phiphi = 0
    chain_Bcstar_to_Bc_to_Bs_to_Bsbar_to_phiphi = 0
    chain_Bcstar_to_Bc_to_Bsstar_to_Bs_to_Bsbar_to_phiphi = 0

    # Bs decays
    if event.the_Bs.pdgId() == 531:
      if decay_chain == []:
        chain_Bs_to_phiphi = 1
      elif decay_chain == [533]:
        chain_Bsstar_to_Bs_to_phiphi = 1
      elif decay_chain == [541]:
        chain_Bc_to_Bs_to_phiphi = 1
      elif decay_chain == [541, 543]:
        chain_Bcstar_to_Bc_to_Bs_to_phiphi = 1
      elif decay_chain == [533, 541]:
        chain_Bc_to_Bsstar_to_Bs_to_phiphi = 1
      elif decay_chain == [533, 541, 543]:
        chain_Bcstar_to_Bc_to_Bsstar_to_Bs_to_phiphi = 1
      elif decay_chain == [-531]:
        chain_Bsbar_to_Bs_to_phiphi = 1
      elif decay_chain == [-531, -533]:
        chain_Bsstarbar_to_Bsbar_to_Bs_to_phiphi = 1
      elif decay_chain == [-531, -541]:
        chain_Bcbar_to_Bsbar_to_Bs_to_phiphi = 1
      elif decay_chain == [-531, -533, -541]:
        chain_Bcbar_to_Bsstarbar_to_Bsbar_to_Bs_to_phiphi = 1
      elif decay_chain == [-531, -541, -543]:
        chain_Bcstarbar_to_Bcbar_to_Bsbar_to_Bs_to_phiphi = 1
      elif decay_chain == [-531, -533, -541, -543]:
        chain_Bcstarbar_to_Bcbar_to_Bsstarbar_to_Bsbar_to_Bs_to_phiphi = 1
      else:
        print "WARNING: unknown decay chain for Bs: {}".format(decay_chain)

    # Bsbar decays
    if event.the_Bs.pdgId() == -531:
      if decay_chain == []:
        chain_Bsbar_to_phiphi = 1
      elif decay_chain == [-533]:
        chain_Bsstarbar_to_Bsbar_to_phiphi = 1
      elif decay_chain == [-541]:
        chain_Bcbar_to_Bsbar_to_phiphi = 1
      elif decay_chain == [-541, -543]:
        chain_Bcstarbar_to_Bcbar_to_Bsbar_to_phiphi = 1
      elif decay_chain == [-533, -541]:
        chain_Bcbar_to_Bsstarbar_to_Bsbar_to_phiphi = 1
      elif decay_chain == [-533, -541, -543]:
        chain_Bcstarbar_to_Bcbar_to_Bsstarbar_to_Bsbar_to_phiphi = 1
      elif decay_chain == [531]:
        chain_Bs_to_Bsbar_to_phiphi = 1
      elif decay_chain == [531, 533]:
        chain_Bsstar_to_Bs_to_Bsbar_to_phiphi = 1
      elif decay_chain == [531, 541]:
        chain_Bc_to_Bs_to_Bsbar_to_phiphi = 1
      elif decay_chain == [531, 533, 541]:
        chain_Bc_to_Bsstar_to_Bs_to_Bsbar_to_phiphi = 1
      elif decay_chain == [531, 541, 543]:
        chain_Bcstar_to_Bc_to_Bs_to_Bsbar_to_phiphi = 1
      elif decay_chain == [531, 533, 541, 543]:
        chain_Bcstar_to_Bc_to_Bsstar_to_Bs_to_Bsbar_to_phiphi = 1
      else:
        print "WARNING: unknown decay chain for Bsbar: {}".format(decay_chain)


    # get Bs daughters
    event.the_phi1 = get_list_daughters(mother=event.the_Bs, pdgid=333)[0]
    event.the_phi2 = get_list_daughters(mother=event.the_Bs, pdgid=333)[1]

    event.the_phis_p4 = event.the_phi1.polarP4() + event.the_phi2.polarP4()
    #Bs_invmass = event.the_phis_p4.mass()
    #Bs_invpt = event.the_phis_p4.pt()
    #print '{}%'.format(((Bs_invmass - 5.367) / 5.367) * 100)
    #print '{}%'.format(((event.the_Bs.p4().mass() - 5.367) / 5.367) * 100)
    #print Bs_invpt

    # daughters of phi1 (kaons)
    phi1_daughters = get_list_daughters(mother=event.the_phi1, pdgid=321)
    #print len(phi1_daughters)
    #for phi1_daughter in phi1_daughters:
    #  print phi1_daughter.pdgId()
    if len(phi1_daughters) == 0:
      print 'WARNING: no kaon found as daughter of phi1'
      print '--> skipping'
      continue
    event.the_k1 = phi1_daughters[0]  
    event.the_k2 = phi1_daughters[1]  
        
    # daughters of phi2 (kaons)
    phi2_daughters = get_list_daughters(mother=event.the_phi2, pdgid=321)
    if len(phi2_daughters) == 0:
      print 'WARNING: no kaon found as daughter of phi2'
      print '--> skipping'
      continue
    event.the_k3 = phi2_daughters[0]  
    event.the_k4 = phi2_daughters[1]  

    # study acceptance
    # with trigger muon
    #if(event.the_k1.pt() > 0.5 and abs(event.the_k1.eta()) < 2.5 and event.the_k2.pt() > 0.5 and abs(event.the_k2.eta()) < 2.5 and event.the_k3.pt() > 0.5 and abs(event.the_k3.eta()) < 2.5 and event.the_k4.pt() > 0.5 and abs(event.the_k4.eta()) < 2.5 and the_muons[0].pt() > 7 and abs(the_muons[0].eta()) < 1.5):

    # without trigger muon
    if(event.the_k1.pt() > 0.5 and abs(event.the_k1.eta()) < 2.5 and event.the_k2.pt() > 0.5 and abs(event.the_k2.eta()) < 2.5 and event.the_k3.pt() > 0.5 and abs(event.the_k3.eta()) < 2.5 and event.the_k4.pt() > 0.5 and abs(event.the_k4.eta()) < 2.5):

    #if((event.the_k1.pt() > 0.5 and abs(event.the_k1.eta()) < 2.5 and event.the_k2.pt() > 0.5 and abs(event.the_k2.eta()) < 2.5) or (event.the_k3.pt() > 0.5 and abs(event.the_k3.eta()) < 2.5 and event.the_k4.pt() > 0.5 and abs(event.the_k4.eta()) < 2.5)):
    #if(event.the_k1.pt() > 0.5 and abs(event.the_k1.eta()) < 2.5 and event.the_k2.pt() > 0.5 and abs(event.the_k2.eta()) < 2.5 and event.the_k3.pt() > 0.5 and abs(event.the_k3.eta()) < 2.5 and event.the_k4.pt() > 0.5 and abs(event.the_k4.eta()) < 2.5):
    #if(event.the_k1.pt() > 0.5 and abs(event.the_k1.eta()) < 2.5 and event.the_k2.pt() > 0.5 and abs(event.the_k2.eta()) < 2.5):
    #if(event.the_k3.pt() > 0.5 and abs(event.the_k3.eta()) < 2.5 and event.the_k4.pt() > 0.5 and abs(event.the_k4.eta()) < 2.5):
    #if(event.the_k1.pt() > 0.7 and abs(event.the_k1.eta()) < 2.4 and event.the_k2.pt() > 0.7 and abs(event.the_k2.eta()) < 2.4 and event.the_k3.pt() > 0.7 and abs(event.the_k3.eta()) < 2.4 and event.the_k4.pt() > 0.7 and abs(event.the_k4.eta()) < 2.4):
      count_acceptance += 1


    #TODO search for triggering muon on the other side?
    #  if event.the_hn.lep.pt() > 6.8 and abs(event.the_hn.lep.eta()) < 1.55: event.the_hn.lep.satisfies_BParkHLT_cond = 1
    # checking that there is at least one triggering muon per event
    #event.atleast_1triggeringmuon_otherside = -1
    #if event.number_triggering_muon_hnlside == 0:
    #  the_trigger_muon = [ip for ip in event.genP if abs(ip.pdgId())==13 and ip.pt()>6.8 and abs(ip.eta())<1.55 and ip.isLastCopy()] 
    #  if len(the_trigger_muon): event.atleast_1triggeringmuon_otherside = 1
    #  else: event.atleast_1triggeringmuon_otherside = 0


    #TODO get displacement info?
    ## identify the primary vertex
    ## for that, needs the B muon
    #if len(the_pls):
    #  event.the_hn.the_pv = event.the_pl.vertex()
    #
    ## identify the secondary vertex
    #if len(the_lep_daughters):
    #  event.the_hn.the_sv = event.the_hn.lep.vertex()
    #
    ## 2D transverse and 3D displacement, Pythagoras
    #if len(the_pls) and len(the_lep_daughters):
    #  event.Lxy  = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
    #                       (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2) * 10 # in mm

    #  event.Lxyz = np.sqrt((event.the_hn.the_pv.x() - event.the_hn.the_sv.x())**2 + \
    #                       (event.the_hn.the_pv.y() - event.the_hn.the_sv.y())**2 + \
    #                       (event.the_hn.the_pv.z() - event.the_hn.the_sv.z())**2) * 10 # in mm

    ## per event ct, as derived from the flight distance and Lorentz boost  
    #event.the_hn.beta  = event.the_hn.p4().Beta()
    #event.the_hn.gamma = event.the_hn.p4().Gamma()
    #
    ## we get the lifetime from the kinematics 
    #if len(the_pls) and len(the_lep_daughters):
    #  event.the_hn.ct_reco = event.Lxyz / (event.the_hn.beta * event.the_hn.gamma)
    #  #print 'hnl ct reco: {a}'.format(a=event.the_hn.ct_reco)


    #event.Lxyz_pl = np.sqrt((event.the_b_mother.vx() - event.the_pl.vx())**2 + \
    #                    (event.the_b_mother.vy() - event.the_pl.vy())**2 + \
    #                    (event.the_b_mother.vz() - event.the_pl.vz())**2) * 10 # in mm
    #
    ## get the lifetime of the B
    #event.the_b_mother.beta  = event.the_b_mother.p4().Beta()
    #event.the_b_mother.gamma = event.the_b_mother.p4().Gamma()
   
    #event.the_b_mother.ct_reco = event.Lxyz_pl / (event.the_b_mother.beta * event.the_b_mother.gamma)

    # reset before filling
    for k, v in tofill.iteritems(): tofill[k] = -99. # initialise before filling

    tofill['run'] = event.eventAuxiliary().run()
    tofill['lumi'] = event.eventAuxiliary().luminosityBlock()
    tofill['event'] = event.eventAuxiliary().event()
     
    tofill['Bs_pt'   ] = event.the_Bs.pt()     
    tofill['Bs_eta'  ] = event.the_Bs.eta()    
    tofill['Bs_phi'  ] = event.the_Bs.phi()    
    tofill['Bs_mass' ] = event.the_Bs.mass()   
    tofill['Bs_q'    ] = event.the_Bs.charge()   
    tofill['Bs_pdgid'] = event.the_Bs.pdgId()   
    #tofill['Bs_ct_reco'] = event.the_Bs.ct_reco   
     
    tofill['phi1_pt'     ] = event.the_phi1.pt()     
    tofill['phi1_eta'    ] = event.the_phi1.eta()    
    tofill['phi1_phi'    ] = event.the_phi1.phi()    
    tofill['phi1_mass'   ] = event.the_phi1.mass()   
    tofill['phi1_q'      ] = event.the_phi1.charge()
    tofill['phi1_pdgid'  ] = event.the_phi1.pdgId()   

    tofill['phi2_pt'     ] = event.the_phi2.pt()     
    tofill['phi2_eta'    ] = event.the_phi2.eta()    
    tofill['phi2_phi'    ] = event.the_phi2.phi()    
    tofill['phi2_mass'   ] = event.the_phi2.mass()   
    tofill['phi2_q'      ] = event.the_phi2.charge()
    tofill['phi2_pdgid'  ] = event.the_phi2.pdgId()   

    tofill['k1_pt'     ] = event.the_k1.pt()     
    tofill['k1_eta'    ] = event.the_k1.eta()    
    tofill['k1_phi'    ] = event.the_k1.phi()    
    tofill['k1_mass'   ] = event.the_k1.mass()   
    tofill['k1_q'      ] = event.the_k1.charge()
    tofill['k1_pdgid'  ] = event.the_k1.pdgId()   

    tofill['k2_pt'     ] = event.the_k2.pt()     
    tofill['k2_eta'    ] = event.the_k2.eta()    
    tofill['k2_phi'    ] = event.the_k2.phi()    
    tofill['k2_mass'   ] = event.the_k2.mass()   
    tofill['k2_q'      ] = event.the_k2.charge()
    tofill['k2_pdgid'  ] = event.the_k2.pdgId()   

    tofill['k3_pt'     ] = event.the_k3.pt()     
    tofill['k3_eta'    ] = event.the_k3.eta()    
    tofill['k3_phi'    ] = event.the_k3.phi()    
    tofill['k3_mass'   ] = event.the_k3.mass()   
    tofill['k3_q'      ] = event.the_k3.charge()
    tofill['k3_pdgid'  ] = event.the_k3.pdgId()   

    tofill['k4_pt'     ] = event.the_k4.pt()     
    tofill['k4_eta'    ] = event.the_k4.eta()    
    tofill['k4_phi'    ] = event.the_k4.phi()    
    tofill['k4_mass'   ] = event.the_k4.mass()   
    tofill['k4_q'      ] = event.the_k4.charge()
    tofill['k4_pdgid'  ] = event.the_k4.pdgId()   

    # invariant masses   
    tofill['invmass_phi1phi2'] = (event.the_phi1.p4() + event.the_phi2.p4()).mass()
    tofill['invmass_k1k2k3k4'] = (event.the_k1.p4() + event.the_k2.p4() + event.the_k3.p4() + event.the_k4.p4()).mass()
    tofill['invmass_k1k2'] = (event.the_k1.p4() + event.the_k2.p4()).mass()
    tofill['invmass_k3k4'] = (event.the_k3.p4() + event.the_k4.p4()).mass()

    tofill['charge_k1k2'] = event.the_k1.charge() + event.the_k2.charge()
    tofill['charge_k3k4'] = event.the_k3.charge() + event.the_k4.charge()
    tofill['charge_phi1phi2'] = event.the_phi1.charge() + event.the_phi2.charge()

    # decay chains
    # unmixed Bs decays
    tofill['chain_Bs_to_phiphi'] = chain_Bs_to_phiphi
    tofill['chain_Bsstar_to_Bs_to_phiphi'] = chain_Bsstar_to_Bs_to_phiphi
    tofill['chain_Bc_to_Bs_to_phiphi'] = chain_Bc_to_Bs_to_phiphi
    tofill['chain_Bcstar_to_Bc_to_Bs_to_phiphi'] = chain_Bcstar_to_Bc_to_Bs_to_phiphi
    tofill['chain_Bc_to_Bsstar_to_Bs_to_phiphi'] = chain_Bc_to_Bsstar_to_Bs_to_phiphi
    tofill['chain_Bcstar_to_Bc_to_Bsstar_to_Bs_to_phiphi'] = chain_Bcstar_to_Bc_to_Bsstar_to_Bs_to_phiphi

    # mixed Bs decays
    tofill['chain_Bsbar_to_Bs_to_phiphi'] = chain_Bsbar_to_Bs_to_phiphi
    tofill['chain_Bsstarbar_to_Bsbar_to_Bs_to_phiphi'] = chain_Bsstarbar_to_Bsbar_to_Bs_to_phiphi
    tofill['chain_Bcbar_to_Bsbar_to_Bs_to_phiphi'] = chain_Bcbar_to_Bsbar_to_Bs_to_phiphi
    tofill['chain_Bcbar_to_Bsstarbar_to_Bsbar_to_Bs_to_phiphi'] = chain_Bcbar_to_Bsstarbar_to_Bsbar_to_Bs_to_phiphi
    tofill['chain_Bcstarbar_to_Bcbar_to_Bsbar_to_Bs_to_phiphi'] = chain_Bcstarbar_to_Bcbar_to_Bsbar_to_Bs_to_phiphi
    tofill['chain_Bcstarbar_to_Bcbar_to_Bsstarbar_to_Bsbar_to_Bs_to_phiphi'] = chain_Bcstarbar_to_Bcbar_to_Bsstarbar_to_Bsbar_to_Bs_to_phiphi

    # unmixed Bsbar decays
    tofill['chain_Bsbar_to_phiphi'] = chain_Bsbar_to_phiphi
    tofill['chain_Bsstarbar_to_Bsbar_to_phiphi'] = chain_Bsstarbar_to_Bsbar_to_phiphi
    tofill['chain_Bcbar_to_Bsbar_to_phiphi'] = chain_Bcbar_to_Bsbar_to_phiphi
    tofill['chain_Bcstarbar_to_Bcbar_to_Bsbar_to_phiphi'] = chain_Bcstarbar_to_Bcbar_to_Bsbar_to_phiphi
    tofill['chain_Bcbar_to_Bsstarbar_to_Bsbar_to_phiphi'] = chain_Bcbar_to_Bsstarbar_to_Bsbar_to_phiphi
    tofill['chain_Bcstarbar_to_Bcbar_to_Bsstarbar_to_Bsbar_to_phiphi'] = chain_Bcstarbar_to_Bcbar_to_Bsstarbar_to_Bsbar_to_phiphi

    # mixed Bsbar decays
    tofill['chain_Bs_to_Bsbar_to_phiphi'] = chain_Bs_to_Bsbar_to_phiphi
    tofill['chain_Bsstar_to_Bs_to_Bsbar_to_phiphi'] = chain_Bsstar_to_Bs_to_Bsbar_to_phiphi
    tofill['chain_Bc_to_Bs_to_Bsbar_to_phiphi'] = chain_Bc_to_Bs_to_Bsbar_to_phiphi
    tofill['chain_Bc_to_Bsstar_to_Bs_to_Bsbar_to_phiphi'] = chain_Bc_to_Bsstar_to_Bs_to_Bsbar_to_phiphi
    tofill['chain_Bcstar_to_Bc_to_Bs_to_Bsbar_to_phiphi'] = chain_Bcstar_to_Bc_to_Bs_to_Bsbar_to_phiphi
    tofill['chain_Bcstar_to_Bc_to_Bsstar_to_Bs_to_Bsbar_to_phiphi'] = chain_Bcstar_to_Bc_to_Bsstar_to_Bs_to_Bsbar_to_phiphi

    #TODO add dxy?   
    #tofill['mu_fromB_dxy'     ] = get_dxy(event.the_pl)
    
    ntuple.Fill(array('f',tofill.values()))

  if do_trgmu_study:
    print 'list of mother pdgid:'
    print list_mu_mother_pdgid

    print 'percentage of events with triggering muon: {}%'.format((count_mu_trg / float(count_tot)) * 100)
    print 'percentage of events with triggering muon from B meson: {}%'.format((count_B_meson / float(count_mu_trg)) * 100)
    print 'percentage of events with triggering muon from B baryon: {}%'.format((count_B_baryon / float(count_mu_trg)) * 100)
    print 'percentage of events with triggering muon from other particle: {}%'.format((count_other / float(count_mu_trg)) * 100)

  acceptance = count_acceptance / float(count_tot) * 100
  print 'acceptance = {} / {} = {}%'.format(count_acceptance, count_tot, acceptance)

  outfile.cd()
  ntuple.Write()
  outfile.Close()


if __name__ == "__main__":
  #version_label = 'test_v0'
  #version_label = 'crab_102X_crab_v0'
  version_label = '102X_crab_trgmu_filter'
  #version_label = 'test_trgmu_v0'
  #user = 'anlyon'
  #lepton = 'all'

  #if lepton not in ['muon', 'electron', 'all']:
  #  raise RuntimeError("Lepton not known. Choose among ['muon', 'electron', 'all']")

  #indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/{}/BHNLsGen/{}'.format(user, version_label)
  #indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/{}/CPVGen/{}'.format(user, version_label)
  #indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/102X_crab_v0/BsToPhiPhiTo4K/crab_102X_crab_v0_BsToPhiPhiTo4K_20250130_141241/250130_131709/0000'

  #indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/102X_crab_v0/BsToPhiPhiTo4K/crab_102X_crab_v0_BsToPhiPhiTo4K_20250130_141241/250130_131709/0000'
  #indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/102X_crab_v0/BsToPhiPhiTo4K/crab_102X_crab_v0_BsToPhiPhiTo4K_20250203_151733/250203_141746/0000'
  #indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/test_trgmu_v0'
 
  #indirectory = '/pnfs/psi.ch/cms/trivcat/store/user/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/crab_102X_crab_trgmu_filter_BsToPhiPhiTo4K_20250212_225911/250212_215924/*'


  indirectory = '/eos/user/a/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/crab_102X_crab_trgmu_filter_BsToPhiPhiTo4K_20250212_225911/250212_215924/000*/'

  # get all the subdirectories (signal points)
  #pointdirs = [f for f in glob.glob('{}/*'.format(indirectory))]

  #for pointdir in pointdirs:
  print '----------------------'
  print ' Analysing {}'.format(indirectory)
  print '----------------------'

  #if 'mass2.0_ctau10.0' not in pointdir: continue

  #pointname = pointdir[pointdir.rfind('/')+1:].replace('.', 'p')

  infiles = '{}/step1*root'.format(indirectory)
  outdir = './outputfiles/{}'.format(version_label)
  outfilename = 'genTree.root'

  if not path.exists(outdir):
    os.system('mkdir -p {}'.format(outdir))

  #treeProducer(infiles=infiles, outdir=outdir, outfilename=outfilename, lepton=lepton)
  treeProducer(infiles=infiles, outdir=outdir, outfilename=outfilename)
