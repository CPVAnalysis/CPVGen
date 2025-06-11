import os
import sys
import glob
import ROOT
from DataFormats.FWLite import Events, Handle
from collections import OrderedDict

def get_list_daughters(mother, pdgid):
    list_daughters = [mother.daughter(jj) for jj in range(mother.numberOfDaughters()) if abs(mother.daughter(jj).pdgId()) == pdgid]
    sorted_list_daughters = sorted([ii for ii in list_daughters], key = lambda x : x.pt(), reverse=True)

    return sorted_list_daughters


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


def analyse(infiles):
  # get the files
  print 'loading the file ...'
  files = glob.glob(infiles)
  events = Events(files)
  print 'nevts: ',events.size()
  print '... done!'

  # get the handles
  handles = OrderedDict()
  handles['genP'] = ('genParticles' , Handle('std::vector<reco::GenParticle>'))

  decay_chains_Bs_mixed = []
  decay_chains_Bs_unmixed = []
  decay_chains_Bsbar_mixed = []
  decay_chains_Bsbar_unmixed = []

  for i, event in enumerate(events):
    #if float(i)<6000: continue
    #print '\n\n Event {}'.format(i)

    # access the handles
    for k, v in handles.iteritems():
        event.getByLabel(v[0], v[1])
        setattr(event, k, v[1].product())

    if i%1000==0:
        percentage = float(i)/events.size()*100.
        print '\t===> processing %d / %d event \t completed %.1f%s' %(i, events.size(), percentage, '%')

    # get the Bs mesons
    the_Bs_mesons = [ip for ip in event.genP if abs(ip.pdgId())==531 and ip.isLastCopy()]

    # find Bs meson that has two phi daughters
    idx_Bs = -99
    count_Bs_has_two_phi_daughters = 0
    for iBs, the_Bs_meson in enumerate(the_Bs_mesons):
      number_phi_daughters = len(get_list_daughters(mother=the_Bs_meson, pdgid=333))
      if number_phi_daughters == 2:
        idx_Bs = iBs
        count_Bs_has_two_phi_daughters += 1

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

    if event.the_Bs.pdgId() == 531 and len(decay_chain) == 0:
      if decay_chain not in decay_chains_Bs_unmixed: decay_chains_Bs_unmixed.append(decay_chain)

    elif event.the_Bs.pdgId() == 531 and decay_chain[0] != -531:
      if decay_chain not in decay_chains_Bs_unmixed: decay_chains_Bs_unmixed.append(decay_chain)

    elif event.the_Bs.pdgId() == 531 and decay_chain[0] == -531:
      if decay_chain not in decay_chains_Bs_mixed: decay_chains_Bs_mixed.append(decay_chain)

    elif event.the_Bs.pdgId() == -531 and len(decay_chain) == 0:
      if decay_chain not in decay_chains_Bsbar_unmixed: decay_chains_Bsbar_unmixed.append(decay_chain)

    elif event.the_Bs.pdgId() == -531 and decay_chain[0] != 531:
      if decay_chain not in decay_chains_Bsbar_unmixed: decay_chains_Bsbar_unmixed.append(decay_chain)

    elif event.the_Bs.pdgId() == -531 and decay_chain[0] == 531:
      if decay_chain not in decay_chains_Bsbar_mixed: decay_chains_Bsbar_mixed.append(decay_chain)

  decay_chains_Bs_unmixed_sorted = sorted(decay_chains_Bs_unmixed, key=len)
  decay_chains_Bs_mixed_sorted = sorted(decay_chains_Bs_mixed, key=len)
  decay_chains_Bsbar_unmixed_sorted = sorted(decay_chains_Bsbar_unmixed, key=len)
  decay_chains_Bsbar_mixed_sorted = sorted(decay_chains_Bsbar_mixed, key=len)

  return decay_chains_Bs_unmixed_sorted, decay_chains_Bs_mixed_sorted, decay_chains_Bsbar_unmixed_sorted, decay_chains_Bsbar_mixed_sorted


def get_name(decay_chain, isBs):
  dictionary = {}
  dictionary[531] = 'Bs_to_'
  dictionary[-531] = 'Bsbar_to_'
  dictionary[533] = 'Bsstar_to_'
  dictionary[-533] = 'Bsstarbar_to_'
  dictionary[541] = 'Bc_to_'
  dictionary[-541] = 'Bcbar_to_'
  dictionary[543] = 'Bcstar_to_'
  dictionary[-543] = 'Bcstarbar_to_'

  decay_chain.reverse()
  decay_chain_name = 'chain_'
  for pdgid in decay_chain:
    decay_chain_name += dictionary[pdgid]

  decay_chain_name += '{bs}_to_phiphi'.format(bs='Bs' if isBs else 'Bsbar')
  decay_chain.reverse()

  return decay_chain_name


def write_branches(outfile, decay_chains, isBs, isMixed):
  #text = '\n\nbranches = ['
  text = '\n\n    # {mx} {bs} decays'.format(mx='mixed' if isMixed else 'unmixed', bs='Bs' if isBs else 'Bsbar')
  for decay_chain in decay_chains:
    text += "\n    '{}',".format(get_name(decay_chain, isBs))
    chain_Bsbar_to_phiphi = 0

  #text += '\n]'
  outfile.write(text)
  print text

  
def write_initialisation(outfile, decay_chains, isBs, isMixed):
  text = '\n\n    # {mx} {bs} decays'.format(mx='mixed' if isMixed else 'unmixed', bs='Bs' if isBs else 'Bsbar')
  for decay_chain in decay_chains:
    text += "\n    {} = 0".format(get_name(decay_chain, isBs))
  outfile.write(text)
  print text


def write_cond(outfile, decay_chains, isBs):
  text = '\n\n    # {bs} decays'.format(bs='Bs' if isBs else 'Bsbar')
  text += '\n    if event.the_Bs.pdgId() == {}:'.format('531' if isBs else '-531')
  for idx, decay_chain in enumerate(decay_chains):
    name = get_name(decay_chain, isBs)
    if idx == 0:
      text += '\n'.join([
          '\n      if decay_chain == {}:'.format(decay_chain),
          '        {} = 1'.format(name),
          ])
    else:
      text += '\n'.join([
          '\n      elif decay_chain == {}:'.format(decay_chain),
          '        {} = 1'.format(name),
          ])

  text += '\n      else:'
  if isBs:
    text += '\n        print "WARNING: unknown decay chain for Bs: {}".format(decay_chain)'
  else:
    text += '\n        print "WARNING: unknown decay chain for Bsbar: {}".format(decay_chain)'
  outfile.write(text)
  print text


def write_fill(outfile, decay_chains, isBs, isMixed):
  text = '\n\n    # {mx} {bs} decays'.format(mx='mixed' if isMixed else 'unmixed', bs='Bs' if isBs else 'Bsbar')
  for decay_chain in decay_chains:
    name = get_name(decay_chain, isBs)
    text += "\n    tofill['{n}'] = {n}".format(n=name)

  outfile.write(text)
  print text


def process(infiles):
  decay_chains_Bs_unmixed, decay_chains_Bs_mixed, decay_chains_Bsbar_unmixed, decay_chains_Bsbar_mixed = analyse(infiles=infiles)
  #decay_chains_Bs_unmixed = [[], [533], [-531, -533], [-531], [-531, -533, -541, -543], [541, 543], [-531, -533, -541]]

  outfile = open('decay_chains.txt', 'w+')
  outfile.write('\n\n ########## Branches ##########')
  write_branches(outfile, decay_chains_Bs_unmixed, isBs=True, isMixed=False)
  write_branches(outfile, decay_chains_Bs_mixed, isBs=True, isMixed=True)
  write_branches(outfile, decay_chains_Bsbar_unmixed, isBs=False, isMixed=False)
  write_branches(outfile, decay_chains_Bsbar_mixed, isBs=False, isMixed=True)

  outfile.write('\n\n ########## Initialisation ##########')
  write_initialisation(outfile, decay_chains_Bs_unmixed, isBs=True, isMixed=False)
  write_initialisation(outfile, decay_chains_Bs_mixed, isBs=True, isMixed=True)
  write_initialisation(outfile, decay_chains_Bsbar_unmixed, isBs=False, isMixed=False)
  write_initialisation(outfile, decay_chains_Bsbar_mixed, isBs=False, isMixed=True)

  outfile.write('\n\n ########## If condition ##########')
  write_cond(outfile, decay_chains_Bs_unmixed+decay_chains_Bs_mixed, isBs=True)
  write_cond(outfile, decay_chains_Bsbar_unmixed+decay_chains_Bsbar_mixed, isBs=False)

  outfile.write('\n\n ########## Fill ##########')
  write_fill(outfile, decay_chains_Bs_unmixed, isBs=True, isMixed=False)
  write_fill(outfile, decay_chains_Bs_mixed, isBs=True, isMixed=True)
  write_fill(outfile, decay_chains_Bsbar_unmixed, isBs=False, isMixed=False)
  write_fill(outfile, decay_chains_Bsbar_mixed, isBs=False, isMixed=True)

  outfile.close()
  print '-> decay_chains.txt created'


if __name__ == "__main__":
  version_label = '102X_crab_trgmu_filter'
  indirectory = '/eos/user/a/anlyon/CPVGen/102X_crab_trgmu_filter/BsToPhiPhiTo4K/crab_102X_crab_trgmu_filter_BsToPhiPhiTo4K_20250212_225911/250212_215924/000*/'

  print '----------------------'
  print ' Analysing {}'.format(indirectory)
  print '----------------------'

  infiles = '{}/step1*root'.format(indirectory)

  process(infiles=infiles)
