import ROOT
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os
from os import path
from matplotlib.lines import Line2D


class Quantity(object):
  def __init__(self, name='', title='', label='', nbins=0, bin_min=0., bin_max=0.):
    self.name = name
    self.title = title
    self.label = label
    self.nbins = nbins
    self.bin_min = bin_min
    self.bin_max = bin_max


class DecayChain(object):
  def __init__(self, name, title):
    self.name = name
    self.title = title


def get_legend_name(branch_name):
  dictionary = {}
  dictionary['to'] = r'$\to$'
  dictionary['Bs'] = r'$B_s$'
  dictionary['Bsstar'] = r'$B_s^*$'
  dictionary['Bc'] = r'B$_c$'
  dictionary['Bcstar'] = r'$B_c^*$'
  dictionary['Bsbar'] = r'$\overline{B}_s$'
  dictionary['Bsstarbar'] = r'$\overline{B}_s^*$'
  dictionary['Bcbar'] = r'$\overline{B}_c$'
  dictionary['Bcstarbar'] = r'$\overline{B}_c^*$'
  dictionary['phiphi'] = r'$\phi\phi$'
  parts = branch_name.split('_')
  leg_name = ''
  for part in parts:
    if 'chain' in part: continue
    leg_name += dictionary[part]

  return leg_name


def plot_decay_channels(version_label, isBs):
  print '\n'
  print '----------------------'
  print 'Will plot {} decay channels for ntuples {}'.format('Bs' if isBs else 'Bsbar', version_label)
  print '----------------------'

  filename = './outputfiles/{}/genTree.root'.format(version_label)

  f = ROOT.TFile.Open(filename, 'READ')
  tree = f.Get('tree')

  branches = tree.GetListOfBranches()
  decay_chains = []
  for branch in branches:
    name = branch.GetName()
    if 'chain_' not in name: continue
    if isBs and 'Bs_to_phiphi' not in name: continue
    if not isBs and 'Bsbar_to_phiphi' not in name: continue
    decay_chains.append(DecayChain(name, get_legend_name(name)))

  channels = {}
  channels_tmp = []

  mean_tot = 0
  for decay_chain in decay_chains:
    hist_name = 'hist_{}'.format(decay_chain.name)
    hist = ROOT.TH1D(hist_name, hist_name, 50, 0, 1.1)
    tree.Project(hist_name, decay_chain.name, '')
    hist.SetDirectory(0)
    mean_channel =  hist.GetMean()
    mean_tot += mean_channel
    channels_tmp.append(mean_channel)

  #if do_normalise:
  channels[''] = [channel/mean_tot for channel in channels_tmp] # trick of dict to remove label on y-axis
  #else:
  #  channels[''] = channel # trick of dict to remove label on y-axis

  print 'mean_tot = {}'.format(mean_tot)
  
  # plot the rates
  the_df = pd.DataFrame(channels, index=[f.title for f in decay_chains])
  the_df = the_df.transpose()

  ax = the_df.plot.barh(stacked=True, colormap='PiYG' if isBs else 'PRGn', figsize=(10, 7))#'Set3')#'Pastel1')#'OrRd')#,colormap='PiYG')
  if isBs:
    ax.set_title(r'$B_s \to \phi\phi$ decay channels')
  else:
    ax.set_title(r'$\overline{B}_s \to \phi\phi$ decay channels')
  ax.set_xlim(0, 1)
  ax.margins(y=0)

  ax.set_ylim(-0.3, 0.3)

  handles, labels = ax.get_legend_handles_labels()

  group1_header = Line2D([0], [0], color='none', label='Unmixed decays')
  group2_header = Line2D([0], [0], color='none', label='\nMixed decays')

  grouped_handles = [group1_header] + handles[:6] + [group2_header] + handles[6:]
  grouped_labels = ['Unmixed decays'] + labels[:6] + ['\nMixed decays'] + labels[6:]

  lgd = ax.legend(grouped_handles, grouped_labels, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)

  if isBs:
    figname = 'channels_Bs'
  else:
    figname = 'channels_Bsbar'

  ax.figure.savefig('myPlots/mixing/{}/{}.png'.format(version_label, figname),bbox_extra_artists=(lgd,), bbox_inches='tight')
  ax.figure.savefig('myPlots/mixing/{}/{}.pdf'.format(version_label, figname),bbox_extra_artists=(lgd,), bbox_inches='tight')

  print '--> myPlots/mixing/{}/{}.png created'.format(version_label, figname)


def plot_mixing_channels(version_label):
  print '\n'
  print '----------------------'
  print 'Will plot mixing channels for ntuples {}'.format(version_label)
  print '----------------------'

  filename = './outputfiles/{}/genTree.root'.format(version_label)

  f = ROOT.TFile.Open(filename, 'READ')
  tree = f.Get('tree')

  decay_chains_Bs_unmixed = []
  decay_chains_Bs_mixed = []
  decay_chains_Bsbar_unmixed = []
  decay_chains_Bsbar_mixed = []

  branches = tree.GetListOfBranches()
  for branch in branches:
    name = branch.GetName()
    if 'chain_' not in name: continue

    if 'Bs_to_phiphi' in name and 'Bsbar_to_Bs_to_phiphi' not in name:
      decay_chains_Bs_unmixed.append(DecayChain(name, get_legend_name(name)))
    elif 'Bs_to_phiphi' in name and 'Bsbar_to_Bs_to_phiphi' in name:
      decay_chains_Bs_mixed.append(DecayChain(name, get_legend_name(name)))
    elif 'Bsbar_to_phiphi' in name and 'Bs_to_Bsbar_to_phiphi' not in name:
      decay_chains_Bsbar_unmixed.append(DecayChain(name, get_legend_name(name)))
    elif 'Bsbar_to_phiphi' in name and 'Bs_to_Bsbar_to_phiphi' in name:
      decay_chains_Bsbar_mixed.append(DecayChain(name, get_legend_name(name)))
    else:
      print 'WARNING - unknown decay chain: {}'.format(name)

  channels = {}
  #channels_tmp = []

  mean_Bs_unmixed_tot = 0
  for decay_chain in decay_chains_Bs_unmixed:
    hist_name = 'hist_{}'.format(decay_chain.name)
    hist = ROOT.TH1D(hist_name, hist_name, 50, 0, 1.1)
    tree.Project(hist_name, decay_chain.name, '')
    hist.SetDirectory(0)
    mean_channel =  hist.GetMean()
    mean_Bs_unmixed_tot += mean_channel

  mean_Bs_mixed_tot = 0
  for decay_chain in decay_chains_Bs_mixed:
    hist_name = 'hist_{}'.format(decay_chain.name)
    hist = ROOT.TH1D(hist_name, hist_name, 50, 0, 1.1)
    tree.Project(hist_name, decay_chain.name, '')
    hist.SetDirectory(0)
    mean_channel =  hist.GetMean()
    mean_Bs_mixed_tot += mean_channel

  mean_Bsbar_unmixed_tot = 0
  for decay_chain in decay_chains_Bsbar_unmixed:
    hist_name = 'hist_{}'.format(decay_chain.name)
    hist = ROOT.TH1D(hist_name, hist_name, 50, 0, 1.1)
    tree.Project(hist_name, decay_chain.name, '')
    hist.SetDirectory(0)
    mean_channel =  hist.GetMean()
    mean_Bsbar_unmixed_tot += mean_channel

  mean_Bsbar_mixed_tot = 0
  for decay_chain in decay_chains_Bsbar_mixed:
    hist_name = 'hist_{}'.format(decay_chain.name)
    hist = ROOT.TH1D(hist_name, hist_name, 50, 0, 1.1)
    tree.Project(hist_name, decay_chain.name, '')
    hist.SetDirectory(0)
    mean_channel =  hist.GetMean()
    mean_Bsbar_mixed_tot += mean_channel

  channels[''] = [mean_Bs_unmixed_tot, mean_Bs_mixed_tot, mean_Bsbar_unmixed_tot, mean_Bsbar_mixed_tot]

  print 'mean_tot = {} + {} + {} + {} = {}'.format(mean_Bs_unmixed_tot, mean_Bs_mixed_tot, mean_Bsbar_unmixed_tot, mean_Bsbar_mixed_tot, mean_Bs_unmixed_tot + mean_Bs_mixed_tot + mean_Bsbar_unmixed_tot + mean_Bsbar_mixed_tot)
  
  # plot the rates
  the_df = pd.DataFrame(channels, index=[r'unmixed $B_s \to \phi\phi$ decays', r'mixed $B_s \to \phi\phi$ decays', r'unmixed $\overline{B}_s \to \phi\phi$ decays', r'mixed $\overline{B}_s \to \phi\phi$ decays'])
  the_df = the_df.transpose()

  ax = the_df.plot.barh(stacked=True, colormap='Pastel1', figsize=(10, 7))#'Set3')#'Pastel1')#'OrRd')#,colormap='PiYG')
  ax.set_title(r'')
  ax.set_xlim(0, 1)
  ax.margins(y=0)
  ax.set_ylim(-0.3, 0.3)

  ax.text(0.125, 0, round(the_df[r'unmixed $B_s \to \phi\phi$ decays'], 2), ha='center', va='center',color='black')
  ax.text(0.375, 0, round(the_df[r'mixed $B_s \to \phi\phi$ decays'], 2), ha='center', va='center',color='black')
  ax.text(0.625, 0, round(the_df[r'unmixed $\overline{B}_s \to \phi\phi$ decays'], 2), ha='center', va='center',color='black')
  ax.text(0.875, 0, round(the_df[r'mixed $\overline{B}_s \to \phi\phi$ decays'], 2), ha='center', va='center',color='black')

  lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)

  ax.figure.savefig('myPlots/mixing/{}/channels_mixing.png'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')
  ax.figure.savefig('myPlots/mixing/{}/channels_mixing.pdf'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')

  print '--> myPlots/mixing/{}/channels_mixing.png created'.format(version_label)


def plot_trgmu_mother(version_label):
  print '----------------------'
  print 'Will plot trigger muon mother for ntuples {}'.format(version_label)
  print '----------------------'

  trgmu_mother_pdgid = Quantity('fabs(trgmu_mother_pdgid)', 'trgmu_mother_pdgid', 'trgmu_mother_pdgid', 1000000, 0, 1000000)
  #[511, 521, 443, 411, 431, 15, 531, 13, 421, 5122, 5332, 553, 5132, 5232, 223, 4122, 541, 100443, 4132, 22, 4232, 221]
  # B mesons
  # 511: B0
  # 521: B+
  # 531: Bs
  # 541: Bc

  # B baryons
  # 5122: Lambda_b
  # 5332: Omega_b
  # 5132: Xi_b-
  # 5232: Xi_b0

  # bbbar meson
  # 553: Upsilon(1S)

  # D mesons
  # 411: D+
  # 431: Ds+
  # 421: D0

  # D baryons
  # 4122: Lambda_c
  # 4132: Lambda_c0
  # 4232: Lambda_c+

  # ccbar mesons charmonium
  # 443: J/psi
  # 100443: psi(2S)

  # light flavoured mesons
  # 223: omega(782)
  # 221: eta

  # other
  # 22: gamma
  # 13: mu
  
  filename = './outputfiles/{}/genTree.root'.format(version_label)

  f = ROOT.TFile.Open(filename, 'READ')
  tree = f.Get('tree')

  trgmu_mother_species = {}

  hist_name = 'hist_trgmu_mother_pdgid'
  hist = ROOT.TH1D(hist_name, hist_name, trgmu_mother_pdgid.nbins, trgmu_mother_pdgid.bin_min, trgmu_mother_pdgid.bin_max)
  tree.Project(hist_name, 'fabs(trgmu_mother_pdgid)', '')
  hist.SetDirectory(0)

  is_B_meson = (hist.GetBinContent(hist.GetXaxis().FindBin(511)) + hist.GetBinContent(hist.GetXaxis().FindBin(521)) + hist.GetBinContent(hist.GetXaxis().FindBin(531)) + hist.GetBinContent(hist.GetXaxis().FindBin(541))) / hist.GetEntries()
  is_B_baryon = (hist.GetBinContent(hist.GetXaxis().FindBin(5122)) + hist.GetBinContent(hist.GetXaxis().FindBin(5332)) + hist.GetBinContent(hist.GetXaxis().FindBin(5132)) + hist.GetBinContent(hist.GetXaxis().FindBin(5232))) / hist.GetEntries()
  is_bbbar_meson = (hist.GetBinContent(hist.GetXaxis().FindBin(552))) / hist.GetEntries()
  is_D_meson = (hist.GetBinContent(hist.GetXaxis().FindBin(411)) + hist.GetBinContent(hist.GetXaxis().FindBin(431)) + hist.GetBinContent(hist.GetXaxis().FindBin(421))) / hist.GetEntries()
  is_D_baryon = (hist.GetBinContent(hist.GetXaxis().FindBin(4122)) + hist.GetBinContent(hist.GetXaxis().FindBin(4132)) + hist.GetBinContent(hist.GetXaxis().FindBin(4232))) / hist.GetEntries()
  is_ccbar_meson = (hist.GetBinContent(hist.GetXaxis().FindBin(443)) + hist.GetBinContent(hist.GetXaxis().FindBin(100443))) / hist.GetEntries()
  is_light_meson = (hist.GetBinContent(hist.GetXaxis().FindBin(223)) + hist.GetBinContent(hist.GetXaxis().FindBin(221))) / hist.GetEntries()
  is_other = 1 - (is_B_meson + is_B_baryon + is_bbbar_meson + is_D_meson + is_D_baryon + is_ccbar_meson + is_light_meson)

  trgmu_mother_species[''] = [is_B_meson, is_B_baryon, is_bbbar_meson, is_D_meson, is_D_baryon, is_ccbar_meson, is_light_meson, is_other]

  # plot the rates
  the_df = pd.DataFrame(trgmu_mother_species, index=['B mesons', 'B baryons', r'$b\overline{b}$ mesons', 'D mesons', 'D baryons', r'$c\overline{c}$ mesons', 'Light flavoured mesons', 'Other'])
  the_df = the_df.transpose()

  ax = the_df.plot.barh(stacked=True, colormap='Paired', figsize=(10, 7))#'Set3')#'Pastel1')#'OrRd')#,colormap='PiYG')
  lgd = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  ax.set_title(r'Trigger muon mother')
  ax.set_xlim(0, 1)
  ax.margins(y=0)
  ax.set_ylim(-0.3, 0.3)

  ax.figure.savefig('myPlots/mixing/{}/trgmu_mother_species.png'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')
  ax.figure.savefig('myPlots/mixing/{}/trgmu_mother_species.pdf'.format(version_label),bbox_extra_artists=(lgd,), bbox_inches='tight')

  print '--> myPlots/mixing/{}/trgmu_mother_species.png created'.format(version_label)



if __name__ == "__main__":

  #version_label = '102X_crab_trgmu_filter'
  #version_label = 'test_filter_phi_v2' 
  version_label = 'test_fragment_v2' 

  outdirname = './myPlots/mixing/{}'.format(version_label)
  if not path.exists(outdirname):
    os.system('mkdir -p {}'.format(outdirname))

  #plot_decay_channels(version_label=version_label, isBs=True)
  #plot_decay_channels(version_label=version_label, isBs=False)
  #plot_mixing_channels(version_label=version_label)
  plot_trgmu_mother(version_label=version_label)
