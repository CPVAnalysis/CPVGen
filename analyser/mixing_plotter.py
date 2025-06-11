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



if __name__ == "__main__":

  #version_label = '102X_crab_trgmu_filter'
  version_label = 'test_filter_phi_v2' 

  outdirname = './myPlots/mixing/{}'.format(version_label)
  if not path.exists(outdirname):
    os.system('mkdir -p {}'.format(outdirname))

  plot_decay_channels(version_label=version_label, isBs=True)
  plot_decay_channels(version_label=version_label, isBs=False)
  plot_mixing_channels(version_label=version_label)
