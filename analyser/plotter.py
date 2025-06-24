import os
from os import path
import ROOT
import numpy as np
import math


class Quantity(object):
  def __init__(self, name='', title='', label='', nbins=0, bin_min=0., bin_max=0., do_log=False):
    self.name = name
    self.title = title
    self.label = label
    self.nbins = nbins
    self.bin_min = bin_min
    self.bin_max = bin_max
    self.do_log = do_log


def plot_kaons(inputfilename, outdir, quantity):

  f = ROOT.TFile.Open(inputfilename, 'READ') 
  tree = f.Get('tree')

  canv_name = 'canv'
  canv = ROOT.TCanvas(canv_name, canv_name, 800, 700)
  if quantity.do_log:
    canv.SetLogy()
  canv.SetLeftMargin(0.15)
  canv.SetBottomMargin(0.15)

  leg = ROOT.TLegend(0.65, 0.6, 0.85, 0.8)
  leg.SetTextSize(0.05)
  leg.SetLineColor(0)
  leg.SetFillColorAlpha(0, 0)
  leg.SetBorderSize(0)

  nbins = quantity.nbins
  bin_min = quantity.bin_min
  bin_max = quantity.bin_max
      
  hist1_name = 'hist1'
  hist1 = ROOT.TH1D(hist1_name, hist1_name, nbins, bin_min, bin_max)
  hist1.SetLineWidth(3)
  hist1.SetLineColor(6)

  hist2_name = 'hist2'
  hist2 = ROOT.TH1D(hist2_name, hist2_name, nbins, bin_min, bin_max)
  hist2.SetLineWidth(3)
  hist2.SetLineColor(4)

  hist3_name = 'hist3'
  hist3 = ROOT.TH1D(hist3_name, hist3_name, nbins, bin_min, bin_max)
  hist3.SetLineWidth(3)
  hist3.SetLineColor(2)

  hist4_name = 'hist4'
  hist4 = ROOT.TH1D(hist4_name, hist4_name, nbins, bin_min, bin_max)
  hist4.SetLineWidth(3)
  hist4.SetLineColor(8)

  tree.Project(hist1_name, 'k1_{}'.format(quantity.name))
  tree.Project(hist2_name, 'k2_{}'.format(quantity.name))
  tree.Project(hist3_name, 'k3_{}'.format(quantity.name))
  tree.Project(hist4_name, 'k4_{}'.format(quantity.name))

  leg.AddEntry(hist1, 'k_{1}')
  leg.AddEntry(hist2, 'k_{2}')
  leg.AddEntry(hist3, 'k_{3}')
  leg.AddEntry(hist4, 'k_{4}')

  hist1.SetTitle('')
  hist1.GetXaxis().SetTitle(quantity.title)
  hist1.GetXaxis().SetLabelSize(0.04)
  hist1.GetXaxis().SetTitleSize(0.047)
  hist1.GetXaxis().SetTitleOffset(1.)
  hist1.GetYaxis().SetTitle('Events / Bin')
  hist1.GetYaxis().SetLabelSize(0.04)
  hist1.GetYaxis().SetTitleSize(0.047)
  hist1.GetYaxis().SetTitleOffset(1.5)
  #hist1.GetYaxis().SetRangeUser(0.5, 1e5)

  hist1.Draw('histo')
  hist2.Draw('histo same')
  hist3.Draw('histo same')
  hist4.Draw('histo same')

  leg.Draw()
      
  ROOT.gStyle.SetOptStat(0)

  fig_name = 'comparison_kaons_{}'.format(quantity.label)
  if quantity.do_log: fig_name += '_log'

  canv.SaveAs(outdir + '/' + fig_name + '.png')
  canv.SaveAs(outdir + '/' + fig_name + '.pdf')


def plot_phis(inputfilename, outdir, quantity):
  f = ROOT.TFile.Open(inputfilename, 'READ') 
  tree = f.Get('tree')

  canv_name = 'canv'
  canv = ROOT.TCanvas(canv_name, canv_name, 800, 700)
  if quantity.do_log:
    canv.SetLogy()
  canv.SetLeftMargin(0.15)
  canv.SetBottomMargin(0.15)

  leg = ROOT.TLegend(0.65, 0.6, 0.85, 0.8)
  leg.SetTextSize(0.05)
  leg.SetLineColor(0)
  leg.SetFillColorAlpha(0, 0)
  leg.SetBorderSize(0)

  nbins = quantity.nbins
  bin_min = quantity.bin_min
  bin_max = quantity.bin_max
      
  hist1_name = 'hist1'
  hist1 = ROOT.TH1D(hist1_name, hist1_name, nbins, bin_min, bin_max)
  hist1.SetLineWidth(3)
  hist1.SetLineColor(6)

  hist2_name = 'hist2'
  hist2 = ROOT.TH1D(hist2_name, hist2_name, nbins, bin_min, bin_max)
  hist2.SetLineWidth(3)
  hist2.SetLineColor(4)

  hist3_name = 'hist3'
  hist3 = ROOT.TH1D(hist3_name, hist3_name, nbins, bin_min, bin_max)
  hist3.SetLineWidth(3)
  hist3.SetLineColor(2)

  tree.Project(hist1_name, 'phi1_{}'.format(quantity.name))
  tree.Project(hist2_name, 'phi2_{}'.format(quantity.name))
  tree.Project(hist3_name, 'Bs_{}'.format(quantity.name))

  leg.AddEntry(hist1, '#phi_{1}')
  leg.AddEntry(hist2, '#phi_{2}')
  leg.AddEntry(hist3, 'B_{s}')

  hist1.SetTitle('')
  hist1.GetXaxis().SetTitle(quantity.title)
  hist1.GetXaxis().SetLabelSize(0.04)
  hist1.GetXaxis().SetTitleSize(0.047)
  hist1.GetXaxis().SetTitleOffset(1.)
  hist1.GetYaxis().SetTitle('Events / Bin')
  hist1.GetYaxis().SetLabelSize(0.04)
  hist1.GetYaxis().SetTitleSize(0.047)
  hist1.GetYaxis().SetTitleOffset(1.5)
  #hist1.GetYaxis().SetRangeUser(0.5, 1e5)

  hist1.Draw('histo')
  hist2.Draw('histo same')
  hist3.Draw('histo same')

  leg.Draw()
      
  ROOT.gStyle.SetOptStat(0)

  fig_name = 'comparison_phis_{}'.format(quantity.label)
  if quantity.do_log: fig_name += '_log'

  canv.SaveAs(outdir + '/' + fig_name + '.png')
  canv.SaveAs(outdir + '/' + fig_name + '.pdf')


def plot_mixing(inputfilename, outdir, quantity):
  f = ROOT.TFile.Open(inputfilename, 'READ') 
  tree = f.Get('tree')

  canv_name = 'canv_mixing_' + quantity.label
  canv = ROOT.TCanvas(canv_name, canv_name, 800, 700)
  if quantity.do_log:
    canv.SetLogy()
  canv.SetLeftMargin(0.15)
  canv.SetBottomMargin(0.15)

  leg = ROOT.TLegend(0.65, 0.6, 0.75, 0.8)
  leg.SetTextSize(0.05)
  leg.SetLineColor(0)
  leg.SetFillColorAlpha(0, 0)
  leg.SetBorderSize(0)

  nbins = quantity.nbins
  bin_min = quantity.bin_min
  bin_max = quantity.bin_max
      
  hist1_name = 'hist1_mixing_' + quantity.label
  hist1 = ROOT.TH1D(hist1_name, hist1_name, nbins, bin_min, bin_max)
  hist1.SetLineWidth(3)
  hist1.SetLineColor(6)

  hist2_name = 'hist2_mixing_' + quantity.label
  hist2 = ROOT.TH1D(hist2_name, hist2_name, nbins, bin_min, bin_max)
  hist2.SetLineWidth(3)
  hist2.SetLineColor(4)

  tree.Project(hist1_name, quantity.name, 'is_unmixed==1')
  tree.Project(hist2_name, quantity.name, 'is_mixed==1')

  leg.AddEntry(hist1, 'unmixed')
  leg.AddEntry(hist2, 'mixed')

  frame = hist1.Clone('frame')
  frame.SetTitle('')
  frame.GetXaxis().SetTitle(quantity.title)
  frame.GetXaxis().SetLabelSize(0.04)
  frame.GetXaxis().SetTitleSize(0.047)
  frame.GetXaxis().SetTitleOffset(1.)
  frame.GetYaxis().SetTitle('Events / Bin')
  frame.GetYaxis().SetLabelSize(0.04)
  frame.GetYaxis().SetTitleSize(0.047)
  frame.GetYaxis().SetTitleOffset(1.5)

  frame.Draw()
  hist1.Draw('histo same')
  hist2.Draw('histo same')

  leg.Draw()
      
  ROOT.gStyle.SetOptStat(0)

  fig_name = 'comparison_mixing_{}'.format(quantity.label)
  if quantity.do_log:
    fig_name += '_log'

  canv.SaveAs(outdir + '/' + fig_name + '.png')
  canv.SaveAs(outdir + '/' + fig_name + '.pdf')


def plot_pdgid(inputfilename, outdir, quantity):
  f = ROOT.TFile.Open(inputfilename, 'READ') 
  tree = f.Get('tree')

  canv_name = 'canv_pdgid_' + quantity.label
  canv = ROOT.TCanvas(canv_name, canv_name, 800, 700)
  if quantity.do_log:
    canv.SetLogy()
  canv.SetLeftMargin(0.15)
  canv.SetBottomMargin(0.15)

  leg = ROOT.TLegend(0.65, 0.6, 0.75, 0.8)
  leg.SetTextSize(0.05)
  leg.SetLineColor(0)
  leg.SetFillColorAlpha(0, 0)
  leg.SetBorderSize(0)

  nbins = quantity.nbins
  bin_min = quantity.bin_min
  bin_max = quantity.bin_max
      
  hist1_name = 'hist1_pdgid_' + quantity.label
  hist1 = ROOT.TH1D(hist1_name, hist1_name, nbins, bin_min, bin_max)
  hist1.SetLineWidth(3)
  hist1.SetLineColor(6)

  hist2_name = 'hist2_pdgid_' + quantity.label
  hist2 = ROOT.TH1D(hist2_name, hist2_name, nbins, bin_min, bin_max)
  hist2.SetLineWidth(3)
  hist2.SetLineColor(4)

  tree.Project(hist1_name, quantity.name, 'Bs_pdgid>0')
  tree.Project(hist2_name, quantity.name, 'Bs_pdgid<0')

  leg.AddEntry(hist1, 'B_{s}')
  leg.AddEntry(hist2, '#bar{B}_{s}')

  frame = hist1.Clone('frame')
  frame.SetTitle('')
  frame.GetXaxis().SetTitle(quantity.title)
  frame.GetXaxis().SetLabelSize(0.04)
  frame.GetXaxis().SetTitleSize(0.047)
  frame.GetXaxis().SetTitleOffset(1.)
  frame.GetYaxis().SetTitle('Events / Bin')
  frame.GetYaxis().SetLabelSize(0.04)
  frame.GetYaxis().SetTitleSize(0.047)
  frame.GetYaxis().SetTitleOffset(1.5)

  frame.Draw()
  hist1.Draw('histo same')
  hist2.Draw('histo same')

  leg.Draw()
      
  ROOT.gStyle.SetOptStat(0)

  fig_name = 'comparison_pdgid_{}'.format(quantity.label)
  if quantity.do_log:
    fig_name += '_log'

  canv.SaveAs(outdir + '/' + fig_name + '.png')
  canv.SaveAs(outdir + '/' + fig_name + '.pdf')


def plot_asymmetry_graph_mixing(inputfilename, outdir):

  ct = Quantity('Bs_ct_3D_ps', 'ct 3D (ps)', 'Bs_ct_3D_ps', 50, 0, 7)
  #ct = Quantity('Bs_ct_ps_alasin2b', 'ct (ps)', 'Bs_ct_ps', 50, 0, 7)
  step_ct = 0.1
  ct_ranges = np.linspace(0, 7, 7 / step_ct + 1)

  graph = ROOT.TGraphAsymmErrors()

  for ict, ct_range in enumerate(ct_ranges):
    if ict == len(ct_ranges)-1: break
    ct_range_min = ct_ranges[ict]
    ct_range_max = ct_ranges[ict+1]
    #print '{} {}'.format(ct_range_min, ct_range_max)
    f = ROOT.TFile.Open(inputfilename, 'READ') 
    tree = f.Get('tree')

    nbins = ct.nbins
    bin_min = ct.bin_min
    bin_max = ct.bin_max
        
    hist1_name = 'hist1_mixing_{}'.format(ct_range)
    hist1 = ROOT.TH1D(hist1_name, hist1_name, nbins, bin_min, bin_max)
    hist1.SetLineWidth(3)
    hist1.SetLineColor(6)

    hist2_name = 'hist2_mixing_{}'.format(ct_range)
    hist2 = ROOT.TH1D(hist2_name, hist2_name, nbins, bin_min, bin_max)
    hist2.SetLineWidth(3)
    hist2.SetLineColor(4)

    range_selection = '{ct} > {mi} && {ct} < {ma}'.format(ct=ct.name, mi=ct_range_min, ma=ct_range_max)
    tree.Project(hist1_name, ct.name, 'is_unmixed==1 && {}'.format(range_selection))
    tree.Project(hist2_name, ct.name, 'is_mixed==1 && {}'.format(range_selection))

    #err_unmixed = ROOT.double(0.)
    #n_unmixed = hist1.IntegralAndError(0, 1300000, err_unmixed)
    n_unmixed = hist1.GetEntries()
    err_unmixed = math.sqrt(n_unmixed)

    #err_mixed = ROOT.double(0.)
    #n_mixed = hist2.IntegralAndError(0, 1300000, err_mixed)
    n_mixed = hist2.GetEntries()
    err_mixed = math.sqrt(n_unmixed)

    asymmetry = (n_unmixed - n_mixed) / (n_unmixed + n_mixed)
    err_asymmetry = abs(asymmetry) * (err_unmixed / n_unmixed + err_mixed / n_mixed)

    point = graph.GetN()
    graph.SetPoint(point, (ct_range_min + ct_range_max)/2., asymmetry)
    graph.SetPointError(point, 0, 0, err_asymmetry, err_asymmetry) 
    #graph.SetPointError(point, 0, 0, 0, 0) 


  canv_name = 'canv'
  canv = ROOT.TCanvas(canv_name, canv_name, 800, 700)
  #canv.SetLogy()
  canv.SetLeftMargin(0.15)
  canv.SetBottomMargin(0.15)

  graph.SetLineColor(1)
  graph.SetMarkerColor(2)
  graph.SetMarkerStyle(20)
  graph.GetXaxis().SetTitle('ct (ps)')
  graph.GetXaxis().SetLabelSize(0.037)
  graph.GetXaxis().SetTitleSize(0.042)
  graph.GetXaxis().SetTitleOffset(1.1)
  graph.GetYaxis().SetTitle('Asymmetry')
  graph.GetYaxis().SetLabelSize(0.037)
  graph.GetYaxis().SetTitleSize(0.042)
  graph.GetYaxis().SetTitleOffset(1.5)

  graph.Draw('AP')

  fig_name = 'asymmetry_mixing'
  canv.SaveAs(outdir + '/' + fig_name + '.png')
  canv.SaveAs(outdir + '/' + fig_name + '.pdf')


def plot_asymmetry_graph_pdgid(inputfilename, outdir):

  ct = Quantity('Bs_ct_3D_ps', 'ct 3D (ps)', 'Bs_ct_3D_ps', 10000, 0, 7)
  step_ct = 0.1
  ct_ranges = np.linspace(0, 7, 7 / step_ct + 1)
  #ct_ranges = np.linspace(0, 4, 4 / step_ct + 1)

  graph = ROOT.TGraphAsymmErrors()

  for ict, ct_range in enumerate(ct_ranges):
    if ict == len(ct_ranges)-1: break
    ct_range_min = ct_ranges[ict]
    ct_range_max = ct_ranges[ict+1]
    print '{} {}'.format(ct_range_min, ct_range_max)
    f = ROOT.TFile.Open(inputfilename, 'READ') 
    tree = f.Get('tree')

    nbins = ct.nbins
    bin_min = ct.bin_min
    bin_max = ct.bin_max
        
    hist1_name = 'hist1_mixing_{}'.format(ct_range)
    hist1 = ROOT.TH1D(hist1_name, hist1_name, nbins, bin_min, bin_max)
    hist1.SetLineWidth(3)
    hist1.SetLineColor(6)

    hist2_name = 'hist2_mixing_{}'.format(ct_range)
    hist2 = ROOT.TH1D(hist2_name, hist2_name, nbins, bin_min, bin_max)
    hist2.SetLineWidth(3)
    hist2.SetLineColor(4)

    range_selection = '{ct} > {mi} && {ct} < {ma}'.format(ct=ct.name, mi=ct_range_min, ma=ct_range_max)
    tree.Project(hist1_name, ct.name, 'Bs_pdgid>0 && {}'.format(range_selection))
    tree.Project(hist2_name, ct.name, 'Bs_pdgid<0 && {}'.format(range_selection))

    #err_unmixed = ROOT.double(0.)
    #n_unmixed = hist1.IntegralAndError(0, 1300000, err_unmixed)
    n_unmixed = hist1.GetEntries()
    err_unmixed = math.sqrt(n_unmixed)

    #err_mixed = ROOT.double(0.)
    #n_mixed = hist2.IntegralAndError(0, 1300000, err_mixed)
    n_mixed = hist2.GetEntries()
    err_mixed = math.sqrt(n_unmixed)

    asymmetry = (n_unmixed - n_mixed) / (n_unmixed + n_mixed)
    err_asymmetry = abs(asymmetry) * (err_unmixed / n_unmixed + err_mixed / n_mixed)

    point = graph.GetN()
    graph.SetPoint(point, (ct_range_min + ct_range_max)/2., asymmetry)
    graph.SetPointError(point, 0, 0, err_asymmetry, err_asymmetry) 
    #graph.SetPointError(point, 0, 0, 0, 0) 


  canv_name = 'canv'
  canv = ROOT.TCanvas(canv_name, canv_name, 800, 700)
  #canv.SetLogy()
  canv.SetLeftMargin(0.15)
  canv.SetBottomMargin(0.15)

  graph.SetLineColor(1)
  graph.SetMarkerColor(2)
  graph.SetMarkerStyle(20)
  graph.GetXaxis().SetTitle('ct (ps)')
  graph.GetXaxis().SetLabelSize(0.037)
  graph.GetXaxis().SetTitleSize(0.042)
  graph.GetXaxis().SetTitleOffset(1.1)
  graph.GetYaxis().SetTitle('Asymmetry')
  graph.GetYaxis().SetLabelSize(0.037)
  graph.GetYaxis().SetTitleSize(0.042)
  graph.GetYaxis().SetTitleOffset(1.5)

  graph.Draw('AP')

  fig_name = 'asymmetry_pdgid'
  canv.SaveAs(outdir + '/' + fig_name + '.png')
  canv.SaveAs(outdir + '/' + fig_name + '.pdf')


def study_scale(quantity):

  #inputfilename1 = './outputfiles/test_fragment_v2/genTree.root' # 102X, scale1p0, pt>0.5GeV, pdl2014
  inputfilename1 = './outputfiles/test_scale1p0_2018/genTree.root' # 106X, scale1p0, pt>0.3GeV, pdl2020
  inputfilename2 = './outputfiles/test_scale5p0/genTree.root' # 106X, scape5p0, pt>0.3GeV, pdl2020

  f1 = ROOT.TFile.Open(inputfilename1, 'READ') 
  f2 = ROOT.TFile.Open(inputfilename2, 'READ') 

  tree1 = f1.Get('tree')
  tree2 = f2.Get('tree')

  canv_name = 'canv'
  canv = ROOT.TCanvas(canv_name, canv_name, 800, 700)

  pad_up = ROOT.TPad("pad_up","pad_up",0,0.25,1,1)
  pad_up.Draw()
  pad_down = ROOT.TPad("pad_down","pad_down",0,0,1,0.25)
  pad_down.SetBottomMargin(0.25)
  pad_down.Draw()


  pad_up.cd()

  if quantity.do_log:
    pad_up.SetLogy()
  canv.SetLeftMargin(0.15)
  canv.SetBottomMargin(0.15)

  leg = ROOT.TLegend(0.65, 0.6, 0.85, 0.8)
  leg.SetTextSize(0.05)
  leg.SetLineColor(0)
  leg.SetFillColorAlpha(0, 0)
  leg.SetBorderSize(0)

  nbins = quantity.nbins
  bin_min = quantity.bin_min
  bin_max = quantity.bin_max
      
  hist1_name = 'hist1'
  hist1 = ROOT.TH1D(hist1_name, hist1_name, nbins, bin_min, bin_max)
  hist1.SetLineWidth(3)
  hist1.SetLineColor(6)

  hist2_name = 'hist2'
  hist2 = ROOT.TH1D(hist2_name, hist2_name, nbins, bin_min, bin_max)
  hist2.SetLineWidth(3)
  hist2.SetLineColor(4)

  tree1.Project(hist1_name, format(quantity.name))
  tree2.Project(hist2_name, format(quantity.name))

  hist1.Scale(1./hist1.Integral())
  hist2.Scale(1./hist2.Integral())

  leg.AddEntry(hist1, 'scale = 1.0')
  leg.AddEntry(hist2, 'scale = 5.0')

  hist1.SetTitle('')
  hist1.GetXaxis().SetTitle(quantity.title)
  hist1.GetXaxis().SetLabelSize(0.04)
  hist1.GetXaxis().SetTitleSize(0.047)
  hist1.GetXaxis().SetTitleOffset(1.)
  hist1.GetYaxis().SetTitle('Normalised to unity')
  hist1.GetYaxis().SetLabelSize(0.04)
  hist1.GetYaxis().SetTitleSize(0.047)
  hist1.GetYaxis().SetTitleOffset(1.1)
  #hist1.GetYaxis().SetRangeUser(0.5, 1e5)

  hist1.Draw('histo')
  hist2.Draw('histo same')

  leg.Draw()

  pad_down.cd()

  hist_ratio = hist1.Clone('hist_ratio')
  hist_ratio.Divide(hist2)

  hist_ratio.SetLineWidth(2)
  hist_ratio.SetMarkerStyle(20)
  hist_ratio.SetTitle('')

  hist_ratio.GetXaxis().SetLabelSize(0.11)
  hist_ratio.GetXaxis().SetTitleSize(0.13)
  hist_ratio.GetXaxis().SetTitleOffset(0.8)
  hist_ratio.GetYaxis().SetTitle('Ratio')
  hist_ratio.GetYaxis().SetLabelSize(0.11)
  hist_ratio.GetYaxis().SetTitleSize(0.13)
  hist_ratio.GetYaxis().SetTitleOffset(0.345)
  hist_ratio.GetYaxis().SetNdivisions(6)
  hist_ratio.GetYaxis().SetRangeUser(0, 2)


  hist_ratio.Draw('PE')

  line = ROOT.TLine(quantity.bin_min, 1, quantity.bin_max, 1)
  line.SetLineColor(4)
  line.SetLineWidth(2)
  line.Draw('same')
      
  canv.cd()
  ROOT.gStyle.SetOptStat(0)

  outdir = './myPlots/study_scale'
  if not path.exists(outdir):
    os.system('mkdir -p {}'.format(outdir))
  fig_name = quantity.label

  if quantity.do_log: fig_name += '_log'

  canv.SaveAs(outdir + '/' + fig_name + '.png')
  canv.SaveAs(outdir + '/' + fig_name + '.pdf')



if __name__ == "__main__":

  do_plot_kaons = False
  do_plot_phis = False
  do_plot_mixing = False
  do_plot_pdgid = False
  do_plot_asymmetry_mixing = False
  do_plot_asymmetry_pdgid = False
  do_study_scale = True

  ROOT.gROOT.SetBatch(True)
  ROOT.TH1.SetDefaultSumw2()

  #version_label = 'crab_102X_crab_v0'
  version_label = '102X_crab_trgmu_filter'

  #inputfilename = './outputfiles/test_trgmu_v0/genTree.root'
  #inputfilename = './outputfiles/crab_102X_crab_v0/genTree.root'
  #inputfilename = './outputfiles/102X_crab_trgmu_filter/genTree_tot.root'
  inputfilename = './outputfiles/{}/genTree_tot.root'.format(version_label)

  outdir = './myPlots/{}'.format(version_label)
  if not path.exists(outdir):
    os.system('mkdir -p {}'.format(outdir))

  if do_plot_kaons: 
    quantities = [
      Quantity('pt', 'p_{T} (GeV)', 'pt', 50, 0, 50, do_log=True),
      Quantity('eta', '#eta', 'eta', 100, -7, 7, do_log=True),
      Quantity('phi', '#phi', 'phi', 100, -4, 4, do_log=False),
      Quantity('q', 'charge', 'charge', 5, -2, 2, do_log=False),
      Quantity('mass', 'mass (GeV)', 'mass', 50, 0, 1, do_log=False),
    ]

    for quantity in quantities:
      plot_kaons(inputfilename=inputfilename, outdir=outdir, quantity=quantity)

  if do_plot_phis: 
    quantities = [
      Quantity('pt', 'p_{T} (GeV)', 'pt', 50, 0, 50, do_log=True),
      Quantity('eta', '#eta', 'eta', 100, -7, 7, do_log=True),
      Quantity('phi', '#phi', 'phi', 100, -4, 4, do_log=False),
      Quantity('q', 'charge', 'charge', 5, -2, 2, do_log=False),
      Quantity('mass', 'mass (GeV)', 'mass', 100, 0.6, 5.5, do_log=False),
    ]

    for quantity in quantities:
      plot_phis(inputfilename=inputfilename, outdir=outdir, quantity=quantity)

  quantities = [
    Quantity('Bs_ct_3D_ps', 'ct (ps)', 'Bs_ct_3D_ps', 100, 0, 7, do_log=False),
    Quantity('Bs_ct_3D_ps', 'ct (ps)', 'Bs_ct_3D_ps', 100, 0, 7, do_log=True),
    #Quantity('Bs_ct_2D_ps', 'ct 2D (ps)', 'Bs_ct_2D_ps', 50, 0, 7, do_log=False),
    #Quantity('Bs_ct_2D_ps', 'ct 2D (ps)', 'Bs_ct_2D_ps', 50, 0, 7, do_log=True),
    Quantity('Bs_ct_3D_cm', 'ct (cm)', 'Bs_ct_3D_cm', 100, 0, 0.3, do_log=True),
    Quantity('Lxy', 'L_{xy} (cm)', 'Lxy', 60, 0, 1, do_log=True),
    Quantity('Lxyz', 'L_{xyz} (cm)', 'Lxyz', 60, 0, 2, do_log=True),
    Quantity('Bs_beta', 'B_{s} #beta', 'Bs_beta', 60, 0, 1, do_log=False),
    Quantity('Bs_gamma', 'B_{s} #gamma', 'Bs_gamma', 60, 0, 30, do_log=False),
    Quantity('Bs_pt', 'B_{s} p_{T} (GeV)', 'Bs_pt', 60, 0, 50, do_log=True),
    Quantity('Bs_eta', 'B_{s} #eta', 'Bs_eta', 60, -5, 5, do_log=False),
    Quantity('Bs_phi', 'B_{s} #phi', 'Bs_phi', 60, -4, 4, do_log=False),
    Quantity('Bs_q', 'B_{s} charge', 'Bs_charge', 5, -2, 2, do_log=False),
    Quantity('Bs_pdgid', 'B_{s} pdgid', 'Bs_pdgid', 100, -800, 800, do_log=False),
    Quantity('Bs_mass', 'B_{s} mass (GeV)', 'Bs_mass', 60, 5.2, 5.5, do_log=False),
    Quantity('phi1_pt', '#phi_{1} p_{T} (GeV)', 'phi1_pt', 60, 0, 50, do_log=True),
    Quantity('phi1_eta', '#phi_{1} #eta', 'phi1_eta', 60, -5, 5, do_log=False),
    Quantity('phi1_phi', '#phi_{1} #phi', 'phi1_phi', 60, -4, 4, do_log=False),
    Quantity('phi1_q', '#phi_{1} charge', 'phi1_charge', 5, -2, 2, do_log=False),
    Quantity('phi1_pdgid', '#phi_{1} pdgid', 'phi1_pdgid', 100, -800, 800, do_log=False),
    Quantity('phi1_mass', '#phi_{1} mass (GeV)', 'phi1_mass', 60, 0.98, 1.05, do_log=False),
    Quantity('phi2_pt', '#phi_{2} p_{T} (GeV)', 'phi2_pt', 60, 0, 50, do_log=True),
    Quantity('phi2_eta', '#phi_{2} #eta', 'phi2_eta', 60, -5, 5, do_log=False),
    Quantity('phi2_phi', '#phi_{2} #phi', 'phi2_phi', 60, -4, 4, do_log=False),
    Quantity('phi2_q', '#phi_{2} charge', 'phi2_charge', 5, -2, 2, do_log=False),
    Quantity('phi2_pdgid', '#phi_{2} pdgid', 'phi2_pdgid', 100, -800, 800, do_log=False),
    Quantity('phi2_mass', '#phi_{2} mass (GeV)', 'phi2_mass', 60, 0.98, 1.05, do_log=False),
    Quantity('k1_pt', 'k_{1} p_{T} (GeV)', 'k1_pt', 60, 0, 50, do_log=True),
    Quantity('k1_eta', 'k_{1} #eta', 'k1_eta', 60, -5, 5, do_log=False),
    Quantity('k1_phi', 'k_{1} #phi', 'k1_phi', 60, -4, 4, do_log=False),
    Quantity('k1_q', 'k_{1} charge', 'k1_charge', 5, -2, 2, do_log=False),
    Quantity('k1_pdgid', 'k_{1} pdgid', 'k1_pdgid', 100, -800, 800, do_log=False),
    Quantity('k1_mass', 'k_{1} mass (GeV)', 'k1_mass', 60, 0., 1., do_log=False),
    Quantity('k2_pt', 'k_{2} p_{T} (GeV)', 'k2_pt', 60, 0, 50, do_log=True),
    Quantity('k2_eta', 'k_{2} #eta', 'k2_eta', 60, -5, 5, do_log=False),
    Quantity('k2_phi', 'k_{2} #phi', 'k2_phi', 60, -4, 4, do_log=False),
    Quantity('k2_q', 'k_{2} charge', 'k2_charge', 5, -2, 2, do_log=False),
    Quantity('k2_pdgid', 'k_{2} pdgid', 'k2_pdgid', 100, -800, 800, do_log=False),
    Quantity('k2_mass', 'k_{2} mass (GeV)', 'k2_mass', 60, 0., 1., do_log=False),
    Quantity('k3_pt', 'k_{3} p_{T} (GeV)', 'k3_pt', 60, 0, 50, do_log=True),
    Quantity('k3_eta', 'k_{3} #eta', 'k3_eta', 60, -5, 5, do_log=False),
    Quantity('k3_phi', 'k_{3} #phi', 'k3_phi', 60, -4, 4, do_log=False),
    Quantity('k3_q', 'k_{3} charge', 'k3_charge', 5, -2, 2, do_log=False),
    Quantity('k3_pdgid', 'k_{3} pdgid', 'k3_pdgid', 100, -800, 800, do_log=False),
    Quantity('k3_mass', 'k_{3} mass (GeV)', 'k3_mass', 60, 0., 1., do_log=False),
    Quantity('k4_pt', 'k_{4} p_{T} (GeV)', 'k4_pt', 60, 0, 50, do_log=True),
    Quantity('k4_eta', 'k_{4} #eta', 'k4_eta', 60, -5, 5, do_log=False),
    Quantity('k4_phi', 'k_{4} #phi', 'k4_phi', 60, -4, 4, do_log=False),
    Quantity('k4_q', 'k_{4} charge', 'k4_charge', 5, -2, 2, do_log=False),
    Quantity('k4_pdgid', 'k_{4} pdgid', 'k4_pdgid', 100, -800, 800, do_log=False),
    Quantity('k4_mass', 'k_{4} mass (GeV)', 'k4_mass', 60, 0., 1., do_log=False),
    Quantity('invmass_phi1phi2', 'm(#phi_{1}#phi_{2}) (GeV)', 'invmass_phi1phi2', 60, 5.2, 5.5, do_log=False),
    Quantity('invmass_phi1phi2', 'm(#phi_{1}#phi_{2}) (GeV)', 'invmass_phi1phi2', 60, 5.2, 5.5, do_log=True),
    Quantity('invmass_k1k2k3k4', 'm(k_{1}k_{2}k_{3}k_{4}) (GeV)', 'invmass_k1k2k3k4', 60, 5.2, 5.5, do_log=False),
    Quantity('invmass_k1k2k3k4', 'm(k_{1}k_{2}k_{3}k_{4}) (GeV)', 'invmass_k1k2k3k4', 60, 5.2, 5.5, do_log=True),
    Quantity('invmass_k1k2', 'm(k_{1}k_{2}) (GeV)', 'invmass_k1k2', 60, 0.98, 1.05, do_log=False),
    Quantity('invmass_k1k2', 'm(k_{1}k_{2}) (GeV)', 'invmass_k1k2', 60, 0.98, 1.05, do_log=True),
    Quantity('invmass_k3k4', 'm(k_{3}k_{4}) (GeV)', 'invmass_k3k4', 60, 0.98, 1.05, do_log=False),
    Quantity('invmass_k3k4', 'm(k_{3}k_{4}) (GeV)', 'invmass_k3k4', 60, 0.98, 1.05, do_log=True),
    Quantity('charge_phi1phi2', 'charge(#phi_{1}#phi_{2})', 'charge_phi1phi2', 5, -2, 2, do_log=False),
    Quantity('charge_k1k2', 'charge(k_{1}k_{2})', 'charge_k1k2', 5, -2, 2, do_log=False),
    Quantity('charge_k3k4', 'charge(k_{3}k_{4})', 'charge_k3k4', 5, -2, 2, do_log=False),
  ]

  if do_plot_mixing:
    for quantity in quantities:
      plot_mixing(inputfilename=inputfilename, outdir=outdir, quantity=quantity)

  if do_plot_pdgid:
    for quantity in quantities:
      plot_pdgid(inputfilename=inputfilename, outdir=outdir, quantity=quantity)

  if do_plot_asymmetry_mixing:
    plot_asymmetry_graph_mixing(inputfilename=inputfilename, outdir=outdir)

  if do_plot_asymmetry_pdgid:
    plot_asymmetry_graph_pdgid(inputfilename=inputfilename, outdir=outdir)

  if do_study_scale:
    for quantity in quantities:
      study_scale(quantity=quantity)


