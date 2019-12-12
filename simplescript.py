#!/usr/bin/env python

import ROOT
import rootlogon
ROOT.gROOT.SetBatch(True)

gaus = ROOT.TF1('gaus', 'gaus')
gaus.SetParameters(1.0, 0.5, 0.2)

land = ROOT.TF1('land', 'landau')
land.SetParameters(1.0, 0.2, 0.1)

exp = ROOT.TF1('exp', 'exp(-3*x)')

funcs = [ gaus, land, exp ]

hists = []
for f in funcs:
    h = ROOT.TH1F('%s_hist' % f.GetName(), f.GetName(), 40, 0.0, 1.0)
    for i in xrange(1000):
        h.Fill(f.GetRandom())
    hists.append(h)

for h in hists: 
    c = ROOT.TCanvas()
    h.Draw()
    #c.SaveAs('%s.eps' % h.GetName())
    c.SaveAs('%s.pdf' % h.GetName())
    
