#!/usr/bin/env python

from ROOT import *
import sys, os
sys.path.append("../python")
from NtupleAnalyzer import *

ans = {}
ans['ttll'] = NtupleAnalyzer(
    inFileNames = "../../ProcessEvent/test/f.root",
    outFileName = "ttll.root",
    modName = "ntuple", treeName = "event", eventCounterName = "hNEvent",
)

for an in ans:
    an = ans[an]
    an.addH1("l1_pt", "leptons_pt[0]", "Leading lepton p_{T};p_{T}^{l, 1st lead} [GeV];Events/20 GeV", 100, 0, 200)
    an.addH1("l2_pt", "leptons_pt[1]", "Trailing lepton p_{T};p_{T}^{l, 2nd lead} [GeV];Events/20 GeV", 100, 0, 200)

    an.addCutStep("S1", "", "l1_pt,l2_pt,l1_eta,l2_eta")
    an.addCutStep("S2", "abs(mLL-91.2) > 15", "l1_pt,l2_pt,l1_eta,l2_eta")

an.process()
