#!/bin/python2.7
"""
@author: Rafal Maselek
This script conducts simple binary classification for one of the two hypothesis:
1) LOR with greatest possibility to be true is true
2) LOR with greatest possibility to contain one prompt and one 511 keV is false
"""
import plotter
import data_loader as dl
import classification as cf
import numpy as np
import timeit

folder511 = "data/NEMA/"
folder_prompt = "data/NEMA/"
# lower cut on edep
edep_cut = 0.06
# True to analyze detector-scattered and accidental hits
goja_event_analysis = True
# If False, hypothesis 1) is used, otherwise 2)
use_prompt = False
# loop step sets the number of files that will be used (by default from range [1,101) )
loop_step = 8
# use sophisticated classificator
sophisticated = True
#name of the output file
filename = "lors_class_sophi_GOJA.png"

TPR = []
FPR = []
PPV = []
SPC = []

for ii in range(1, 100, loop_step):
    events_true, phantom_scatt, detector_scatt, accidential = dl.load_data([folder511+"anni{}".format(ii)+".root"],
                                                                           [folder_prompt+"prompt{}".format(ii)+".root"],
                                                                           edep_cut,
                                                                           goja_event_analysis)
    events = events_true+phantom_scatt+detector_scatt+accidential
    # name of the histogram file, remeber to use proper file for GOJA-like and non-GOJA analysis!
    if goja_event_analysis:
        histograms = np.loadtxt("histogramGOJA.txt")
    else:
        histograms = np.loadtxt("histogram.txt")
    lors, lors_from_annihilation, lors_with_prompt, lors_true_anni = dl.find_lors(events, histograms, verbose=False)
    lors_pairs = cf.remove_farthest_lor(lors)
    if sophisticated:
        tpr, spc, ppv, fpr = cf.binary_classification_probability_sophisticated(lors_pairs, verbose=False)
    else:
        tpr, spc, ppv, fpr = cf.binary_classification_probability(lors_pairs, use_prompt=use_prompt, verbose=False)
    print("[CLASSIFICATION: TPR={} SPC={} PPV={} FPR={}]".format(tpr, spc, ppv, fpr))
    TPR.append(tpr)
    SPC.append(spc)
    PPV.append(ppv)
    FPR.append(fpr)

r = range(1, 100, loop_step)
print("ALL VALUES:")
print("TPR: ",TPR)
print("PPV: ",PPV)
print("SPC: ",SPC)
print("FPR: ",FPR)
plotter.plot_classification_plots(TPR, PPV, FPR, r, filename)
