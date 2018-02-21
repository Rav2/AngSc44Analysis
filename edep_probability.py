#!/bin/python2.7
"""
@author: Rafal Maselek
This script creates an edep histogram to estimate the probability that a given hit originates from 511 keV annihilation process.
"""
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import data_loader as dl
import numpy as np


def save_hist(X, Y, file_name="out.txt"):
    X = np.array(X[0:-1]).flatten()
    Y = np.array(Y)
    if Y.ndim == 1:
        print("Y is 1-dim")
        arr_to_save = [X, Y]
    else:
        arr_to_save = [X]
        for hist in Y:
            arr_to_save.append(np.array(hist).flatten())
        arr_to_save = np.array(arr_to_save)
    np.savetxt(fname=file_name, X=arr_to_save)


def make_histogram(file_511, file_prompt, file_name_end = "0", use_goja=False, smear=False):
    # load date
    events_true, phantom_scatt, detector_scatt, accidential = dl.load_data([file_511], [file_prompt], edep_cut=0.06, use_goja_event_analysis=use_goja)
    events = events_true + phantom_scatt + detector_scatt + accidential

    true_511 = []
    phantom_511 = []
    scintillator_511 = []
    true_prompt = []
    phantom_prompt = []
    scintillator_prompt = []
    for event in events:
        for ii in range(len(event)):
            if ii < 2:
                if event[ii].nCrystalCompton == 1:
                    if event[ii].nPhantomCompton == 0:
                        true_511.append(event[ii].edep)
                    else:
                        phantom_511.append(event[ii].edep)
                else:
                    scintillator_511.append(event[ii].edep)
            else:
                if event[ii].nCrystalCompton == 1:
                    if event[ii].nPhantomCompton == 0:
                        true_prompt.append(event[ii].edep)
                    else:
                        phantom_prompt.append(event[ii].edep)
                else:
                    scintillator_prompt.append(event[ii].edep)
    binNo = 200
    total_count_no = len(true_511) + len(phantom_511) + len(scintillator_511) + len(true_prompt) + len(phantom_prompt) + len(scintillator_prompt)
    bin_width = 1.2 / float(binNo)
    edep = [true_511, phantom_511, scintillator_511, true_prompt, phantom_prompt, scintillator_prompt]
    labels = ["511_true", "511_phantom", "511_scint", "prompt_true", "prompt_phantom", "prompt_scint"]

    # experiment-derived smearing
    new_edep = []
    if smear:
        for arr in edep:
            new_arr = []
            w = []
            for em in arr:
                new_em = np.random.normal(em, 0.044*np.sqrt(em))
                if new_em < 0.0:
                    new_em = 0.0
                new_arr.append(new_em)
            new_edep.append(new_arr)
        edep = new_edep
    # printing some info
    print("Width of single bin: {}".format(bin_width))
    print("Sum of all entries: {}".format(total_count_no))
    for ii in range(len(labels)):
        print(labels[ii]+": {}".format(len(edep[ii])))
    # making the histogram
    weights = [np.ones_like(np.array(edep_arr))/float(total_count_no) for edep_arr in edep]
    distr, bins, trash = plt.hist(edep, binNo, (0.0, 1.2), histtype='step', fill=False, stacked=False, normed=False, label=labels, weights=weights)

    print("PDF INTEGRAL = {}".format(sum([sum(x) for x in distr])))
    true_511_pdf = distr[0]
    print("PDF INTEGRAL FOR TRUE 511 KEV= {}".format(sum(true_511_pdf)))
    plt.legend(loc='upper right')
    plt.title("Deposited energy spectra, normalized as PDF")
    plt.xlabel("deposited energy [MeV]")
    plt.ylabel("probability density [1]")
    plt.savefig("results/edep_hist_norm"+file_name_end+".png")
    return distr, bins
    #plt.show()

# plt.clf()
# distr, bins = make_histogram("data/NEMA/anni50.root", "data/NEMA/prompt50.root", str(50), True, False)
# save_hist(bins, distr, "histogram.txt")
folder511 = "data/nema511_1_res/"
folder_prompt = "data/nemaprompt_1_res/"
for ii in range(1, 101, 25):
    plt.clf()
    distr, bins = make_histogram(folder511+"anni{}".format(ii)+".root", folder_prompt+"prompt{}".format(ii)+".root", str(ii), True, False)
    save_hist(bins, distr, "histogram.txt")