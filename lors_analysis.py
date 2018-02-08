#!/bin/python2.7
from __future__ import print_function
import matplotlib.pyplot as plt
import data_loader as dl
import plotter
import classification as cf
import numpy as np
import math
from main import *


def analyse_lors(files_511, files_prompt, goja_event_analysis, loop_step):
    d_min = []
    d_mid = []
    d_max = []
    loop_step = 10
    for nn in range(0, nloops, loop_step):
        if (short_run):
            events_true, phantom_scatt, detector_scatt, accidential = dl.load_data([files_511[-1]],
                                                                                   [files_prompt[-1]], edep_cut,
                                                                                   goja_event_analysis)
        elif nloops > 1:
            events_true, phantom_scatt, detector_scatt, accidential = dl.load_data([files_511[nn]],
                                                                                   [files_prompt[nn]], edep_cut,
                                                                                   goja_event_analysis)
        else:
            events_true, phantom_scatt, detector_scatt, accidential = dl.load_data(files_511, files_prompt,
                                                                                   edep_cut, goja_event_analysis)
        # event_check(events)

        events = events_true + phantom_scatt
        print(len(events_true), ' ', len(phantom_scatt), ' ', len(events))

        lors, lors_from_annihilation, lors_with_prompt = dl.find_lors(events)
        d1, d2, d3 = dl.count_sorted_lors(lors)
        d_min.append(d1)
        d_mid.append(d2)
        d_max.append(d3)
    return d_min, d_mid, d_max


print('[START]')
edep_cut = 0.06
short_run = False 
goja_event_analysis = False
nloops = 100
data_folder = "data/NEMA"
loop_step = 10

nloops2 = 100
data_folder2 = "data/NEMA"
use_second_data = True
loop_step2 = 10

file_list_511, file_list_prompt = prepare_fname_lists("anni", "prompt", data_folder, nloops)
dmin1, dmid1, dmax1 = analyse_lors(file_list_511, file_list_prompt, goja_event_analysis, loop_step)
if use_second_data:
    file_list_511, file_list_prompt = prepare_fname_lists("anni", "prompt", data_folder2, nloops2)
    dmin2, dmid2, dmax2 = analyse_lors(file_list_511, file_list_prompt, goja_event_analysis, loop_step2)

if nloops > 1:
    r =range(1, nloops+1, loop_step)
    f = plt.figure()
    ax1 = plotter.plot_lors_fractions(dmin1, dmid1, dmax1, r, filename="fractions_of_lors")
    if use_second_data:
        r2 = range(1, nloops + 1, loop_step2)
        ax2 = plotter.plot_lors_fractions(dmin2, dmid2, dmax2, r2, labels=["d_min2", "d_mid2", "d_max2"], style="-c-m-y",\
                                          filename="fractions_of_lors")
        f.axes.append(ax1)
        f.axes.append(ax2)
        plt.savefig("TWO_PLOTS.png")


