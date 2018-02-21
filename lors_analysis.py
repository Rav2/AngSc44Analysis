#!/bin/python2.7
"""
@author: Rafal Maselek
This script analysis lors: it selects true lors and organizes them by their distance to the origin, relative amongst all
three lors
"""
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import data_loader as dl
import plotter
import classification as cf
import numpy as np
import math
from main import *


def analyse_lors(files_511, files_prompt, goja_event_analysis, loop_step):
    """
    Counts the number of true lors that have the smallest distance from the origin out of three, or the highest distance, 
    or not the smallest nor the highest (middle one).
    :param files_511: List of strings with file names of files with 511 keV data or single string in case of a single file.
    :param files_prompt: List of strings with file names of files with prompt data or single string in case of a single file.
    :param goja_event_analysis: If true, GOJA-like analysis will be performed.
    :param loop_step: Loop step. E.g. if it is equal to 5, one out of 5 files will be analyzed.
    :return: Arrays d_min, d_mid, d_max containing numbers of lors in given category. Each entry corresponds to different
    file.
    """
    d_min = []
    d_mid = []
    d_max = []
    for nn in range(0, nloops, loop_step):
        if short_run:
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

        lors, lors_from_annihilation, lors_with_prompt, lors_true_anni = dl.find_lors(events)
        d1, d2, d3 = dl.count_sorted_lors(lors)
        d_min.append(d1)
        d_mid.append(d2)
        d_max.append(d3)
    return d_min, d_mid, d_max

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# SCRIPT PARAMETERS #
print('[START]')
edep_cut = 0.06
short_run = False 
goja_event_analysis = False
nloops = 100
data_folder = "data/NEMA"
loop_step = 10

nloops2 = 100
data_folder2 = "data/NEMA"
use_second_data = True # if yes, then second set of data will be analyzed and plotted on the same plot
loop_step2 = 10
# LOADING DATA #
file_list_511, file_list_prompt = prepare_fname_lists("anni", "prompt", data_folder, nloops)
dmin1, dmid1, dmax1 = analyse_lors(file_list_511, file_list_prompt, goja_event_analysis, loop_step)
# Loading and analyzing second set of data
if use_second_data:
    file_list_511, file_list_prompt = prepare_fname_lists("anni", "prompt", data_folder2, nloops2)
    dmin2, dmid2, dmax2 = analyse_lors(file_list_511, file_list_prompt, goja_event_analysis, loop_step2)

# DRAWING PLOT(S) #
if nloops > 1:
    r =range(1, nloops+1, loop_step)
    if use_second_data:
        r2 = range(1, nloops + 1, loop_step2)
        plotter.plot_lors_fractions2(dmin1, dmid1, dmax1, r, dmin2, dmid2, dmax2, r2)
    else:
        plotter.plot_lors_fractions(dmin1, dmid1, dmax1, r, filename="fractions_of_lors")



