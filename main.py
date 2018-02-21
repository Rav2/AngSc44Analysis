#!/bin/python2.7
"""
@author: Rafal Maselek
Main file of the scipt package
"""
from __future__ import print_function
import plotter
import data_loader as dl
import classification as cf
import numpy as np
import math


def event_check(events):
    """
    Checks if distance between two points is larger than 400 mm
    :param events: List of Events
    :return: nothing (prints INFO if distance is lower than 400 mm)
    """
    counter = 0
    wrong_events = []
    for event in events:
        x1 = event[0].posX
        y1 = event[0].posY
        x2 = event[1].posX
        y2 = event[1].posY
        d = math.sqrt((x1-x2)**2+(y1-y2)**2)
        if d < 400:
            print('distance between two points caused by 511keV gammas is too small!')
            counter += 1
            print(d)
            print(event[0].allFields)
            print(event[1].allFields)
            print()
            wrong_events.append(event)
    print('suspicious events=', counter)
    print('all events', len(events))
    plotter.plot_hits(wrong_events)


def prepare_fname_lists(annihilation_fname, prompt_fname, folder, N):
    """
    Creates list of file names numerated by integers. 
    :param annihilation_fname: Core name of the file with 511 keV data.
    :param prompt_fname: Core name of the file with prompt data.
    :param folder: Path to the folder containing files.
    :param N: Maximal number of files, used to enumerate them. 
    :return: Two lists: list of file names for 511 keV/prompt data
    """
    file_list_511 =[]
    file_list_prompt = []
    for ii in range(1, N+1):
        file_list_511.append(folder+"/"+annihilation_fname+"{}".format(ii)+".root")
        file_list_prompt.append(folder+"/"+prompt_fname+"{}".format(ii)+".root")
    return file_list_511, file_list_prompt

def main():
    """
    Main function of the program.
    :return: nothing
    """
    print('[START]')
    file_list_511 = ['anni50.root', ]
    file_list_prompt = ['prompt50.root']
    data_folder_511 = "data/nema511_1_res"
    data_folder_prompt = "data/nemaprompt_1_res"
    edep_cut = 0.06
    short_run = True # Set True if using only one file for each type of data. Otherwise set False and data from all files will be loaded.
    goja_event_analysis = True
    file_list_511 = [data_folder_511+'/'+name for name in file_list_511]
    file_list_prompt = [data_folder_prompt+'/'+name for name in file_list_prompt]

    if(short_run):
        events_true, phantom_scatt, detector_scatt, accidential = dl.load_data([file_list_511[-1]], [file_list_prompt[-1]], edep_cut, goja_event_analysis)
    else:
        events_true, phantom_scatt, detector_scatt, accidential = dl.load_data(file_list_511, file_list_prompt, edep_cut, goja_event_analysis)
    # event_check(events)

    events = events_true+phantom_scatt+detector_scatt+accidential

    # dl.write_goja_output(events_true+phantom_scatt+detector_scatt+accidential)
    histograms = np.loadtxt("histogram.txt")
    lors, lors_from_annihilation, lors_with_prompt, lors_true_anni= dl.find_lors(events, histograms)
    # lors_pairs = cf.remove_farthest_lor(lors)
    # tpr, spc, ppv, fpr = cf.binary_classification_probability(lors_pairs, use_prompt=False)
    # print("511 kev: tpr={} spc={} ppv={} fpr={}".format(tpr,spc,ppv,fpr))

    plotter.plot_edep_distribution(events, filename='edep_distribution')
    plotter.plot_position_and_time_distribution(events, only_511keV = True)
    plotter.plot_d_distribution(lors, filename='d_distribution')
    plotter.plot_d_distribution(lors_from_annihilation, filename='d_distribution_annihilation_lors')
    plotter.plot_d_distribution(lors_true_anni, filename='d_distribution_true_annihilation_lors')
    plotter.plot_d_distribution(lors_with_prompt, filename='d_distribution_lors_with_prompt')

    # d_tresholds = np.linspace(0, 437.3, 101)
    # TPR, SPC, PPV, FPR = cf.binary_classification_simple(lors, d_tresholds)
    # if len(TPR):
        # plotter.plot_classification_plots(TPR, PPV, FPR, d_tresholds)
    print('[EXIT]')


# execute main only if not imported to other script
if __name__ == "__main__":
    main()
