#!/bin/python2.7
from __future__ import print_function
import data_loader as dl
import plotter
import classification as cf
import numpy as np
import math


def event_check(events):
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
    file_list_511 =[]
    file_list_prompt = []
    for ii in range(1, N+1):
        file_list_511.append(folder+"/"+annihilation_fname+"{}".format(ii)+".root")
        file_list_prompt.append(folder+"/"+prompt_fname+"{}".format(ii)+".root")
    return file_list_511, file_list_prompt

def main():
    print('[START]')
    file_list_511 = ['511keV_1.root', ]
    file_list_prompt = ['prompt_1.root']
    data_folder = "old_data"
    edep_cut = 0.06
    short_run = False
    goja_event_analysis = False
    file_list_511 = [data_folder+'/'+name for name in file_list_511]
    file_list_prompt = [data_folder+'/'+name for name in file_list_prompt]

    if(short_run):
        events_true, phantom_scatt, detector_scatt, accidential = dl.load_data([file_list_511[-1]], [file_list_prompt[-1]], edep_cut, goja_event_analysis)
    else:
        events_true, phantom_scatt, detector_scatt, accidential = dl.load_data(file_list_511, file_list_prompt, edep_cut, goja_event_analysis)
    # event_check(events)

    events = events_true+phantom_scatt
    print(len(events_true), ' ', len(phantom_scatt), ' ', len(events))
    dl.write_goja_output(events_true+phantom_scatt+detector_scatt+accidential)
    lors, lors_from_annihilation, lors_with_prompt = dl.find_lors(events)

    plotter.plot_edep_distribution(events, filename='edep_distribution')
    plotter.plot_position_and_time_distribution(events, only_511keV = True)
    plotter.plot_d_distribution(lors, filename='d_distribution')
    plotter.plot_d_distribution(lors_from_annihilation, filename='d_distribution_annihilation_lors')
    plotter.plot_d_distribution(lors_with_prompt, filename='d_distribution_lors_with_prompt')

    d_tresholds = np.linspace(0, 437.3, 101)
    TPR, SPC, PPV, FPR = cf.binary_classification_simple(lors, d_tresholds)
    if(len(TPR)):
        plotter.plot_classification_plots(TPR, PPV, FPR, d_tresholds)
    print('[EXIT]')


if __name__ == "__main__":
    main()
