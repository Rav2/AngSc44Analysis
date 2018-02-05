#!/bin/python2.7
from __future__ import print_function
import data_loader as dl
import plotter
import classification as cf
import numpy as np
import math
from main import *


print('[START]')
file_list_511 = []
file_list_prompt = []
edep_cut = 0.06
short_run = False 
goja_event_analysis = False
nloops = 100
data_folder = "data/NEMA"
file_list_511 = [data_folder+'/'+name for name in file_list_511]
file_list_prompt = [data_folder+'/'+name for name in file_list_prompt]

if nloops > 1:
	file_list_511, file_list_prompt = prepare_fname_lists("anni", "prompt", data_folder, 100)
	d_min = []
	d_mid = []
	d_max = []


loop_step = 1
for nn in range(0, nloops, loop_step):
	if(short_run):
		events_true, phantom_scatt, detector_scatt, accidential = dl.load_data([file_list_511[-1]], [file_list_prompt[-1]], edep_cut, goja_event_analysis)
	elif nloops>1:
		events_true, phantom_scatt, detector_scatt, accidential = dl.load_data([file_list_511[nn]], [file_list_prompt[nn]], edep_cut, goja_event_analysis)
	else:
		events_true, phantom_scatt, detector_scatt, accidential = dl.load_data(file_list_511, file_list_prompt, edep_cut, goja_event_analysis)
	# event_check(events)

	events = events_true+phantom_scatt
	print(len(events_true), ' ', len(phantom_scatt), ' ', len(events))

	lors, lors_from_annihilation, lors_with_prompt = dl.find_lors(events)
	d1, d2, d3 = dl.count_sorted_lors(lors)
	d_min.append(d1)
	d_mid.append(d2)
	d_max.append(d3)

if nloops > 1:
	r =range(1, nloops+1, loop_step)
	plotter.plot_lors_fractions(d_min, d_mid, d_max, r, "fractions_of_lors")

