import data_loader as dl
import classification as cf
import plotter
import numpy as np

folder511 = "data/nema511_1_res/"
folder_prompt = "data/nemaprompt_1_res/"
edep_cut = 0.06
goja_event_analysis = False
use_prompt = False
loop_step = 15


TPR = []
FPR = []
PPV = []
SPC = []

for ii in range(1, 101, loop_step):
    events_true, phantom_scatt, detector_scatt, accidential = dl.load_data([folder511+"anni{}".format(ii)+".root"],
                                                                           [folder_prompt+"prompt{}".format(ii)+".root"],
                                                                           edep_cut,
                                                                           goja_event_analysis)
    events = events_true+phantom_scatt+detector_scatt+accidential
    histograms = np.loadtxt("histogram.txt")
    lors, lors_from_annihilation, lors_with_prompt = dl.find_lors(events, histograms, True)
    lors_pairs = cf.remove_farthest_lor(lors)
    tpr, spc, ppv, fpr = cf.binary_classification_probability(lors_pairs, use_prompt=use_prompt, verbose=False)
    if not use_prompt:
        print("511 kev file: {}; tpr={} spc={} ppv={} fpr={}".format("anni{}".format(ii), tpr, spc, ppv, fpr))
    TPR.append(tpr)
    SPC.append(spc)
    PPV.append(ppv)
    FPR.append(fpr)

r = range(1, 101, loop_step)
plotter.plot_classification_plots(TPR, PPV, FPR, r, "lors_classification.png")
