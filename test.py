import timeit
setup_str = """
import numpy as np
import plotter
import data_loader as dl
import classification as cf


folder511 = "data/dummy/" #nema511_1_res/"
folder_prompt = "data/dummy/" #nemaprompt_1_res/"
edep_cut = 0.06
goja_event_analysis = False
use_prompt = False
TPR = []
FPR = []
PPV = []
SPC = []
ii = 1
"""
snippet="""
dl.load_data([folder511+"anni{}".format(ii)+".root"],
	                    [folder_prompt+"prompt{}".format(ii)+".root"],
	                    edep_cut,
	                    goja_event_analysis)

"""

print timeit.timeit(setup = setup_str,
                    stmt = snippet,
                    number = 1)