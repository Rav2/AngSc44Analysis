#from data_loader.py import LOR
import sys

def calculate_binary_coeff(TP, FP, TN, FN):
	if TP+FN == 0:
		TPR = 1
	else:
		TPR =float(TP)/float(TP+FN)
	if TN+FP == 0:
		SPC = 1
	else:
		SPC = TN/float(TN+FP)
	if TP+FP == 0:
		PPV = 1
	else:
		PPV = TP/float(TP+FP)
	if TN+FP == 0:
		FPR = 1
	else:
		FPR = 1.0-TN/float(TN+FP)
	return TPR, SPC, PPV, FPR



def binary_classification_simple(lors, d_tresholds):
    #calculation of efficiency, purity etc.
	TPR = []
	FPR = []
	PPV = []
	SPC = []
	for d_t in d_tresholds:
		FP=FN=TP=TN=0
		for lor in lors:
			if lor.d > d_t:
				if lor.is_from_annihilation:
					FN += 1
				else:
					TN += 1
			else:  
				if lor.is_from_annihilation:
					TP += 1
				else:
					FP += 1
		tpr, spc, ppv, fpr = calculate_binary_coeff(TP, FP, TN, FN)
		TPR.append(tpr)
		FPR.append(fpr)
		PPV.append(ppv)
		SPC.append(spc)
	print("[BINARY CLASSIFICATION DONE]")
	return TPR, FPR, PPV, SPC

def binary_table(events, lors, edep_threshold, dist_threshold):
	if 3*len(events) != len(lors):
		raise Exception("Number of lors doesn't match the number of events!")
	# hypothesis: lor is true iff conditions from both classifications are fullfiled
	FP=FN=TP=TN=0
	proper_events = []
	proper_anni_lors = []
	for ii in range(len(events)):
		for jj in range(3):
			cond_edep = events[ii][jj].edep < edep_threshold and events[ii][(jj+1)%3].edep < edep_threshold
			cond_dist = lors[3*ii+jj].d < dist_threshold
			if cond_edep and cond_dist:
				if lors[3*ii+jj].is_from_annihilation:
					TP += 1
					proper_anni_lors.append(lors[3*ii+jj])
					proper_events.append(events[ii])
				else:
					FP += 1
			elif lors[3*ii+jj].is_from_annihilation:
				FN += 1
			else:
				TN += 1
	# TODO: finish it
	return calculate_binary_coeff(TP, FP, TN, FN), proper_events, proper_anni_lors