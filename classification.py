"""
@author: Rafal Maselek
This file contains functions for statistical classification of events/lors.
"""
#from data_loader.py import LOR
import sys
import numpy as np

def calculate_binary_coeff(TP, FP, TN, FN, verbose=False):
    """
    Calculates TPR, PPV, FPR and SPC
    :param TP: Number of true positives.
    :param FP: Number of false positives.
    :param TN: Number of true negatives.
    :param FN: Number of false negatives.
    :return: TPR, SPC, PPV, FPR
    """
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
    if verbose:
        print("TP={}; TN={}; FP={}; FN={}".format(TP,TN,FP,FN))
        print("TPR={}; SPC={}; PPV={}; FPR={}".format(TPR, SPC, PPV, FPR))
    return TPR, SPC, PPV, FPR


def binary_classification_simple(lors, d_thresholds):
    """
    Calculates TPR, PPV, FPR and SPC for many values of d_threshold
    :param lors: Array of LORs.
    :param d_thresholds: Array of threshold values of LOR distance.
    :return: Arrays TPR, FPR, PPV, SPC for all values of thresholds.
    """
    TPR = []
    FPR = []
    PPV = []
    SPC = []
    for d_t in d_thresholds:
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


def remove_farthest_lor(lors):
    grouped_lors = []
    for ii in range(len(lors)//3):
        two_lors = []
        d = []
        for jj in range(3):
            d.append(lors[3*ii+jj].d)
        max_d = max(d)
        for jj in range(3):
            if lors[3*ii+jj].d < max_d:
                two_lors.append(lors[3*ii+jj])
        if len(two_lors) == 2:
            grouped_lors.append(two_lors)
    print("[FARTHEST LORS REMOVED]")
    return np.array(grouped_lors)


def binary_classification_probability(pairs_of_lors, use_prompt=False, verbose=False):
    FP = FN = TP = TN = 0
    for pair in pairs_of_lors:
        p = []
        if not use_prompt:
            for lor in pair:
                p.append(lor.annihilation_p)
            for lor in pair:
                if lor.annihilation_p == max(p):
                    if lor.is_from_annihilation:
                        TP += 1
                    else:
                        FP += 1
                else:
                    if lor.is_from_annihilation:
                        FN += 1
                    else:
                        TN += 1
        else:
            for lor in pair:
                p.append(lor.p_prompt_and_511)
            for lor in pair:
                if lor.p_prompt_and_511 == max(p):
                    if lor.is_prompt:
                        TP += 1
                    else:
                        FP += 1
                else:
                    if lor.is_prompt:
                        FN += 1
                    else:
                        TN += 1
    tpr, spc, ppv, fpr = calculate_binary_coeff(TP, FP, TN, FN, verbose)
    print("[BINARY CLASSIFICATION DONE]")
    return tpr, spc, ppv, fpr


def binary_classification_probability_sophisticated(pairs_of_lors, verbose=False):
    FP = FN = TP = TN = 0
    # TODO:FINISH
    for pair in pairs_of_lors:
        p = []
        for ii in range(2):
            lor = pair[ii]
            p.append(lor.annihilation_p * (1.0-lor.p_prompt_and_511))
            if pair[ii].d < pair[(ii+1)%2]:
                p *= 0.9
            else:
                p *= 0.1
        ii_max_p = 3
        if p[0] > p[1]:
            ii_max_p = 0
        else:
            ii_max_p = 1

        for ii in range(2):
            if ii == ii_max_p:
                if lor.is_from_annihilation:
                    TP += 1
                else:
                    FP += 1
            else:
                if lor.is_from_annihilation:
                    FN += 1
                else:
                    TN += 1

    tpr, spc, ppv, fpr = calculate_binary_coeff(TP, FP, TN, FN, verbose)
    print("[BINARY CLASSIFICATION DONE]")
    return tpr, spc, ppv, fpr
