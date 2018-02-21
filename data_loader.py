"""
@author: Rafal Maselek
This file contains functions that load and sort data.
It contains definitions of Lor and Hit classes.
"""
from ROOT import gROOT, TCanvas, TH1, TH2, TTree, TFile
import math
import sys
import numpy as np
from enum import Enum
from matplotlib import pyplot as plt


class CoincType(Enum):
    """
    Enum type for coincidences (like in P. Kowalski's program 'GOJA')
    """
    kUnspecified = 0
    kTrue = 1
    kPhantomScattered = 2
    kDetectorScattered = 3
    kAccidental = 4


class Hit:
    """
    Class representing a signle hit, like in GATE output.
    """
    def __init__(self, tree):
        self.PDGEncoding = tree.PDGEncoding
        self.trackID = tree.trackID
        self.parentID = tree.parentID
        self.time = tree.time
        self.edep = tree.edep
        self.posX = tree.posX
        self.posY = tree.posY
        self.posZ = tree.posZ
        self.baseID = tree.baseID
        self.photonID = tree.photonID
        self.nPhantomCompton = tree.nPhantomCompton
        self.nCrystalCompton = tree.nCrystalCompton
        self.nPhantomRayleigh = tree.nPhantomRayleigh
        self.nCrystalRayleigh = tree.nCrystalRayleigh
        self.primaryID = tree.primaryID
        self.sourcePosX = tree.sourcePosX
        self.sourcePosY = tree.sourcePosY
        self.sourcePosZ = tree.sourcePosZ
        self.sourceID = tree.sourceID
        self.eventID = tree.eventID
        self.volumeID = tree.volumeID
        self.processName = tree.processName
        self.comptVolName = tree.comptVolName
        self.RayleighVolName = tree.RayleighVolName
        self.coincType = CoincType.kUnspecified
        self.allFields = [self.PDGEncoding, self.trackID, self.parentID, self.time,self.edep, self.posX, self.posY,
                            self.posZ, self.nPhantomCompton, self.baseID, self.photonID, 
                            self.nCrystalCompton, self.nPhantomRayleigh, self.nCrystalRayleigh, self.primaryID,
                            self.sourcePosX, self.sourcePosY, self.sourcePosZ,
                            self.sourceID, self.eventID, self.volumeID,
                            self.processName, self.comptVolName, self.RayleighVolName]


class LOR:
    """Class represnting a single Line of response in 2D."""
    def __init__(self, d, is_from_annihilation, is_prompt, annihilation_p, prompt_p):
        self.d = d # distance from the origin.
        self.is_from_annihilation = is_from_annihilation # true if LOR connects two hits by un-scattered 511 keV gammas
        self.is_prompt = is_prompt
        self.annihilation_p = annihilation_p
        self.prompt_p = prompt_p


def is_proper_hit(hit, edep_cut, use_goja_event_analysis=False):
    """
    Checks if hit passes some basic selections.
    :param hit: Hit object.
    :param edep_cut: Value of lower cut on deposited energy.
    :param use_goja_event_analysis: 
    :return: True if hit passes all selections.
    """
    val = hit.nCrystalRayleigh==0 and hit.nPhantomRayleigh==0 \
        and hit.edep >= edep_cut and hit.PDGEncoding == 22 \
        and (hit.processName[:-1].lower().strip() == 'compton' or hit.processName[:-1].lower().strip() ==  'compt')
    # if use_goja_event_analysis==False then check additional constraint
    if not use_goja_event_analysis:
        return val and hit.nCrystalCompton==1 #and hit.nPhantomCompton==0
    else:
        return val


def goja_event_analysis(hits):
    """
    Determine the coincidence type of two hits originating from two annihilation gammas.
    :param hits: Array or tuple of at least two hits, where two first come from e-e+ annihilation.
    :return: Array of two hits with coincidence types set.
    """
    h1 = hits[0]
    h2 = hits[1]

    if h1.eventID == h2.eventID:
        if h1.nPhantomCompton==0 and h2.nPhantomCompton==0:
            if h1.nCrystalCompton==1 and h2.nCrystalCompton==1:
                t = CoincType.kTrue
            else:
                 t = CoincType.kDetectorScattered
        else:
            t = CoincType.kPhantomScattered
    else:
        t = CoincType.kAccidental # this might be used with time window separation, with eventID separation it's useless
    h1.coincType = t
    h2.coincType = t
    return [h1, h2]


def find_coincidences(hits, edep_cut, use_goja_event_analysis=False):
    """
    Checks if hits with the same eventID form a valid 2-gamma event.
    :param hits: Array of hits with the same eventID.
    :param edep_cut: Lower selection cut on deposited energy.
    :param use_goja_event_analysis: If true, only hits that were scattered only once are proper. Also allows an event to
    contain more than 2 proper hits.
    :return: Array with hits forming a coincidence.
    """
    new_hits = []
    for hit in hits:
        if is_proper_hit(hit, edep_cut, use_goja_event_analysis):
            new_hits.append(hit)
    if len(new_hits)>2 and not use_goja_event_analysis:
        print('event with more than 2 proper gammas, not in goja mode')
    elif len(new_hits) == 2:
        new_hits = goja_event_analysis(new_hits)
        return new_hits

    return []


def load_data(file_list_511, file_list_prompt, edep_cut = 0.06, use_goja_event_analysis=False):
    """
    Loads data from files with 511 keV and prompt data.
    :param file_list_511: List of files of 511 keV data.
    :param file_list_prompt: List of files of prompt data.
    :param edep_cut: Lower selection value for deposited energy.
    :param use_goja_event_analysis: If true, GOJA-like analysis is performed.
    :return: Four lists of events: true, phantom-scattered, detector-scattered, accidental.
    """

    print('\n[LOADING 511 KEV DATA...]')
    proper_hits_511 = []
    proper_hits_prompt = []
    add_to_proper_hits_511 = proper_hits_511.append
    add_to_proper_hits_prompt = proper_hits_prompt.append

    # loading of 511 keV data
    for f511_file_name in file_list_511:
        print("[LOADING: "+f511_file_name+"]")
        f_511 = TFile(f511_file_name)
        bufor = []
        # TODO: Optimize reading from the tree, because it takes way to much time
        for event in f_511.Hits:
            if len(bufor) == 0 or bufor[0].eventID == event.eventID:
                bufor.append(Hit(event))
            else:
                for proper_hit in find_coincidences(bufor, edep_cut, use_goja_event_analysis):
                    add_to_proper_hits_511(proper_hit)
                bufor = [Hit(event)]
        f_511.Close()


    print('[511 KEV DATA LOADED. LOADING PROMPT DATA...]')
    # loading prompt data
    for f_prompt_file_name in file_list_prompt:
        print("[LOADING: "+f_prompt_file_name+"]")
        f_prompt = TFile(f_prompt_file_name)
        for event in f_prompt.Hits :
            if is_proper_hit(event, edep_cut, use_goja_event_analysis):
                add_to_proper_hits_prompt(Hit(event))
            if 2*len(proper_hits_prompt) > len(proper_hits_511):
                break
        f_prompt.Close()


    # truncating the data to equal number of events      
    n_511 = len(proper_hits_511)
    n_prompt = len(proper_hits_prompt)
    if n_511//2 > n_prompt:
        n_events = n_prompt
    else:
        n_events = n_511//2

    print('[NO OF EVENTS: {}]'.format(n_events))

    events_true = []
    events_phantom_scattered = []
    events_detector_scattered = []
    events_accidental = []
    add_to_true = events_true.append
    add_to_phantom_scattered = events_phantom_scattered.append
    add_to_detector_scattered = events_detector_scattered.append
    add_to_accidental = events_accidental.append
    # classifying events based on the coincidence type of the two annihilation gammas
    for ii in range(n_events):
        ev = (proper_hits_511[2*ii], proper_hits_511[2*ii+1], proper_hits_prompt[ii])
        t = proper_hits_511[2*ii].coincType
        if t == CoincType.kTrue:
            add_to_true(ev)
        elif t == CoincType.kPhantomScattered:
            add_to_phantom_scattered(ev)
        elif t == CoincType.kDetectorScattered:
            add_to_detector_scattered(ev)
        elif t == CoincType.kAccidental:
            add_to_accidental(ev)
    print('[NO OF TRUE EVENTS: {}, NO OF PHANTOM-SCATTERED EVENTS: {}, NO OF DETECTOR-SCATTERED EVENTS: {}, NO OF ACCIDENTAL EVENTS: {}]'\
        .format(len(events_true), len(events_phantom_scattered), len(events_detector_scattered), len(events_accidental)))
    print('[DATA LOADED]')
    return events_true, events_phantom_scattered, events_detector_scattered, events_accidental


def find_lors(events, histograms = [], verbose=False):
    """
    Finds LOR parameters: distance from the origin and angle between OX axis and vector from the origin to the center of LOR.
    Projections of 3D LORs onto OXY plane are analysed.
    :param events: List of events.
    :return: Lists of LOR objects: all LORs, LORs coming from annihilation that weren't scattered, LORs containg prompt hit.
    """
    lors = []
    lors_true_annihilation = []
    lors_annihilation = []
    lors_with_prompt = []
    probs511 = []
    probsPrompt = []
    for event in events:
        for ii in range(3):
            edep1 = event[ii].edep
            edep2 = event[(ii+1)%3].edep
            p1 = np.array([event[ii].posX, event[ii].posY, event[ii].posZ])
            p2 = np.array([event[(ii+1)%3].posX, event[(ii+1)%3].posY, event[(ii+1)%3].posZ])
            # construct vector representing the direction of LOR
            l = p1-p2
            # find s, which idicates intersection point on LOR with the normal that passes through the origin
            s = -(p1[0]*l[0]+p1[1]*l[1]+p1[2]*l[2])/(l[0]**2+l[1]**2+l[2]**2)
            p_intersection = p1 + s*l
            # d is the distance between the origin (point (0,0,0)) and LOR (intersection point)
            d = np.sqrt(p_intersection[0]**2+p_intersection[1]**2+p_intersection[2]**2)

            p511 = 1.0
            p_prompt = 1.0
            if len(histograms) > 0:
                # find the probability that LOR is true
                # find the bin for given edep1
                index = 0
                while histograms[0][index] < edep1:
                    index += 1
                p511 *= histograms[1][index-1]
                 # find the bin for given edep2
                index = 0
                while histograms[0][index] < edep2:
                    index += 1
                p511 *= histograms[1][index-1]
                # find the probability that LOR contains prompt
                while histograms[0][index] < edep1:
                    index += 1
                p_prompt *= histograms[4][index-1]
                index = 0
                while histograms[0][index] < edep2:
                    index += 1
                p_prompt *= histograms[4][index-1]

            # for LORs with at least one annihilation hit
            if ii == 0:
                is_true_annihilation = (event[ii].coincType == CoincType.kTrue)
                lors.append(LOR(d, is_true_annihilation, False, p511, p_prompt)) # True for back-to-back emission
                lors_annihilation.append(LOR(d, is_true_annihilation, False, p511, p_prompt))
                if is_true_annihilation:
                    lors_true_annihilation.append(LOR(d, is_true_annihilation, False, p511, p_prompt))
            # for LORs with prompt hit
            else:
                lors.append(LOR(d, False, True, p511, p_prompt))
                lors_with_prompt.append(LOR(d, False, True, p511, p_prompt))
            probs511.append(p511)
            probsPrompt.append(p_prompt)
    if verbose:
        print('[Verbose mode is ON. Histogram of LORs probabilities will be saved to a file.]')
        plt.hist([probs511, probsPrompt], 100, histtype='step', fill=False, stacked=False, normed=False, label=["511", "prompt"])
        plt.legend(loc='upper right')
        plt.savefig("p_plots.png")
    print('[LORS FOUND]')
    return lors, lors_annihilation, lors_with_prompt, lors_true_annihilation


def count_sorted_lors(lors):
    """
    True lors are classified by their relative distance to the origin: the smallest of three lors, the middle one, or the greatest.
    :param lors: List of LOR objects grouped to events.
    :return: Fraction of LORs that are: closest, middle, farthest from the origin in the whole data sample.
    """
    d_min = 0
    d_mid = 0
    d_max = 0
    print(len(lors), len(lors)/3)
    for ii in range(len(lors)/3):
        index = -10
        for jj in range(3):
            if lors[3*ii+jj].is_from_annihilation:
                index = jj
                break
        if index == -10:
            continue
        if lors[3*ii+index].d < lors[3*ii+(index+1)%3].d and lors[3*ii+index].d < lors[3*ii+(index+2)%3].d:
            d_min += 1
        elif lors[3*ii+index].d < lors[3*ii+(index+1)%3].d or lors[3*ii+index].d < lors[3*ii+(index+2)%3].d:
            d_mid += 1
        else:
            d_max += 1
    annihilation_lors_no = float(d_min+d_mid+d_max)
    # avoid division by 0 error:
    if annihilation_lors_no == 0:
        return 0,0,0
    # divide by the number of aniihilation lors to get fractions
    print("[ALL LORS: {}]".format(len(lors)))
    print("[ANNIHILATION LORS FOUND: {}]".format(annihilation_lors_no))
    return d_min/annihilation_lors_no, d_mid/annihilation_lors_no, d_max/annihilation_lors_no


def write_goja_output(coincidences, filename='goja_output.txt'):
    """
    Writes GOJA-like output to a text file.
    :param coincidences: List of coincidences (pairs of Hits).
    :param filename: Name of the output file.
    :return: nothing
    """
    with open(filename, 'w') as f:
        for coinc in coincidences:
            if len(coinc) < 2:
                raise Exception("A coincidence provided with less than two hits!")
            h1 = coinc[0]
            h2 = coinc[1]
            line = np.array([h1.posX/10.0, h1.posY/10.0, h1.posZ/10.0, h1.time*10**12, \
                            h2.posX/10.0, h2.posY/10.0, h2.posZ/10.0, h2.time*10**12,\
                        h1.volumeID[1], h2.volumeID[1], h1.edep*1000, h2.edep*1000,\
                         h1.coincType.value, h1.sourcePosX/10.0, h1.sourcePosY/10.0, h1.sourcePosZ/10.0])
            fmt = '%.2f\t%.2f\t%.2f\t%.1f\t%.2f\t%.2f\t%.2f\t%.1f\t%.1f\t%.1f\t%.2f\t%.2f\t%d\t%.2f\t%.2f\t%.2f'
            np.savetxt(f, line.reshape(1, line.shape[0]), fmt)
