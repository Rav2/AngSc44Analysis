from ROOT import gROOT, TCanvas, TH1, TH2, TTree, TFile
import math
import sys
import numpy as np
from enum import Enum

class CoincType(Enum):
	kUnspecified = 0
	kTrue = 1
	kPhantomScattered = 2
	kDetectorScattered = 3
	kAccidental = 4


class Hit:
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
	def __init__(self, d, theta, is_from_annihilation):
		self.d = d
		self.theta = theta
		self.is_from_annihilation = is_from_annihilation


def is_proper_hit(hit, edep_cut, use_goja_event_analysis=False):
	val = hit.nCrystalRayleigh==0 and hit.nPhantomRayleigh==0 \
		and hit.edep >= edep_cut and hit.PDGEncoding == 22 \
		and (hit.processName[:-1].lower().strip() == 'compton' or hit.processName[:-1].lower().strip() ==  'compt')
	if not use_goja_event_analysis:
		return val and hit.nCrystalCompton==1 #and hit.nPhantomCompton==0
	else:
		return val


def goja_event_analysis(hits):
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
	print('[LOADING 511 KEV DATA...]')
	proper_hits_511 = []
	proper_hits_prompt = []
	add_to_proper_hits_511 = proper_hits_511.append
	add_to_proper_hits_prompt = proper_hits_prompt.append

	for f511_file_name in file_list_511:
		print(f511_file_name)
		f_511 = TFile(f511_file_name)
		bufor = []
		for event in f_511.Hits:
			if len(bufor) == 0 or bufor[0].eventID == event.eventID:
				bufor.append(Hit(event))
			else:
				for proper_hit in find_coincidences(bufor, edep_cut, use_goja_event_analysis):
					add_to_proper_hits_511(proper_hit)
				bufor = [Hit(event)]
		#f_511.Close()


	print('[511 KEV DATA LOADED. LOADING PROMPT DATA...]')
	for f_prompt_file_name in file_list_prompt:
		f_prompt = TFile(f_prompt_file_name)
		for event in f_prompt.Hits :
			if is_proper_hit(event, edep_cut, False):
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
	print('[NO OF TRUE EVENTS: {}, NO OF PHANTOM-SCATTERED EVENTS: {}, NO OF DETECTOR-SCATTERED EVENTS: {}, NO OF ACCIDENTIAL EVENTS: {}]'\
		.format(len(events_true), len(events_phantom_scattered), len(events_detector_scattered), len(events_accidental)))
	print('[DATA LOADED]')
	return events_true, events_phantom_scattered, events_detector_scattered, events_accidental


def find_lors(events):
	lors = []
	lors_annihilation = []
	lors_with_prompt = []
	for event in events:        
		for ii in range(3):
			middleX = 0.5*(event[ii].posX+event[(ii+1)%3].posX)
			middleY = 0.5*(event[ii].posY+event[(ii+1)%3].posY)
			A = -1.0*event[(ii+1)%3].posY + event[ii].posY
			B = event[(ii+1)%3].posX - event[ii].posX
			C = -1.0*event[ii].posY*B + event[ii].posX*(-1.0 * A)
			d = math.fabs(C)/math.sqrt(A*A+B*B)
			theta = math.atan2(middleY, middleX)

			if ii ==0:
				is_true_annihilation = event[ii].coincType == CoincType.kTrue
				lors.append(LOR(d, theta, is_true_annihilation)) # True for back-to-back emission
				lors_annihilation.append(LOR(d, theta, is_true_annihilation))
			else:
				lors.append(LOR(d, theta, False))
				lors_with_prompt.append(LOR(d, theta, False))
	print('[LORS FOUND]')
	return lors, lors_annihilation, lors_with_prompt


def count_sorted_lors(lors):
	d_min = 0
	d_mid = 0
	d_max = 0
	for ii in range(len(lors)/3):
		index = -10
		for jj in range(3):
			if lors[3*ii+jj].is_from_annihilation:
				index = jj
				break
		if lors[3*ii+index].d < lors[3*ii+(index+1)%3].d and lors[3*ii+index].d < lors[3*ii+(index+2)%3].d:
			d_min += 1
		elif lors[3*ii+index].d < lors[3*ii+(index+1)%3].d or lors[3*ii+index].d < lors[3*ii+(index+2)%3].d:
			d_mid += 1
		else:
			d_max += 1
	annihilation_lors_no = float(d_min+d_mid+d_max)
	# divide by the number of aniihilation lors to get fractions
	return d_min/annihilation_lors_no, d_mid/annihilation_lors_no, d_max/annihilation_lors_no


def write_goja_output(coincidences, filename='goja_output.txt'):
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