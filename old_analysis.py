from ROOT import gROOT, TCanvas, TH1, TH2, TTree, TFile
from matplotlib import pyplot as plt
import math
import numpy as np
gROOT.Reset()


def manage_bufor(array):
    new_buf = []
    #print(len(array))
    for buf in array:
        if buf[-2] == 1 and buf[4] >= 0.01: # check nCrystalCompton number and edep
            new_buf.append(buf)
    # check if there are 2 hits  from different strips
    if len(new_buf) == 2 and new_buf[0][5] != new_buf[1][5]:
        dist = 0.5*math.sqrt((new_buf[0][0]-new_buf[1][0])**2+(new_buf[0][1]-new_buf[1][1])**2)
        #if dist < 300:
            #print('dist=', dist)
            #print('X, Y, Z, T, Edep, volumeID, nCrystalCmpton, eventID')
            #print(new_buf[0])
            #print(new_buf[1])
        return new_buf
    else:
        #print('invalid bufor length=', len(new_buf))
        return []
    
    
# reading data
coords_511 = []
coords_prompt = []

f_511 = TFile("511keV_1.root")
bufor = []
for event in f_511.Hits:
    if len(bufor) == 0 or bufor[0][-1] == event.eventID:
        bufor.append((event.posX, event.posY, event.posZ, event.time, event.edep, event.volumeID[1], event.nCrystalCompton, event.eventID))
    else:
        for arr in manage_bufor(bufor):
            coords_511.append(arr)
        bufor = [(event.posX, event.posY, event.posZ, event.time, event.edep, event.volumeID[1], event.nCrystalCompton, event.eventID)]
        
      
f_prompt = TFile("prompt_1.root")
for event in f_prompt.Hits :
    if event.nCrystalCompton == 1 and event.edep >= 0.01:
      coords_prompt.append((event.posX, event.posY, event.posZ, event.time, event.edep, event.volumeID[1], event.nCrystalCompton, event.eventID))

# truncating the data to equal number of events      
n_511 = len(coords_511)
n_prompt = len(coords_prompt)

if n_511//2 > n_prompt:
    n_events = n_prompt
else:
    n_events = n_511//2

print('NUMBER OF EVENTS=', n_events)
    
x_511 = [coord[0] for coord in coords_511[:2*n_events]]
y_511 = [coord[1] for coord in coords_511[:2*n_events]]

x_prompt = [coord[0] for coord in coords_prompt[:n_events]]
y_prompt = [coord[1] for coord in coords_prompt[:n_events]]

# Validation of data
fig, ax = plt.subplots()
#ax.scatter(x_511[:2], y_511[2:4])
print('max x_511=',max(x_511),'max y_511=',max(y_511))
print('min x_511=',min(x_511),'min y_511=',min(y_511))
#plt.show()

#ax.scatter(x_prompt, y_prompt)
print('max x_prompt=',max(x_prompt),'max y_prompt=',max(y_prompt))
print('min x_prompt=',min(x_prompt),'min y_prompt=',min(y_prompt))
#plt.show()

# putting coordinates together into events
events = []
for ii in range(n_events):
    events.append([(x_511[2*ii], y_511[2*ii]), (x_511[2*ii+1], y_511[2*ii+1]), (x_prompt[ii], y_prompt[ii])])
    #dist = 0.5*math.sqrt((x_511[2*ii+1]-x_511[2*ii])**2 + (y_511[2*ii+1]-y_511[2*ii])**2)
    #print('dist=', dist)

# calculating lors' parameters in (d, theta) representation
lors = []
lors_true_d = []
for event in events:        
    for ii in range(3):
        theta = math.atan((event[(ii+1)%3][0] - event[ii][0])/(-event[(ii+1)%3][1] + event[ii][1]))
        d = event[ii][0]*math.cos(theta) + event[ii][1]*math.sin(theta)
        if d<0:
            d = -1.0*d
            if theta >0:
                theta -= math.pi
            else:
                theta += math.pi
        if ii ==0:
            lors.append((d, theta, True)) # True for back-to-back emission
            lors_true_d.append(d)
        else:
            lors.append((d, theta, False))

# plotting the histogram of d    
plt.clf()
#plt.hist([x[1] for x in lors], 220)
#plt.show()
plt.hist([x[0] for x in lors], 220)
plt.savefig('d_distr.png')
plt.clf()

plt.hist(lors_true_d, 220)
plt.savefig('d_distr_true.png')
    
#calculation of efficiency, purity etc.
lors = np.asarray(lors)
d_tresholds = np.linspace(0, 437.3, 400)
TPR = [] #sensitivity/efficiency
SPC = [] #specifity
FPR = []
PPV = [] # purity
for d_t in d_tresholds:
    FP=FN=TP=TN=0
    cond = np.count_nonzero([lors[:, 0]>d_t, lors[:, 2]], axis=0) # each entry corresponds to a LOR. If it is equal 0 then it is FN if 2 then it is FP if 1 then it is TN or TP
    FP = np.count_nonzero(cond==2)
    FN = np.count_nonzero(cond==0)
    cond2 = np.count_nonzero([lors[:, 0]>d_t, lors[:, 2]==False], axis=0)
    TP = np.count_nonzero(cond2 == 2)
    TN = np.count_nonzero(cond2 == 0)    
    TPR.append(float(TP)/float(TP+FN))
    SPC.append(TN/float(TN+FP))
    PPV.append(TP/float(TP+FP))
    FPR.append(1.0-TN/float(TN+FP))

fig, axarr = plt.subplots(1, 2, figsize=(10, 5))
axarr[0].plot(d_tresholds, TPR, d_tresholds, PPV)
axarr[0].legend(['TPR', 'PPV'])
axarr[0].set_xlabel('d treshold [mm]')
axarr[0].set_title('TPR and PPV as a function o d_treshold')

axarr[1].plot(FPR, TPR, FPR, FPR)
axarr[1].set_title('ROC curve for 44Sc')
axarr[1].set_xlabel('FPR')
axarr[1].set_ylabel('TPR')
plt.savefig('plots.png')


