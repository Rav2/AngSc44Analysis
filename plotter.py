"""
@author: Rafal Maselek
This File contains functions for making various plots.
"""

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt


def plot_d_distribution(lors, filename):
    """
    Plots distributions of LOR parameters (d, phi)
    :param lors: List of Lor objects.
    :param filename: Name of the file with the plot.
    :return: nothing
    """
    plt.clf()
    d = []
    for lor in lors:
        d.append(lor.d)
    fig, ax = plt.subplots()
    nbins = 50 # no of bin
    if nbins > len(lors):
        nbins = len(lors)
    ax.hist(d, nbins)
    ax.set_title('distribution of distance from origin to line')
    ax.set_xlabel('d [mm]')
    ax.set_ylabel('entries [1]')
    plt.savefig("results/"+filename)
    print('[D DISTRIBUTION PLOTTED]')


def plot_classification_plots(TPR, PPV, FPR, d_tresholds, filename='classification_plots.png'):
    """
    Plots classification plots, one contains TPR and PPV as a function of d_treshold, the second one contains ROC curve.
    :param TPR: True positive rate/ Sensitivity
    :param PPV: Positive predictive value/ purity
    :param FPR: False positive rate/ Fall-out
    :param d_tresholds: List of values of treshold on the distance between the origin and lor
    :param filename: Name of the file with two plots.
    :return: nothing
    """
    fig, axarr = plt.subplots(1, 2, figsize=(10, 5))
    axarr[0].plot(d_tresholds, TPR, 'r', d_tresholds, PPV, 'b')
    axarr[0].legend(['TPR', 'PPV'])
    axarr[0].set_xlabel('r [mm]')
    axarr[0].set_ylim(0.0, 1.0)
    axarr[0].set_title('TPR and PPV as a function o r')

    axarr[1].plot(FPR, TPR, 'b', FPR, FPR, '#000000')
    axarr[1].set_title('ROC curve for 44Sc')
    axarr[1].set_xlabel('FPR')
    axarr[1].set_xlim(0.0, 1.0)
    axarr[1].set_ylim(0.0, 1.0)
    axarr[1].set_ylabel('TPR')
    plt.savefig('results/'+filename)
    print('[CLASSIFICATION PLOTS PLOTTED]')


def plot_hits(events, filename='hits.png'):
    """
    Plots at most 100 points of hits position.
    :param events: List of events.
    :param filename: Name of the file with the plot.
    :return: nothing
    """
    fig = plt.figure()
    ax = Axes3D(fig)
    if len(events) > 100:
        limit = 100
    else:
        limit = len(events)
    X = []
    Y = []
    Z = []
    for ii in range(limit):
        for kk in range(2):
            X.append(events[ii][kk].posX)
            Y.append(events[ii][kk].posY)
            Z.append(events[ii][kk].posZ)
    ax.scatter(X, Y, Z)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.savefig('results/'+filename)
    plt.show()


def plot_edep_distribution(events, filename="edep_distribution.png"):
    """
    Plots distribution of the energy deposited for annihilation and prompt gammas. 
    :param events: List of events.
    :param filename: Name of the file with the plot.
    :return: nothing
    """
    edep_511 = []
    edep_prompt =[]
    for event in events:
        edep_511.append(event[0].edep)
        edep_511.append(event[1].edep)
        edep_prompt.append(event[2].edep)
    fig, axarr = plt.subplots(1, 2, figsize=(10, 5))
    nbins = 50
    if nbins > len(edep_511):
        nbins = len(edep_511)
    axarr[0].hist(edep_511, nbins)
    axarr[0].set_title('distribution of edep of 511 keV photons')
    axarr[0].set_xlabel('edep [MeV]')
    axarr[0].set_ylabel('entries [1]')
    axarr[1].hist(edep_prompt, nbins)
    axarr[1].set_title('distribution of edep of 1157 keV photons')
    axarr[1].set_xlabel('edep [MeV]')
    axarr[1].set_ylabel('entries [1]')
    plt.savefig('results/'+filename)
    print('[EDEP 511 keV: MIN={}, MAX={}'.format(min(edep_511), max(edep_511)))
    print('[EDEP 1157 keV: MIN={}, MAX={}'.format(min(edep_prompt), max(edep_prompt)))
    print('[EDEP DISTRIBUTION PLOTTED]')


def plot_position_and_time_distribution(events, filename="position_and_time_distribution.png", only_511keV=False):
    """
    Plots distributions of x,y,z,t coordinates of hits.
    :param events: List of events.
    :param filename: Name of the file with the plot.
    :param only_511keV: If true, only distributions for 511 keV photons will be plotted.
    :return: nothing
    """
    x = []
    y =[]
    z = []
    t =[]
    add_to_x = x.append
    add_to_y = y.append
    add_to_z = z.append
    add_to_t = t.append
    for event in events:
        for ii in range(len(event)):
            if only_511keV and ii==2:
                break
            hit = event[ii]
            add_to_x(hit.posX)
            add_to_y(hit.posY)
            add_to_z(hit.posZ)
            add_to_t(hit.time)

    fig, axarr = plt.subplots(2, 2, figsize=(10, 10))
    nbins = 50
    if nbins > len(events):
        nbins = len(events)
    axarr[0,0].hist(x, nbins)
    axarr[0,0].set_title('distribution of x coordinate of hits')
    axarr[0,0].set_xlabel('x [mm]')
    axarr[0,0].set_ylabel('entries [1]')
    axarr[0,1].hist(y, nbins)
    axarr[0,1].set_title('distribution of y coordinate of hits')
    axarr[0,1].set_xlabel('y [mm]')
    axarr[0,1].set_ylabel('entries [1]')
    axarr[1,0].hist(z, nbins)
    axarr[1,0].set_title('distribution of z coordinate of hits')
    axarr[1,0].set_xlabel('z [mm]')
    axarr[1,0].set_ylabel('entries [1]')
    axarr[1,1].hist(t, nbins)
    axarr[1,1].set_title('distribution of interaction time of hits')
    axarr[1,1].set_xlabel('t [?]')
    axarr[1,1].set_ylabel('entries [1]')

    plt.savefig('results/'+filename)
    print('[POSITION&TIME DISTRIBUTION PLOTTED]')


def plot_lors_fractions(d_min, d_mid, d_max, source_pars, labels = ["d_min", "d_mid", "d_max"], style="b-r-g-", filename="lors_fractions.png"):
    """
    Plots the fractions of lors classified by the distance from origin. Every true lor belong to an event containing three
    lors, so the true lor can be either the closest one to the origin, the one in the middle, or the farthest one. 
    :param d_min: Array containing the number of lors that were closest to the origin.
    :param d_mid: Array containing the number of lors that were not the closest nor the farthest.
    :param d_max: Array containing the number of lors that were farthest from the origin.
    :param labels: Array containing labels for three sets of data.
    :param style: 6 char string defining the line style for 3 sets of data, e.g. style="b-r-g-".
    :param source_pars: Array containing values of source's radii. 
    :param filename: Name of the output file.
    :return: Axis object with the plot.
    """
    if len(d_min) != len(d_mid) != len(d_max):
        raise Exception("Improper dimensions od d arrays! Counting was done wrong!")

    plt.title('Fractions of annihilation lors as a function of source dimension')
    plt.xlabel('radius of the source cyllinder [mm]')
    plt.ylabel('fraction [1]')
    plt.plot(source_pars, d_min, style[0:2], label=labels[0])
    plt.plot(source_pars, d_mid, style[2:4], label=labels[1])
    plt.plot(source_pars, d_max, style[4:6], label=labels[2])
    plt.legend()
    plt.savefig("results/"+filename)
    print('[FRACTIONS OF LORS PLOTTED]')
    return plt.gca()


def plot_lors_fractions2(d_min1, d_mid1, d_max1, source_pars1, d_min2, d_mid2, d_max2, source_pars2, filename="lors_fractions2.png"):
    """
    Draws fractions of LORs like plot_lors_fractions but for two sets of data. Making additional function for that was
    necessary due to problems with matplotlib at the NCBJ cluster.
    :param d_min1: Array containing the number of lors that were closest to the origin for data set 1.
    :param d_mid1: Array containing the number of lors that were not the closest nor the farthest for data set 1.
    :param d_max1: Array containing the number of lors that were farthest from the origin for data set 1.
    :param source_pars1: Array containing values of source's radii for data set 1.
    :param d_min2: Array containing the number of lors that were closest to the origin for data set 2.
    :param d_mid2: Array containing the number of lors that were not the closest nor the farthest for data set 2.
    :param d_max2: Array containing the number of lors that were farthest from the origin for data set 2.
    :param source_pars2: Array containing values of source's radii for data set 2.
    :param filename: Name of the output file.
    :return: Axis object with the plot.
    """
    if len(d_min1) != len(d_mid1) != len(d_max1) or len(d_min2) != len(d_mid2) != len(d_max2):
        raise Exception("Improper dimensions od d arrays! Counting was done wrong!")

    plt.title('Fractions of annihilation lors as a function of source dimension')
    plt.xlabel('radius of the source cyllinder [mm]')
    plt.ylabel('fraction [1]')
    labels = ["d_min", "d_mid", "d_max", "d_min_phantom", "d_mid_phantom", "d_max_phantom"]
    plt.plot(source_pars1, d_min1, "b-", label=labels[0])
    plt.plot(source_pars1, d_mid1, "g-", label=labels[1])
    plt.plot(source_pars1, d_max1, "r-", label=labels[2])
    plt.plot(source_pars2, d_min2, "y-", label=labels[3])
    plt.plot(source_pars2, d_mid2, "m-", label=labels[4])
    plt.plot(source_pars2, d_max2, "c-", label=labels[5])
    plt.legend()
    plt.savefig("results/"+filename)
    print('[FRACTIONS OF LORS PLOTTED]')
    return plt.gca()

