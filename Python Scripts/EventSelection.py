#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 10:08:52 2021

@author: sb16165
"""
import sys
import os
import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# custom imports
import Master
from Master import Unwrap, Conditional, Vector3

import SelectionQuantities
from Plots import PlotHist, PlotHist2D, PlotHistComparison

# Dictionary of parameters we can caculate and cut on
paramString = {
    0  : "pandora tag" ,
    1  : "CNN score",
    2  : "nHits" ,
    3  : "energy residual",
    4  : "beam angle",
    5  : "mc angle",
    6  : "start hits",
    7  : "energy",
    8  : "mc_energy",
    9  : "pair seperation",
    10 : "pair angle",
    11 : "pair leading",
    12 : "pair second",
    13 : "invariant mass",
    14 : "true pair angle"
}

paramAxesLabels = {
    0  : "Pandora tag of daughters",
    1  : "CNN scores of beam daughters",
    2  : "Number of collection plane hits",
    3  : "Reconstruted energy residual",
    4  : "Angle between beam track and daughter shower (rad)",
    5  : "Angle between shower and MC parent (rad)",
    6  : "Hits close to the shower start",
    7  : "reconstructed shower energy (GeV)",
    8  : "mc energy of beam daughters (GeV)",
    9  : "Shower pair Separation (cm)",
    10 : "Angle between shower pairs (rad)",
    11 : "Leading shower energy (GeV)",
    12 : "Secondary shower energy (GeV)",
    13 : "Shower pair invariant mass (GeV)",
    14 : "True shower pair angle (rad)"
}


def save(name):
    """
    Saves the last created plot to file.
    ----- Parameters -----
    out                 : file path to save to
    reference_filename  : global variable, is the file name prefix of the plots
    ----------------------
    """
    out = subDirectory + name + reference_filename + ".png"
    plt.savefig(out)
    plt.close()


def InvertConditional(conditional):
    """
    Will invert a provided conditional. Uses Master.Conditional
    """
    if conditional == Conditional.GREATER:
        return Conditional.LESS
    if conditional == Conditional.LESS:
        return Conditional.GREATER
    if conditional == Conditional.EQUAL:
        return Conditional.NOT_EQUAL
    if conditional == Conditional.NOT_EQUAL:
        return Conditional.EQUAL


def GetCondString(conditional):
    """
    Will convert the meaning of a conditional to string format. Uses Master.Conditional
    """
    if conditional == Conditional.NOT_EQUAL:
        return " != "
    if conditional == Conditional.EQUAL:
        return " = "
    if conditional == Conditional.GREATER:
        return " > "
    if conditional == Conditional.LESS:
        return " < "
    else:
        return " ??? "


def StringToCond(string):
    """
    Parses the command line input for --conditional into a variable of type Master.Conditional.
    """
    if string == "ne":
        return Conditional.NOT_EQUAL
    elif string == "e":
        return Conditional.EQUAL
    elif string == "g":
        return Conditional.GREATER
    elif string == "l":
        return Conditional.LESS
    else:
        return -999


def PlotPandoraTagging(tags, label="", title=""):
    """
    Plots Pandora tagging of each daughter events in a bar plot.
    (make sure you unwrap the data first!)
    ----- Parameters -----
    unique      : type of tags in the data, can be -999, 11 or 13
    amount      : number of times a tag occurs in data
    ----------------------
    """
    unique, amount = list(np.unique(tags, return_counts=True)) # get tags and the number of occurances
    
    # switch from numerical tags to words for easy interpretation
    for i in range(len(unique)):
        if unique[i] == -999:
            unique[i] = str(unique[i]) + "\nN/A"
            continue
        if unique[i] == 11:
            unique[i] = str(unique[i]) + "\nshower"
            continue
        if unique[i] == 13:
            unique[i] = str(unique[i]) + "\ntrack"
            continue
    # plot
    plt.bar(unique, amount, label=label)
    plt.xlabel(paramAxesLabels[0])
    plt.ylabel("Number of daughters")
    plt.title(title)
    if label != "": plt.legend()
    plt.tight_layout()


def CalculateParameters(conditional=None):
    """
    will retrieve data in the root file created by the pi0 analyser.
    Then the data is used to calcuate the selection quantities and if
    a conditional is given, will apply a cut on the whole dataset.
    ----- Parameters -----
    data        : object of type Master.Data(), contains functions 
                  capable of retrieving data from the root file
    ----------------------
    """
    print("getting daughter info...")
    start_pos = data.start_pos()
    direction = data.direction()
    energy = data.energy()
    mc_energy = data.mc_energy()
    nHits = data.nHits()
    cnn_em = data.cnn_em()
    cnn_track = data.cnn_track()
    pandoraTag = data.pandoraTag()
    
    if beamData is True:
        print("getting beam info...")
        beam_start_pos = data.beam_start_pos()
        beam_end_pos = data.beam_end_pos()
    
    print("getting daughter truth info...")
    true_start_pos = data.true_start_pos()
    true_end_pos = data.true_end_pos()
    true_momentum = data.true_momentum()
    
    # custom data
    hit_radial = data.hit_radial()
    hit_longitudinal = data.hit_longitudinal()
    
    
    # calculate quantities from selected data
    print("calculate CNN score")
    cnn_score = SelectionQuantities.CNNScore(cnn_em, cnn_track)
    
    print("energy residual")
    energyResidual = (energy - mc_energy)/ mc_energy
    
    if beamData is True:
        print("angle between beam and daughters")
        beam_angle = SelectionQuantities.BeamTrackShowerAngle(beam_start_pos, beam_end_pos, direction)
    else:
        beam_angle = []
        for i in range(len(nHits)):
            evt = []
            for j in range(len(nHits[i])):
                evt.append(-999)
            beam_angle.append(evt)
        beam_angle = np.array(beam_angle, object)
    
    print("angle between daughters and mc particle")
    mc_angle = SelectionQuantities.DaughterRecoMCAngle(true_start_pos, true_end_pos, direction)
    
    print("start hits")
    start_hits = SelectionQuantities.GetShowerStartHits(hit_radial, hit_longitudinal)
    
    ### apply cuts on data that is unpaired ###
    parameters = [pandoraTag, cnn_score, nHits, energyResidual, beam_angle, mc_angle, start_hits, energy, mc_energy]
    
    if cutData is True:
        mask = Master.SelectionMask()
        mask.InitiliseMask(nHits)
        mask.CutMask(parameters[param], cut, conditional)
        
        if param < len(parameters):
            parameters = [mask.Apply(p) for p in parameters]
        else:
            print("cuts on paired values not implemented")
    
        ### cut data needed to calculate paired values ###
        start_pos = mask.Apply(start_pos)
        direction = mask.Apply(direction)
        energy = mask.Apply(energy)
    
    print("shower pairs")
    shower_pairs = SelectionQuantities.GetShowerPairs(start_pos)
    
    print("separation")
    pair_separation = SelectionQuantities.ShowerPairSeparation(start_pos, shower_pairs)
    pair_separation = np.reshape(pair_separation, (len(nHits)))
    
    print("pair angle")
    pair_angle = SelectionQuantities.ShowerPairAngle(shower_pairs, direction)
    
    print("true pair angle")
    true_dir = true_momentum.Normalise()
    true_pair_angle = SelectionQuantities.ShowerPairAngle(shower_pairs, true_dir)
    
    print("pair energy")
    pair_energies, pair_leading, pair_second = SelectionQuantities.ShowerPairEnergy(shower_pairs, energy)
    
    print("Invariant mass")
    inv_mass = SelectionQuantities.InvariantMass(pair_angle, pair_energies)
    
    parameters_pair = [pair_separation, pair_angle, pair_leading, pair_second, inv_mass, true_pair_angle]

    parameters = [Unwrap(p) for p in parameters]
    parameters_pair = [Unwrap(p) for p in parameters_pair]

    return parameters, parameters_pair


def PlotSingle(q, nbins):
    """
    Function holding code which plots all the selection quantities and saves them to the output
    directory. Constantly modified and in principle hard to condense into smaller code blocks.
    """
    print("Pandora Tag")
    PlotPandoraTagging(q[0])
    save("pandora_tag")


    print("CNN score")
    PlotHist(q[1][q[1] != -999], nbins, paramAxesLabels[1])
    save("cnn_score")


    print("collection plane hits")
    nHits_plot = q[2][q[2] != 0]
    PlotHist(nHits_plot, nbins, paramAxesLabels[2])
    save("collection_hits")


    print("hits vs mc angle")
    PlotHist2D(q[2], q[5], nbins, [-999, max(q[2])], [0, np.pi], paramAxesLabels[2], paramAxesLabels[5])
    save("hits_vs_mcAngle")


    if beamData is True:
        print("hits vs beam angle")
        PlotHist2D(q[2], q[4], nbins, [-999, max(q[2])], [0, np.pi], paramAxesLabels[2], paramAxesLabels[4])
        save("hits_vs_beamAngle")


    print("hits vs energy residual")
    PlotHist2D(q[2], q[3], nbins, [-999, max(q[2])], [-1, 1], paramAxesLabels[2], paramAxesLabels[3])
    save("hits_vs_energyRes")


    print("mc_energy vs energy residual")
    PlotHist2D(q[8], q[3], 100, [min(q[8]), max_energy], [-1, 1], paramAxesLabels[8], paramAxesLabels[3])
    save("mcEnergy_vs_energyRes")


    print("reco energy vs energy residual")
    PlotHist2D(q[7], q[3], 100, [min(q[7]), max_energy], [-1, 1], paramAxesLabels[7], paramAxesLabels[3])
    save("recoEnergy_vs_energyRes")


    if beamData is True:
        print("beam angle")
        q[4] = q[4][q[4] != -999]
        PlotHist(q[4], nbins, paramAxesLabels[4])
        save("beam_daughterAngle")


    print("pair separation")
    q[9] = q[9][q[9] > 0]
    PlotHist(q[9][q[9] < 51], nbins, paramAxesLabels[9])
    save("pair_separation")


    print("pair angle")
    PlotHist(q[10][q[10] != -999], nbins, paramAxesLabels[10])
    save("shower_pairAngle")


    print("pair energies")
    PlotHistComparison(q[11][q[11] > 0], q[12][q[12] > 0], nbins, xlabel="Shower energy (GeV)", label_1=paramAxesLabels[11], label_2=paramAxesLabels[12], alpha=0.5, density=False)
    save("daughter_energies")


    print("invariant mass")
    inv_mass = q[13][q[13] > 0]
    PlotHist(inv_mass[inv_mass < 0.5], nbins, paramAxesLabels[13])
    save("inv_mass")


    print("true energy vs reco shower energy")
    PlotHist2D(q[8][q[7] > 0], q[7][q[7] > 0], nbins, [min(q[8]), max_energy], [min(q[7]), max_energy], xlabel=paramAxesLabels[8], ylabel=paramAxesLabels[7])
    save("mc_reco_energy")

    
    print("true pair angle")
    PlotHist(q[14], nbins, paramAxesLabels[14])
    save("true_pair_angle")

    
    print("true pair angle vs pair angle")
    PlotHist2D(q[14], q[10], nbins, xlabel=paramAxesLabels[14], ylabel=paramAxesLabels[10])
    save("true_reco_angle")

    
    print("opening angle residual")
    angle_residual = q[10] - q[14]
    PlotHist(angle_residual, nbins, r"$\theta_{reco} - \theta_{true}$ (rad)")
    save("angle_residual")


    print("start hits")
    print("Clyinder with radius " + str(r) + "cm, length " + str(l_max - l_min) + "cm, shower start is offset from the center by " + str(l_max - (l_max - l_min)/2) + "cm")
    PlotHist(q[6], nbins, paramAxesLabels[6])
    save("start_hits")
    
    print("reco shower energy")
    PlotHist(q[7][q[7] > 0], nbins, paramAxesLabels[7])
    save("reco_energy")


    print("true energy")
    true_energy = q[8][q[8] > 0]
    true_energy = true_energy[true_energy <= max_energy]
    PlotHist(true_energy, nbins, paramAxesLabels[8])
    save("true_energy")


    print("nHits vs true energy")
    PlotHist2D(q[2], q[8], nbins, [-999, max(q[2])], [0, max_energy], xlabel=paramAxesLabels[2], ylabel=paramAxesLabels[8])
    save("nHits_vs_trueEnergy")


    print("nHits vs reco shower energy")
    PlotHist2D(q[2], q[7], nbins, [-999, max(q[2])], [0, max(q[7])], xlabel=paramAxesLabels[2], ylabel=paramAxesLabels[7])
    save("nHits_vs_recoEnergy")
    
    print("opening angles vs start hits?")
    #PlotHist2D(q[10], q[6])


def PlotComparison(param, conditional, conditional_inv, cut, q_1, q_2, nbins):
    """
    Function holding code which plots all the selection quantities from both the selected and
    rejected data, saving them to the output directory.
    Constantly modified and in principle hard to condense into smaller code blocks.
        ----- Parameters -----
    l_1        : label for selected data, describes the cut made
    l_2        : label for rejected data, describes the cut made
    ----------------------
    """
    l_1 = paramString[param] + GetCondString(conditional) + str(cut)
    l_2 = paramString[param] + GetCondString(conditional_inv) + str(cut)
    
    print("Pandora Tag")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotPandoraTagging(q_1[0], title=l_1)
    plt.subplot(122)
    PlotPandoraTagging(q_2[0], title=l_2)
    save("pandora_tag")
    
    
    print("CNN score")
    PlotHistComparison(q_1[1][q_1[1] != -999], q_2[1][q_2[1] != -999], nbins, alpha=0.5, xlabel=paramAxesLabels[1], label_1=l_1, label_2=l_2)
    save("cnn_score")
    
    
    print("collection plane hits")
    nHits_1 = q_1[2][q_1[2] != 0]
    nHits_2 = q_2[2][q_2[2] != 0]
    
    PlotHistComparison(nHits_1[nHits_1 < 101], nHits_2[nHits_2 < 101], nbins, alpha=0.5, xlabel=paramAxesLabels[2], label_1=l_1, label_2=l_2)
    save("collection_hits")
    
    
    print("hits vs mc angle")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[2], q_1[5], nbins, [-999, 510], [0, np.pi], paramAxesLabels[2], paramAxesLabels[5], title=l_1)
    
    plt.subplot(122)
    PlotHist2D(q_2[2], q_2[5], nbins, [-999, 510], [0, np.pi], paramAxesLabels[2], paramAxesLabels[5], title=l_2)
    save("hits_vs_mcAngle")
    
    
    if beamData is True:
        print("hits vs beam angle")
        plt.figure(2, (12, 5))
        plt.subplot(121)
        PlotHist2D(q_1[2], q_1[4], nbins, [-999, 501], [0, np.pi], paramAxesLabels[2], paramAxesLabels[4], title=l_1)
        plt.subplot(122)
        PlotHist2D(q_2[2], q_2[4], nbins, [-999, 501], [0, np.pi], paramAxesLabels[2], paramAxesLabels[4], title=l_2)
        save("hits_vs_beamAngle")
    
    
    print("hits vs energy residual")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[2], q_1[3], nbins, [-999, max(q_1[2])], [-1, 1], paramAxesLabels[2], paramAxesLabels[3], title=l_1)
    plt.subplot(122)
    PlotHist2D(q_2[2], q_2[3], nbins, [-999, max(q_2[2])], [-1, 1], paramAxesLabels[2], paramAxesLabels[3], title=l_2)
    save("hits_vs_energyRes")
    
    
    print("mc energy vs energy residual")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[8], q_1[3], nbins, [min(q_1[8]), max_energy], [-1, 1], paramAxesLabels[8], paramAxesLabels[3], title=l_1)
    plt.subplot(122)
    PlotHist2D(q_2[8], q_2[3], nbins, [min(q_2[8]), max_energy], [-1, 1], paramAxesLabels[8], paramAxesLabels[3], title=l_2)
    save("mc_energy_vs_energyRes")
    
    
    print("energy vs energy residual")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[7], q_1[3], nbins, [min(q_1[7]), max_energy], [-1, 1], paramAxesLabels[7], paramAxesLabels[3], title=l_1)
    plt.subplot(122)
    PlotHist2D(q_2[7], q_2[3], nbins, [min(q_2[7]), max_energy], [-1, 1], paramAxesLabels[7], paramAxesLabels[3], title=l_2)
    save("energy_vs_energyRes")
    
    
    if beamData is True:
        print("beam angle")
        beam_angle_1 = q_1[4][q_1[4] != -999]
        beam_angle_2 = q_2[4][q_2[4] != -999]
        PlotHistComparison(beam_angle_1, beam_angle_2, nbins, alpha=0.5, xlabel=paramAxesLabels[4], label_1=l_1, label_2=l_2)
        save("beam_daughterAngle")
    
    
    print("pair separation")
    pair_separation_1 = q_1[9][q_1[9] > 0]
    pair_separation_2 = q_2[9][q_2[9] > 0]
    PlotHistComparison(pair_separation_1[pair_separation_1 < 51], pair_separation_2[pair_separation_2 < 51], nbins, alpha=0.5, xlabel=paramAxesLabels[9], label_1=l_1, label_2=l_2)
    save("pair_separation")
    
    
    print("pair angle")
    pair_angle_1 = q_1[10][q_1[10] != -999]
    pair_angle_2 = q_2[10][q_2[10] != -999]
    PlotHistComparison(pair_angle_1, pair_angle_2, nbins, alpha=0.5, xlabel=paramAxesLabels[10], label_1=l_1, label_2=l_2)
    save("shower_pairAngle")
    
    
    print("leading shower energy")
    q_1[11] = q_1[11].flatten()
    q_2[11] = q_2[11].flatten()
    PlotHistComparison(q_1[11][q_1[11] > 0], q_2[11][q_2[11] > 0], nbins, alpha=0.5, xlabel=paramAxesLabels[11], label_1=l_1, label_2=l_2)
    save("leading_shower_energies")
    
    
    print("sub-leading shower energy")
    q_1[12] = q_1[12].flatten()
    q_2[12] = q_2[12].flatten()
    PlotHistComparison(q_1[12][q_1[12] > 0], q_2[12][q_2[12] > 0], nbins, alpha=0.5, xlabel=paramAxesLabels[12], label_1=l_1, label_2=l_2)
    save("subleading_shower_energies")
    
    
    print("invariant mass")
    inv_mass_1 = q_1[13][q_1[13] > 0]
    inv_mass_2 = q_2[13][q_2[13] > 0]
    PlotHistComparison(inv_mass_1[inv_mass_1 < 0.5], inv_mass_2[inv_mass_2 < 0.5], nbins, alpha=0.5, xlabel=paramAxesLabels[13], sf=3, label_1=l_1, label_2=l_2)
    save("inv_mass")


    print("true energy vs reco shower energy")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[8][q_1[7] > 0], q_1[7][q_1[7] > 0], nbins, [min(q_1[8]), max_energy], [min(q_1[7]), max_energy], xlabel=paramAxesLabels[8], ylabel=paramAxesLabels[7], title=l_1)
    plt.subplot(122)
    PlotHist2D(q_2[8][q_2[7] > 0], q_2[7][q_2[7] > 0], nbins, [min(q_2[8]), max_energy], [min(q_2[7]), max_energy], xlabel=paramAxesLabels[8], ylabel=paramAxesLabels[7], title=l_2)
    save("mc_reco_energy")
    
    
    print("true pair angle")
    PlotHistComparison(q_1[14], q_2[14], nbins, alpha=0.5, xlabel=paramAxesLabels[14], label_1=l_1, label_2=l_2)
    save("true_pair_angle")

    
    print("true pair angle vs pair angle")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[14], q_1[10], nbins, xlabel=paramAxesLabels[14], ylabel=paramAxesLabels[10], title=l_1)
    plt.subplot(122)
    PlotHist2D(q_2[14], q_2[10], nbins, xlabel=paramAxesLabels[14], ylabel=paramAxesLabels[10], title=l_2)
    save("true_reco_angle")

    
    print("opening angle residual")
    angle_residual_1 = q_1[10] - q_1[14]
    angle_residual_2 = q_2[10] - q_2[14]
    PlotHistComparison(angle_residual_1, angle_residual_2, nbins, alpha=0.5, xlabel=r"$\theta_{reco} - \theta_{true}$ (rad)", label_1=l_1, label_2=l_2)
    save("angle_residual")


    print("start hits")
    print("Clyinder with radius " + str(r) + "cm, length " + str(l_max - l_min) + "cm, shower start is offset from the center by " + str(l_max - (l_max - l_min)/2) + "cm")
    PlotHistComparison(q_1[6], q_2[6], nbins, alpha=0.5, xlabel=paramAxesLabels[6], label_1=l_1, label_2=l_2)
    save("start_hits")


    print("reco shower energy")
    PlotHistComparison(q_1[7][q_1[7] > 0], q_2[7][q_2[7] > 0], nbins, alpha=0.5, xlabel=paramAxesLabels[7], label_1=l_1, label_2=l_2)
    save("reco_energy")


    print("true energy")
    true_energy_1 = q_1[8][q_1[8] > 0]
    true_energy_1 = true_energy_1[true_energy_1 <= max_energy]
    true_energy_2 = q_1[8][q_1[8] > 0]
    true_energy_2 = true_energy_2[true_energy_2 <= max_energy]
    PlotHistComparison(true_energy_1, true_energy_2, nbins, alpha=0.5, xlabel=paramAxesLabels[8], label_1=l_1, label_2=l_2)
    save("true_energy")


    print("nHits vs true energy")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[2], q_1[8], nbins, [-999, max(q_1[2])], [0, max_energy], xlabel=paramAxesLabels[2], ylabel=paramAxesLabels[8], title=l_1)
    plt.subplot(122)
    PlotHist2D(q_2[2], q_2[8], nbins, [-999, max(q_1[2])], [0, max_energy], xlabel=paramAxesLabels[2], ylabel=paramAxesLabels[8], title=l_2)
    save("nHits_vs_trueEnergy")
    

    print("nHits vs reco shower energy")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[2], q_1[7], nbins, [-999, max(q_1[2])], [0, max(q_1[7])], xlabel=paramAxesLabels[2], ylabel=paramAxesLabels[7], title=l_1)
    plt.subplot(122)
    PlotHist2D(q_2[2], q_2[7], nbins, [-999, max(q_2[2])], [0, max(q_2[7])], xlabel=paramAxesLabels[2], ylabel=paramAxesLabels[8], title=l_2)
    save("nHits_vs_recoEnergy")


start_time = time.time()

# handle command line imports
args = sys.argv

if "--help" in args or len(args) == 1:
    print("usage: --file <root file> --outName <file name to append> --outDir <output directory> --parameter <0 - 8> --conditional < l, g, ne, e > --cut <cut value> --beam <1/0> --plot <single/both>")
    exit(0)
if "--file" in args:
    root_filename = args[args.index("--file") + 1]
else:
    print("no file chosen!")
    exit(1)
if "--outName" in args:
    reference_filename = "_" + args[args.index("--outName") + 1]
else:
    reference_filename = "_" + root_filename[19:-5]
if "--outDir" in args:
    subDirectory = args[args.index("--outDir") + 1] + "/"
else:
    subDirectory = root_filename[19:-5] + "/"
if "--beam" in args:
    beamData = bool(int(args[args.index("--beam") + 1]))
else:
    beamData = True
if "--parameter" in args:
    param = int(args[args.index("--parameter") + 1])
    if "--conditional" in args:
        conditional = StringToCond(args[args.index("--conditional") + 1])
        if conditional == -999:
            print("Invalid input check help text:")
            exit(1)
    else:
        print("need to specify a condition!")
    if "--cut" in args:
        cut = float(args[args.index("--cut") + 1])
    else:
        print("need to specify a cut!")
        exit(1)
    cutData = True
    if "--plot" in args:
        plotType = args[args.index("--plot") + 1]
    else:
        plotType = "both"
else:
    cutData = False
    conditional = None


# start hit constants
r = 1 # cm
l_min = -1 # cm
l_max = 4 # cm

max_energy = 0.5 # beam momentum/energy

print("Opening root file: " + root_filename)
os.makedirs(subDirectory, exist_ok=True)
data = Master.Data(root_filename)


q_1 = Unwrap(CalculateParameters(conditional))

if cutData is True:
    conditional_inv = InvertConditional(conditional)
    q_2 = Unwrap(CalculateParameters(conditional_inv))


print("Saving Plots:")
nbins = 100

if cutData is True:
    if plotType == "single":
        print("Selection: " + paramString[param] + GetCondString(conditional) + str(cut))
        PlotSingle(q_1, nbins)
    elif plotType == "both":
        PlotComparison(param, conditional, conditional_inv, cut, q_1, q_2, nbins)

else:
    PlotSingle(q_1, nbins)

print("done!")
print("ran in %s seconds. " % (time.time() - start_time) )

