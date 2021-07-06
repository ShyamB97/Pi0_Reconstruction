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
from Master import Unwrap, Vector3, Conditional

import SelectionQuantities
from Plots import Plot, PlotHist, PlotHist2D, PlotHistComparison


def save(name):
    out = subDirectory + name + reference_filename + ".png"
    plt.savefig(out)
    plt.close()


def InvertConditional(conditional):
    if conditional == Conditional.GREATER:
        return Conditional.LESS
    if conditional == Conditional.LESS:
        return Conditional.GREATER
    if conditional == Conditional.EQUAL:
        return Conditional.NOT_EQUAL
    if conditional == Conditional.NOT_EQUAL:
        return Conditional.EQUAL


def GetCondString(conditional):
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


def CalculateParameters(conditional=None):
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
    
    #reimplement once fixed!
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
    
    print("separation")
    pair_separation = SelectionQuantities.ShowerPairSeparation(start_pos)
    pair_separation = np.reshape(pair_separation, (len(nHits), ))
    
    print("pair angle")
    pair_angle = SelectionQuantities.ShowerPairAngle(start_pos, direction)
    
    print("pair energy")
    pair_energies, pair_leading, pair_second = SelectionQuantities.ShowerPairEnergy(start_pos, energy)
    
    print("Invariant mass")
    inv_mass = SelectionQuantities.InvariantMass(start_pos, direction, energy)
    
    
    parameters_pair = [pair_separation, pair_angle, pair_leading, pair_second, inv_mass]

    parameters = [Unwrap(p) for p in parameters]
    parameters_pair = [Unwrap(p) for p in parameters_pair]

    return parameters, parameters_pair


def PlotSingle(q, nbins):
    print("Pandora Tag")
    unique, amount = np.unique(q[0], return_counts=True)
    unique = list(unique)
    
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

    plt.bar(unique, amount)
    plt.xlabel("Pandora tag of daughters")
    plt.ylabel("Number of daughters")
    plt.tight_layout()
    save("pandora_tag")

    print("CNN score")
    PlotHist(q[1][q[1] != -999], nbins, "CNN scores of beam daughters")
    save("cnn_score")


    print("collection plane hits")
    nHits_plot = q[2][q[2] != 0]
    PlotHist(nHits_plot[nHits_plot < 101], nbins, "Number of collection plane hits")
    save("collection_hits")


    print("hits vs mc angle")
    PlotHist2D(q[2], q[5], nbins, [-999, max(q[2])], [0, np.pi], "Number of collection plane hits", "Angle between shower and MC parent (rad)")
    save("hits_vs_mcAngle")


    if beamData is True:
        print("hits vs beam angle")
        PlotHist2D(q[2], q[4], nbins, [-999, max(q[2])], [0, np.pi], "Number of collection plane hits", "Angle between beam track and daughter shower (rad)")
        save("hits_vs_beamAngle")


    print("hits vs energy residual")
    PlotHist2D(q[2], q[3], nbins, [-999, max(q[2])], [-1, 1], "Number of collection plane hits", "Reconstruted energy residual")
    save("hits_vs_energyRes")


    print("mc_energy vs energy residual")
    PlotHist2D(q[8], q[3], 100, [min(q[8]), max(q[8])], [-1, 1], "mc energy of beam daughters (GeV)", "Reconstruted energy residual")
    save("mcEnergy_vs_energyRes")


    print("reco energy vs energy residual")
    PlotHist2D(q[7], q[3], 100, [min(q[7]), max(q[7])], [-1, 1], "reco energy of beam daughters (GeV)", "Reconstruted energy residual")
    save("recoEnergy_vs_energyRes")


    if beamData is True:
        print("beam angle")
        q[4] = q[4][q[4] != -999]
        PlotHist(q[4], nbins, "Angle between beam track and daughter shower (rad)")
        save("beam_daughterAngle")


    print("pair separation")
    q[9] = q[9][q[9] > 0]
    PlotHist(q[9][q[9] < 51], nbins, "Shower pair Separation (cm)")
    save("pair_separation")


    print("pair angle")
    q[10] = q[10][q[10] != -999]
    PlotHist(q[10], nbins, "Angle between shower pairs (rad)")
    save("shower_pairAngle")


    print("pair energies")
    pair_leading = q[11][q[11] > 0]
    pair_second = q[12][q[12] > 0]
    PlotHistComparison(pair_leading, pair_second, nbins, xlabel="Shower energy (GeV)", label_1="shower with the most energy in a pair", label_2="shower with the least energy in a pair", alpha=0.5, density=False)
    save("daughter_energies")


    print("invariant mass")
    inv_mass = q[13][q[13] > 0]
    PlotHist(inv_mass[inv_mass < 0.5], nbins, "Shower pair invariant mass (GeV)", "")
    save("inv_mass")

    print("true energy vs reco shower energy")
    PlotHist2D(q[8][q[7] > 0], q[7][q[7] > 0], nbins, xlabel="true energy (GeV)", ylabel="reconstructed shower energy(GeV)")
    save("mc_reco_energy")

    """
    NEED TO IMPLEMENT
    print("true opening angle")
    PlotHist(true_opening_angle, nbins, "True shower pair angle (rad)")
    save("true_opening_angle")
    
    print("true pair angle vs pair angle")
    PlotHist2D(true_opening_angle, pair_angle, nbins, xlabel="True shower pair angle (rad)", ylabel="Angle between shower pairs (rad)")
    save("true_reco_angle")
    
    print("opening angle residual")
    angle_residual = pair_angle - true_opening_angle
    PlotHist(angle_residual, nbins, r"$\theta_{reco} - \theta_{true}$ (rad)")
    save("angle_residual")
    """


def PlotComparison(param, conditional, conditional_inv, cut, q_1, q_2, nbins):
    l_1 = paramString[param] + GetCondString(conditional) + str(cut)
    l_2 = paramString[param] + GetCondString(conditional_inv) + str(cut)
    print("Pandora Tag")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    
    unique, amount = list(np.unique(q_1[0], return_counts=True))
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
    plt.bar(unique, amount)
    plt.xlabel("Pandora tag of daughters")
    plt.ylabel("Number of daughters")
    plt.title(l_1)
    plt.tight_layout()
    
    plt.subplot(122)
    
    unique, amount = list(np.unique(q_2[0], return_counts=True))
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
    plt.xlabel("Pandora tag of daughters")
    plt.ylabel("Number of daughters")
    plt.title(l_2)
    plt.tight_layout()
    
    save("pandora_tag")
    
    
    print("CNN score")
    PlotHistComparison(q_1[1][q_1[1] != -999], q_2[1][q_2[1] != -999], nbins, alpha=0.5, xlabel="CNN scores of beam daughters", label_1=l_1, label_2=l_2)
    save("cnn_score")
    
    
    print("collection plane hits")
    nHits_1 = q_1[2][q_1[2] != 0]
    nHits_2 = q_2[2][q_2[2] != 0]
    
    PlotHistComparison(nHits_1[nHits_1 < 101], nHits_2[nHits_2 < 101], nbins, alpha=0.5, xlabel="Number of collection plane hits", label_1=l_1, label_2=l_2)
    save("collection_hits")
    
    
    print("hits vs mc angle")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[2], q_1[5], nbins, [-999, 510], [0, np.pi], "Number of collection plane hits", "Angle between shower and MC parent (rad)", title=l_1)
    
    plt.subplot(122)
    PlotHist2D(q_2[2], q_2[5], nbins, [-999, 510], [0, np.pi], "Number of collection plane hits", "Angle between shower and MC parent (rad)", title=l_2)
    save("hits_vs_mcAngle")
    
    
    if beamData is True:
        print("hits vs beam angle")
        plt.figure(2, (12, 5))
        plt.subplot(121)
        PlotHist2D(q_1[2], q_1[4], nbins, [-999, 501], [0, np.pi], "Number of collection plane hits", "Angle between beam track and daughter shower (rad)", title=l_1)
        plt.subplot(122)
        PlotHist2D(q_2[2], q_2[4], nbins, [-999, 501], [0, np.pi], "Number of collection plane hits", "Angle between beam track and daughter shower (rad)", title=l_2)
        save("hits_vs_beamAngle")
    
    
    print("hits vs energy residual")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[2], q_1[3], nbins, [-999, max(q_1[2])], [-1, 1], "Number of collection plane hits", "Reconstruted energy residual", title=l_1)
    plt.subplot(122)
    PlotHist2D(q_2[2], q_2[3], nbins, [-999, max(q_2[2])], [-1, 1], "Number of collection plane hits", "Reconstruted energy residual", title=l_2)
    save("hits_vs_energyRes")
    
    
    print("mc energy vs energy residual")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[8], q_1[3], nbins, [-999, max(q_1[10])], [-1, 1], "mc energy of beam daughters (GeV)", "Reconstruted energy residual", title=l_1)
    plt.subplot(122)
    PlotHist2D(q_2[8], q_2[3], nbins, [-999, max(q_2[10])], [-1, 1], "mc energy of beam daughters (GeV)", "Reconstruted energy residual", title=l_2)
    save("mc_energy_vs_energyRes")
    
    
    print("energy vs energy residual")
    plt.figure(2, (12, 5))
    plt.subplot(121)
    PlotHist2D(q_1[7], q_1[3], nbins, [-999, max(q_1[11])], [-1, 1], "reco energy of beam daughters (GeV)", "Reconstruted energy residual", title=l_1)
    plt.subplot(122)
    PlotHist2D(q_2[7], q_2[3], nbins, [-999, max(q_2[11])], [-1, 1], "reco energy of beam daughters (GeV)", "Reconstruted energy residual", title=l_2)
    save("energy_vs_energyRes")
    
    
    if beamData is True:
        print("beam angle")
        beam_angle_1 = q_1[4][q_1[4] != -999]
        beam_angle_2 = q_2[4][q_2[4] != -999]
        PlotHistComparison(beam_angle_1, beam_angle_2, nbins, alpha=0.5, xlabel="Angle between beam track and daughter shower (rad)", label_1=l_1, label_2=l_2)
        save("beam_daughterAngle")
    
    
    print("pair separation")
    pair_separation_1 = q_1[9][q_1[9] > 0]
    pair_separation_2 = q_2[9][q_2[9] > 0]
    PlotHistComparison(pair_separation_1[pair_separation_1 < 51], pair_separation_2[pair_separation_2 < 51], nbins, alpha=0.5, xlabel="Shower pair Separation (cm)", label_1=l_1, label_2=l_2)
    save("pair_separation")
    
    
    print("pair angle")
    pair_angle_1 = q_1[10][q_1[10] != -999]
    pair_angle_2 = q_2[10][q_2[10] != -999]
    PlotHistComparison(pair_angle_1, pair_angle_2, nbins, alpha=0.5, xlabel="Angle between shower pairs (rad)", label_1=l_1, label_2=l_2)
    save("shower_pairAngle")
    
    
    print("leading shower energy")
    q_1[11] = q_1[11].flatten()
    q_2[11] = q_2[11].flatten()
    PlotHistComparison(q_1[11][q_1[11] > 0], q_2[11][q_2[11] > 0], nbins, alpha=0.5, xlabel="Shower energy (GeV)", label_1=l_1, label_2=l_2)
    save("leading_shower_energies")
    
    
    print("sub-leading shower energy")
    q_1[8] = q_1[8].flatten()
    q_2[8] = q_2[8].flatten()
    PlotHistComparison(q_1[12][q_1[12] > 0], q_2[12][q_2[12] > 0], nbins, alpha=0.5, xlabel="Shower energy (GeV)", label_1=l_1, label_2=l_2)
    save("subleading_shower_energies")
    
    
    print("invariant mass")
    inv_mass_1 = q_1[13][q_1[13] > 0]
    inv_mass_2 = q_2[13][q_2[13] > 0]
    PlotHistComparison(inv_mass_1[inv_mass_1 < 0.5], inv_mass_2[inv_mass_2 < 0.5], nbins, alpha=0.5, xlabel="Shower pair invariant mass (GeV)", sf=3, label_1=l_1, label_2=l_2)
    save("inv_mass")


start_time = time.time()

# handle command line imports
args = sys.argv

if "--help" in args or len(args) == 1:
    print("usage: --file <root file> --outName <file name to append> --outDir <output directory> --parameter <0 - 8> --conditional < l, g, ne, e > --cut <cut value> --beam <1/0> --plot <single/both>")
    exit(0)
if "--file" in args:
    root_filename = args[args.index("--file") + 1]
    print(root_filename)
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


paramString = {0 :"pandora tag" , 1 : "CNN score", 2 : "nHits" , 3 : "energy residual", 4 : "beam angle", 5 : "mc angle", 6 : "start_hits", 7 : "energy", 8 : "mc_energy",
               9 : "pair seperation", 10 : "pair angle", 11 : "pair leading", 12 : "pair second", 13 : "invariant mass"}

# start hit constants
r = 1 # cm
l_min = -1 # cm
l_max = 4 # cm

print("Opening root file: " + root_filename)
os.makedirs(subDirectory, exist_ok=True)
data = Master.Data(root_filename)
#beamData = False
#param = 1 # depends on the parameters data structure
#cut = 0.6 # where to cut
#conditional = Conditional.GREATER # what to keep


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

