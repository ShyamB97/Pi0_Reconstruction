#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 10:46:49 2021

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
#from Plots import Plot, PlotHist, PlotHist2D
import Plots


def save(name):
    out = subDirectory + "/" + subDirectory + "_" + name + ".png"
    plt.savefig(out)
    plt.close()


def Calculate(file):
    data = Master.Data(file)
    cnn_em = data.cnn_em()
    cnn_track = data.cnn_track()

    # create the mask
    mask = Master.SelectionMask()
    # initilise mask (set the shape and set all values to 1)
    mask.InitiliseMask(cnn_em)
    
    ### Hacked CNN score fix ###
    mask999 = Master.SelectionMask()
    mask999.InitiliseMask(cnn_em)
    mask999.CutMask(cnn_em, -999, Conditional.NOT_EQUAL)
    
    cnn_em = mask999.Apply(cnn_em)
    cnn_track = mask999.Apply(cnn_track)
    
    cnn_score = SelectionQuantities.CNNScore(cnn_em, cnn_track)
    print("beam events: " + str(len(cnn_score)))

    #return Unwrap(cnn_score)
    return cnn_score


def mask(data_1, data_2, reference):
    # make sure data_2 < data_1
    # reference = data_1
    mask = []
    for i in range(len(data_1)):
        x = np.where(data_2 == data_1[i])[0]
        if len(x) > 0:
            mask.append(1)
        else:
            mask.append(0)

    masked = []
    for i in range(len(reference)):
        if mask[i] == 1:
            masked.append(reference[i])
    masked = np.array(masked, dtype=object)
    return masked


def ProcessData(file):
    print("Opening root file: " + file)
    data = Master.Data(file)
    
    print("getting daughter info...")
    start_pos = data.start_pos()
    direction = data.direction()
    energy = data.energy()
    mc_energy = data.mc_energy()
    nHits = data.nHits()
    cnn_em = data.cnn_em()
    cnn_track = data.cnn_track()
    
    print("getting beam info...")
    beam_start_pos = data.beam_start_pos()
    beam_end_pos = data.beam_end_pos()
    
    print("getting daughter truth info...")
    true_start_pos = data.true_start_pos()
    true_end_pos = data.true_end_pos()
    
    ### Hacked CNN score fix ###
    mask999 = Master.SelectionMask()
    mask999.InitiliseMask(cnn_em)
    mask999.CutMask(cnn_em, -999, Conditional.NOT_EQUAL)
    
    cnn_em = mask999.Apply(cnn_em)
    cnn_track = mask999.Apply(cnn_track)
    
    cnn_score = SelectionQuantities.CNNScore(cnn_em, cnn_track)
    print("beam events: " + str(len(cnn_score)))
    cnn_score = Unwrap(cnn_score)
    print("daughter events: " + str(len(cnn_score)))


    # calculate quantities from selected data
    nHits = Unwrap(nHits)
    
    print("energy residual")
    energyResidual = Unwrap( (energy - mc_energy)/ mc_energy )
    
    print("angle between beam and daughters")
    beam_angle = Unwrap(SelectionQuantities.BeamTrackShowerAngle(beam_start_pos, beam_end_pos, direction))
    
    print("angle between daughters and mc particle")
    mc_angle = Unwrap(SelectionQuantities.DaughterRecoMCAngle(true_start_pos, true_end_pos, direction))
    
    print("separation")
    pair_separation = Unwrap(SelectionQuantities.ShowerPairSeparation(start_pos))
    
    print("pair angle")
    pair_angle = Unwrap(SelectionQuantities.ShowerPairAngle(start_pos, direction))
    
    print("pair energy")
    pair_energies, pair_leading, pair_second = SelectionQuantities.ShowerPairEnergy(start_pos, energy)
    
    pair_energies = Unwrap(pair_energies) / 1000
    pair_leading = Unwrap(pair_leading) / 1000
    pair_second = Unwrap(pair_second) / 1000
    
    print("Invariant mass")
    inv_mass = Unwrap(SelectionQuantities.InvariantMass(start_pos, direction, energy))
    
    return [cnn_score, nHits, energyResidual, beam_angle, 
            mc_angle, pair_separation, pair_angle, pair_energies, 
            pair_leading, pair_second, inv_mass]


start_time = time.time()
#file_1 = "ROOTFiles/pi0Test_output_merged_100_cal.root"
file_2 = "ROOTFiles/pduneana_Prod4_1GeV_23_4_21_100.root"
file_1 = "ROOTFiles/pi0Test_output_PDSPProd4_MC_1GeV_SCE_DataDriven_reco_100_28_04_21.root"
#file_2 = "ROOTFiles/pduneana_comp.root"

subDirectory = "ComparisonPlots"
os.makedirs(subDirectory + "/", exist_ok=True)

evt_1 = Master.GetData(file_1, "EventID")
evt_2 = Master.GetData(file_2, "event")


print("Are event ID's equivilant? : " + str(np.array_equiv(np.sort(evt_1), np.sort(evt_2))) )

q_1 = ProcessData(file_1)
q_2 = ProcessData(file_2)


print("Saving Plots:")
nbins = 100
l_1 = "pi0"
l_2 = "PDSP"

print("CNN score")
Plots.PlotHistComparison(q_1[0], q_2[0], nbins, alpha=0.5, xlabel="CNN scores of beam daughters", label_1=l_1, label_2=l_2)
save("cnn_score")


print("collection plane hits")
nHits_1 = q_1[1][q_1[1] != 0]
nHits_2 = q_2[1][q_2[1] != 0]

Plots.PlotHistComparison(nHits_1[nHits_1 < 101], nHits_2[nHits_2 < 101], nbins, alpha=0.5, xlabel="Number of collection plane hits", label_1=l_1, label_2=l_2)
save("collection_hits")


print("hits vs mc angle")
plt.figure(2, (12, 5))
plt.subplot(121)
Plots.PlotHist2D(q_1[1], q_1[4], nbins, [-999, 510], [0, np.pi], "Number of collection plane hits", "Angle between shower and MC parent (rad)", title=l_1)

plt.subplot(122)
Plots.PlotHist2D(q_2[1], q_2[4], nbins, [-999, 510], [0, np.pi], "Number of collection plane hits", "Angle between shower and MC parent (rad)", title=l_2)
save("hits_vs_mcAngle")


print("hits vs beam angle")
plt.figure(2, (12, 5))
plt.subplot(121)
Plots.PlotHist2D(q_1[1], q_1[3], nbins, [-999, 501], [0, np.pi], "Number of collection plane hits", "Angle between beam track and daughter shower (rad)", title=l_1)
plt.subplot(122)
Plots.PlotHist2D(q_2[1], q_2[3], nbins, [-999, 501], [0, np.pi], "Number of collection plane hits", "Angle between beam track and daughter shower (rad)", title=l_2)
save("hits_vs_beamAngle")


print("hits vs energy residual")
plt.figure(2, (12, 5))
plt.subplot(121)
Plots.PlotHist2D(q_1[1], q_1[2], nbins, [-999, 501], [-1, 1], "Number of collection plane hits", "Reconstruted energy residual", title=l_1)
plt.subplot(122)
Plots.PlotHist2D(q_2[1], q_2[2], nbins, [-999, 501], [-1, 1], "Number of collection plane hits", "Reconstruted energy residual", title=l_2)
save("hits_vs_energyRes")


print("beam angle")
beam_angle_1 = q_1[3][q_1[3] != -999]
beam_angle_2 = q_2[3][q_2[3] != -999]
Plots.PlotHistComparison(beam_angle_1, beam_angle_2, nbins, alpha=0.5, xlabel="Angle between beam track and daughter shower (rad)", label_1=l_1, label_2=l_2)
save("beam_daughterAngle")


print("pair separation")
pair_separation_1 = q_1[5][q_1[5] > 0]
pair_separation_2 = q_2[5][q_2[5] > 0]
Plots.PlotHistComparison(pair_separation_1[pair_separation_1 < 51], pair_separation_2[pair_separation_2 < 51], nbins, alpha=0.5, xlabel="Shower pair Separation (cm)", label_1=l_1, label_2=l_2)
save("pair_separation")


print("pair angle")
pair_angle_1 = q_1[6][q_1[6] != -999]
pair_angle_2 = q_2[6][q_2[6] != -999]
Plots.PlotHistComparison(pair_angle_1, pair_angle_2, nbins, alpha=0.5, xlabel="Angle between shower pairs (rad)", label_1=l_1, label_2=l_2)
save("shower_pairAngle")


print("pair energies")
pair_leading_1 = q_1[-3][q_1[-3] > 0]
pair_second_1 = q_1[-2][q_1[-2] > 0]

pair_leading_2 = q_2[-3][q_2[-3] > 0]
pair_second_2 = q_2[-2][q_2[-2] > 0]

plt.figure(2, (12, 5))
plt.subplot(121)
Plots.PlotHistComparison(pair_leading_1, pair_second_1, nbins, alpha=0.5, xlabel="Shower energy (GeV)", label_1="shower with the most energy in a pair", label_2="shower with the least energy in a pair", title=l_1)
plt.subplot(122)
Plots.PlotHistComparison(pair_leading_2, pair_second_2, nbins, alpha=0.5, xlabel="Shower energy (GeV)", label_1="shower with the most energy in a pair", label_2="shower with the least energy in a pair", title=l_2)
plt.legend()
plt.tight_layout()
save("daughter_energies")


print("invariant mass")
inv_mass_1 = q_1[-1][q_1[-1] > 0] / 1000
inv_mass_2 = q_2[-1][q_2[-1] > 0] / 1000
Plots.PlotHistComparison(inv_mass_1[inv_mass_1 < 0.5], inv_mass_2[inv_mass_2 < 0.5], nbins, alpha=0.5, xlabel="Shower pair invariant mass (GeV)", sf=3, label_1=l_1, label_2=l_2)
save("inv_mass")

print("done!")
print("ran in %s seconds. " % (time.time() - start_time) )
