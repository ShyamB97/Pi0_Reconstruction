#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 09:49:07 2021

@author: sb16165
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import Master
from Master import Unwrap, Vector3, Conditional
from Plots import PlotHist, PlotHist2D, PlotHistComparison
import SelectionQuantities
import Plots


#data = Master.Data("ROOTFiles/pduneana_Prod4_1GeV_2_9_21.root")
#data = Master.Data("pi0Test_output_PDSPProd4_MC_1GeV_SCE_DataDriven_reco_1K_1_24_03_21.root")
#data = Master.Data("ROOTFiles/pi0Test_output_PDSPProd4_MC_1GeV_SCE_DataDriven_reco_3p5K_29_04_21.root");
#file = "ROOTFiles/pi0_0p5GeV_100K_5_7_21.root"
file = "ROOTFiles/pi0Test_output.root"
data = Master.Data(file)

evd_ID = Master.GetData(file, "EventID")

r = 1 # cm
l_min = -1 # cm
l_max = 4 # cm

beamData = False
cutData = False
param = 0
cut = 0
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
    #pair_separation = np.reshape(pair_separation, (len(nHits)))
    
    print("pair angle")
    pair_angle = SelectionQuantities.ShowerPairAngle(shower_pairs, direction)
    
    #print("true pair angle")
    #true_dir = true_momentum.Normalise()
    #true_pair_angle = SelectionQuantities.ShowerPairAngle(shower_pairs, true_dir)
    
    print("pair energy")
    pair_energies, pair_leading, pair_second = SelectionQuantities.ShowerPairEnergy(shower_pairs, energy)
    
    print("Invariant mass")
    inv_mass = SelectionQuantities.InvariantMass(pair_angle, pair_energies)
    
    parameters_pair = [pair_separation, pair_angle, pair_leading, pair_second, inv_mass, pair_energies]

    parameters = [Unwrap(p) for p in parameters]
    #parameters_pair = [Unwrap(p) for p in parameters_pair]

    return parameters, parameters_pair


print("getting daughter info...")
start_pos = data.start_pos()
direction = data.direction()
energy = data.energy()
mc_energy = data.mc_energy()
nHits = data.nHits()

print("getting beam info...")
beam_start_pos = data.beam_start_pos()
beam_end_pos = data.beam_end_pos()

print("getting daughter truth info...")
true_start_pos = data.true_start_pos()
true_end_pos = data.true_end_pos()

# custom data
hit_radial = data.hit_radial()
hit_longitudinal = data.hit_longitudinal()
pandoraTag = data.pandoraTag()


#Calculated quantities

print("calculating cnn_score...")
cnn_em = data.cnn_em()
cnn_track = data.cnn_track()
cnn_score = SelectionQuantities.CNNScore(cnn_em, cnn_track)


"""
# apply selection to all the data

cnnSelection = Master.Selection(cnn_score, 0.6, True)

print("Applying selection...")
start_pos = cnnSelection.Apply(start_pos)

direction = cnnSelection.Apply(direction)

energy = cnnSelection.Apply(energy)
mc_energy = cnnSelection.Apply(mc_energy)
nHits = Unwrap(cnnSelection.Apply(nHits))

beam_start_pos = cnnSelection.Apply(beam_start_pos, True)
beam_end_pos = cnnSelection.Apply(beam_end_pos, True)

true_start_pos = cnnSelection.Apply(true_start_pos)
true_end_pos = cnnSelection.Apply(true_end_pos)
"""


print("energy residual")
energyResidual = Unwrap( (energy - mc_energy)/ mc_energy)

#print("angle between beam and daughters")
#beam_angle = Unwrap(SelectionQuantities.BeamTrackShowerAngle(beam_start_pos, beam_end_pos, direction))

print("angle between daughters and mc particle")
mc_angle = Unwrap(SelectionQuantities.DaughterRecoMCAngle(true_start_pos, true_end_pos, direction))

print("shower pairs")
shower_pairs = SelectionQuantities.GetShowerPairs(start_pos)

print("separation")
pair_separation = SelectionQuantities.ShowerPairSeparation(start_pos, shower_pairs)

print("pair angle")
pair_angle = SelectionQuantities.ShowerPairAngle(shower_pairs, direction)

print("pair energy")
pair_energies, pair_leading, pair_second = SelectionQuantities.ShowerPairEnergy(shower_pairs, energy)

#pair_energies = Unwrap(pair_energies)
#pair_leading = Unwrap(pair_leading)
#pair_second = Unwrap(pair_second)

print("start hits")
start_hits = SelectionQuantities.GetShowerStartHits(hit_radial, hit_longitudinal)

print("Invariant mass")
inv_mass = SelectionQuantities.InvariantMass(pair_angle, pair_energies)

#mask = Master.SelectionMask()
#mask.InitiliseMask(pair_angle)
#mask.CutMask(pair_angle, 1.5, Conditional.LESS)

#pairs_cut = mask.Apply(shower_pairs)
#inv_mass = mask.Apply(inv_mass)

#nHits = SelectionQuantities.GetShowerPairValues(nHits, pairs_cut)

#PlotHist2D(Unwrap(inv_mass), Unwrap(nHits), 100)

#inv_mass = Unwrap(inv_mass)
#PlotHist(inv_mass, 100)

"""
start_hits = Unwrap(start_hits)
PlotHist(start_hits, 100, "Shower start hits")
"""


"""
pandoraTag = Unwrap(pandoraTag)
unique, amount = np.unique(pandoraTag, return_counts=True)
unique = list(unique)
unique[0] = str(unique[0]) + "\nN/A"
unique[1] = str(unique[1]) + "\nshower"
unique[2] = str(unique[2]) + "\ntrack"
plt.bar(unique, amount)
"""

"""
cnn_score = Unwrap(cnn_score)
PlotHist(cnn_score[cnn_score != -999], 100, "CNN scores of beam daughters")
"""

"""
pair_separation = pair_separation[pair_separation > 0]
PlotHist(pair_separation[pair_separation < 51], 100, "Shower pair Separation (cm)")
"""

"""
PlotHist2D(nHits, mc_angle, 100, [-999, 510], [0, np.pi], "Number of collection plane hits", "Angle between shower and MC parent (rad)")
"""

"""
PlotHist2D(nHits, beam_angle, 100, [-999, 501], [0, np.pi], "Number of collection plane hits", "Angle between beam track and daughter shower (rad)")
"""

"""
PlotHist2D(nHits, energyResidual, 100, [-999, 501], [-1, 1], "Number of collection plane hits", "Reconstruted energy residual")
"""

"""
inv_mass = inv_mass[inv_mass > 0]
PlotHist(inv_mass[inv_mass < 0.5], 100, "Shower pair invariant mass (GeV)", "", 3)
"""

"""
pair_angle = pair_angle[pair_angle != -999]
PlotHist(pair_angle, 100, "Angle between shower pairs (rad)")
"""

"""
beam_angle = beam_angle[beam_angle != -999]
PlotHist(beam_angle, 100, "Angle between beam track and daughter shower (rad)")
"""

"""
nHits = Unwrap(nHits)
nHits = nHits[nHits != 0]
PlotHist(nHits[nHits < 101], 100, "Number of collection plane hits")
"""


"""pair_leading = pair_leading[pair_leading > 0]
pair_second = pair_second[pair_second > 0]
PlotHistComparison(pair_leading, pair_second, 100, xlabel="Shower energy (GeV)", label_1="shower with the most energy in a pair", label_2="shower with the least energy in a pair", alpha=0.5, density=False)
"""