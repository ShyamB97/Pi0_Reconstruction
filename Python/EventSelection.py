#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 09:49:07 2021

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
from Plots import Plot, PlotHist, PlotHist2D

def save(name):
    out = subDirectory + name + reference_filename + ".png"
    plt.savefig(out)
    plt.close()


# handle command line imports
args = sys.argv

if "--help" in args or len(args) == 1:
    print("usage: --file <root file> --outName <file name to append> --outDir <output directory>")
    exit(0)
if "--file" in args:
    root_filename = args[args.index("--file") + 1]
    print(root_filename)
else:
    print("no file chosen!")
    exit(1)
if "--outName" in args:
    reference_filename = args[args.index("--outName") + 1]
else:
    reference_filename = "_" + root_filename[19:-5]
if "--outDir" in args:
    subDirectory = args[args.index("--outDir") + 1] + "/ "
else:
    subDirectory = root_filename[19:-5] + "/"


start_time = time.time()
print("Opening root file: " + root_filename)
os.makedirs(subDirectory, exist_ok=True)
data = Master.Data(root_filename)


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


# create the mask
mask = Master.SelectionMask()
# initilise mask (set the shape and set all values to 1)
mask.InitiliseMask(nHits)


### Hacked CNN score fix ###
mask999 = Master.SelectionMask()
mask999.InitiliseMask(cnn_em)
mask999.CutMask(cnn_em, -999, Conditional.NOT_EQUAL)

cnn_em = mask999.Apply(cnn_em)
cnn_track = mask999.Apply(cnn_track)

cnn_score = SelectionQuantities.CNNScore(cnn_em, cnn_track)

"""
# select data based on cnn score
mask.CutMask(cnn_score, 0.6, Conditional.GREATER)


# cut the data using the mask
start_pos = mask.Apply(start_pos)

direction = mask.Apply(direction)

energy = mask.Apply(energy)
mc_energy = mask.Apply(mc_energy)
nHits = mask.Apply(nHits)

beam_start_pos = mask.Apply(beam_start_pos, True)
beam_end_pos = mask.Apply(beam_end_pos, True)

true_start_pos = mask.Apply(true_start_pos)
true_end_pos = mask.Apply(true_end_pos)
"""

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


# save plots of data
print("Saving Plots:")
nbins = 100

print("CNN score")
cnn_score = Unwrap(cnn_score)
PlotHist(cnn_score, nbins, "CNN scores of beam daughters")
save("cnn_score")


print("collection plane hits")
nHits_plot = nHits[nHits != 0]
PlotHist(nHits_plot[nHits_plot < 101], nbins, "Number of collection plane hits")
save("collection_hits")


print("hits vs mc angle")
PlotHist2D(nHits, mc_angle, nbins, [-999, 510], [0, np.pi], "Number of collection plane hits", "Angle between shower and MC parent (rad)")
save("hits_vs_mcAngle")


print("hits vs beam angle")
PlotHist2D(nHits, beam_angle, nbins, [-999, 501], [0, np.pi], "Number of collection plane hits", "Angle between beam track and daughter shower (rad)")
save("hits_vs_beamAngle")


print("hits vs energy residual")
PlotHist2D(nHits, energyResidual, nbins, [-999, 501], [-1, 1], "Number of collection plane hits", "Reconstruted energy residual")
save("hits_vs_energyRes")


print("beam angle")
beam_angle = beam_angle[beam_angle != -999]
PlotHist(beam_angle, nbins, "Angle between beam track and daughter shower (rad)")
save("beam_daughterAngle")


print("pair separation")
pair_separation = pair_separation[pair_separation > 0]
PlotHist(pair_separation[pair_separation < 51], nbins, "Shower pair Separation (cm)")
save("pair_separation")


print("pair angle")
pair_angle = pair_angle[pair_angle != -999]
PlotHist(pair_angle, nbins, "Angle between shower pairs (rad)")
save("shower_pairAngle")


print("pair energies")
pair_leading = pair_leading[pair_leading > 0]
pair_second = pair_second[pair_second > 0]
_, edges = PlotHist(pair_second, label="shower with the least energy in a pair", alpha=0.5)
PlotHist(pair_leading, edges, xlabel="Shower energy (GeV)", label="shower with the most energy in a pair", alpha=0.5)
plt.legend()
plt.tight_layout()
save("daughter_energies")


print("invariant mass")
inv_mass = inv_mass[inv_mass > 0] / 1000
PlotHist(inv_mass[inv_mass < 0.5], nbins, "Shower pair invariant mass (GeV)", "", 3)
save("inv_mass")

print("done!")
print("ran in %s seconds. " % (time.time() - start_time) )