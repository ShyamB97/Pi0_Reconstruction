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
from Master import Unwrap, Vector3


def Plot(x, y, xlabel="", ylabel="", title=""):
    """
    Plot line graph.
    """
    plt.plot(x, y, marker="o")
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


def PlotHist(data, bins=100, xlabel="", title="", label="", alpha=1, sf=2):
    """
    Plot histogram of data and axes including bin width.
    ----- Parameters -----
    height      : bin height
    edges       : right edge of the bins
    ----------------------
    """
    height, edges, _ = plt.hist(data, bins, label=label, alpha=alpha)
    binWidth = round((edges[-1] - edges[0]) / len(edges), sf)
    plt.ylabel("Number of events (bin width=" + str(binWidth) + ")")
    plt.xlabel(xlabel)
    plt.title(title)
    plt.tight_layout()
    return height, edges


def PlotHist2D(data_x, data_y, bins=100, x_range=[], y_range=[], xlabel="", ylabel="", title=""):
    """
    Plots two datasets in a 2D histogram.
    """
    # clamp data_x and data_y given the x range
    if len(x_range) == 2:
        data_y = data_y[data_x > x_range[0]] # clamp y before x
        data_x = data_x[data_x > x_range[0]]
        
        data_y = data_y[data_x < x_range[1]]
        data_x = data_x[data_x < x_range[1]]
    
    # clamp data_x and data_y given the y range
    if len(y_range) == 2:
        data_x = data_x[data_y > y_range[0]] # clamp x before y
        data_y = data_y[data_y > y_range[0]]
        
        data_x = data_x[data_y < y_range[1]]
        data_y = data_y[data_y < y_range[1]]

    # plot data with a logarithmic color scale
    plt.hist2d(data_x, data_y, 100, norm=matplotlib.colors.LogNorm())
    plt.colorbar()

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()


def Angle(x_0, y_0, z_0, x_1, y_1, z_1):
    """
    Calculates the angle between two vectors (or list of vectors).
    ----- Parameters -----
    x_n         : vector component x of vector n
    dot         : dot product of the two vectors
    mag_n       : magnitude of vector n
    cos_angle   : cosine of the angle
    ----------------------
    """
    dot = x_0 * x_1 + y_0 * y_1 + z_0 * z_1

    mag_0 = (x_0**2 + y_0**2 + z_0**2)**0.5
    mag_1 = (x_1**2 + y_1**2 + z_1**2)**0.5

    cos_angle = dot / ( mag_0 * mag_1 )
    return np.arccos(cos_angle)


def GetShowerPairs(start_pos):
    """
    Pairs daughter events based on their start vector separation.
    ----- Parameters -----
    array_n     : nth component of start position of each daughter in the event
    dists_n     : nth compnent of the distance between the jth daughter event and all other daughters
    dists       : distance between the jth daughter event and all other daughters
    pair        : index of the daughter pairs per daughter
    pairs       : index of the daughter pairs per event
    evt_pairs   : index of the daughter pairs for all events
    ----------------------
    """
    nEvents = len(start_pos.x) # should be the same for all
    
    # get shower pairs by comparing the distances of each daughter per event
    evt_pairs = []
    for i in range(nEvents):
        array_x = start_pos.x[i]
        array_y = start_pos.y[i]
        array_z = start_pos.z[i]
    
        # remove null data
        array_x = array_x[array_x != -999]
        array_y = array_y[array_y != -999]
        array_z = array_z[array_z != -999]

        pairs = []
        if len(array_x) > 1:
            for j in range(len(array_x)):
                
                #calculate per displacement per direction
                dists_x = array_x[j] - array_x
                dists_y = array_y[j] - array_y
                dists_z = array_z[j] - array_z
                
                dists = np.sqrt( dists_x**2 + dists_y**2 + dists_z**2 ) # distance beteween shower
                
                # distance to self is zero, so set these to something large to avoid self paring
                dists = np.where( dists <= 0, 1E10, dists )
                
                # pair this shower to the one closest to it (index wise)
                pair = [j, np.argmin(dists)]
                pairs.append(pair)
        evt_pairs.append(pairs)


    # remove output of showers with no pairs and duplicates
    evt_pairs = [removeDuplicates(i) for i in evt_pairs]
    return np.array(evt_pairs, object)


def ShowerPairSeparation(start_pos, actual_pairs=True):
    """
    Gets separation of paired daughter events, either all possible pairs or 
    the true pairs (closest separation per daughter).
    ----- Parameters -----
    nEvents         : number of events
    null_value      : value to replace zero when calculating distances, done to exclude distance to self
    array_n         : nth component of start position of each daughter in the event
    dists_n         : nth compnent of the distance between the jth daughter event and all other daughters
    dists           : distance between the jth daughter event and all other daughters
    pair_dists      : distances per daughter pair
    evt_pairs_dist  : distances of the daughter pairs for all events
    x_0, ...        : start positions of the first shower pair
    x_1, ...        : start positions of the second shower pair
    showerPairs     : daughters paired using GetSHowerPairs
    ----------------------
    """
    nEvents = len(start_pos.x) # should be the same for all
    
    null_value = 1E10 # value to set for null data i.e. distances of zero
    
    if actual_pairs is False:
        # gets seperation of all possible shower pairs
        evt_pairs_dist = []
        for i in range(nEvents):
            array_x = start_pos.x[i]
            array_y = start_pos.y[i]
            array_z = start_pos.z[i]
        
            # remove null data
            array_x = array_x[array_x != -999]
            array_y = array_y[array_y != -999]
            array_z = array_z[array_z != -999]
        
            pair_dists = []
            if len(array_x) > 1:
                for j in range(len(array_x)):
                    
                    #calculate per displacement per direction
                    dists_x = array_x[j] - array_x
                    dists_y = array_y[j] - array_y
                    dists_z = array_z[j] - array_z
                    
                    dists = np.sqrt( dists_x**2 + dists_y**2 + dists_z**2 ) # distance beteween shower
    
                    # distance to self is zero, so set these to something large to avoid self paring
                    dists = np.where( dists <= 0, 1E10, dists )
    
                    pair_dists.append(dists)
                    
            evt_pairs_dist.append(pair_dists)
        
        evt_pairs_dist = [removeDuplicates(i) for i in evt_pairs_dist]
        
        # remove null values i.e distance to self
        for i in range(len(evt_pairs_dist)):
            evt = evt_pairs_dist[i]
            for j in range(len(evt)):
                evt[j] = evt[j][evt[j] != null_value]
            evt_pairs_dist[i] = evt

    else:
        evt_pairs_dist = []
        showerPairs = GetShowerPairs(start_pos)
        for i in range(len(showerPairs)):
            pair_dists = []
            if(len(showerPairs[i]) == 0):
                pair_dists.append(-999)
            else:
                for pair in showerPairs[i]:
                    x_0 = start_pos.x[i][pair[0]]
                    x_1 = start_pos.x[i][pair[1]]
                    y_0 = start_pos.y[i][pair[0]]
                    y_1 = start_pos.y[i][pair[1]]
                    z_0 = start_pos.z[i][pair[0]]
                    z_1 = start_pos.z[i][pair[1]]
                    
                    x = x_1 - x_0
                    y = y_1 - y_0
                    z = z_1 - z_0
                    
                    pair_dists.append( np.sqrt( x**2 + y**2 + z**2 ) )

            evt_pairs_dist.append(pair_dists)

    return np.array(evt_pairs_dist, object)


def ShowerPairAngle(start_pos, _dir):
    """
    Gets the angle between daughter pairs.
    ----- Parameters -----
    showerPairs     : daughters paired using GetSHowerPairs
    x_0, ...        : direction vector components of the first shower pair
    x_1, ...        : direction vector components of the second shower pair
    evt_angles      : angles of shower per event
    ----------------------
    """
    showerPairs = GetShowerPairs(start_pos)

    # get the angle between the shower pairs
    angles = []
    for i in range(len(showerPairs)):
        evt_angles = []
        if(len(showerPairs[i]) != 0):
            for pair in showerPairs[i]:
                x_0 = _dir.x[i][pair[0]]
                x_1 = _dir.x[i][pair[1]]
                y_0 = _dir.y[i][pair[0]]
                y_1 = _dir.y[i][pair[1]]
                z_0 = _dir.z[i][pair[0]]
                z_1 = _dir.z[i][pair[1]]

                evt_angles.append(Angle(x_0, y_0, z_0, x_1, y_1, z_1))
        angles.append(np.array(evt_angles))
    return np.array(angles, object)


def ShowerPairEnergy(start_pos, shower_energy):
    """
    Gets the Energy of the daughter pairs.
    ----- Parameters -----
    pairs           : daughters paired using GetSHowerPairs
    energy          : list of energy for pair 0 and pair 1
    pairs_energy    : shower energy paired by daughters for each event
    leading_energy  : shower in pair with the higher energy
    secondary_energy: shower in pair with the lower energy
    energy_min      : shower in pair with the higher energy
    energy_min      : shower in pair with the lower energy
    ----------------------
    """
    pairs = GetShowerPairs(start_pos)

    pairs_energy = []
    leading_energy = []
    secondary_energy = []
    for i in range(len(pairs)):
        evt_pairs = pairs[i]
        evt_energy = shower_energy[i]
    
        energy = []
        energy_min = []
        energy_max = []
        for j in range(len(evt_pairs)):
            pair = evt_pairs[j]
            paired_energy = [evt_energy[pair[0]],  evt_energy[pair[1]]]
 
            energy.append(paired_energy)
            energy_min.append(min(paired_energy))
            energy_max.append(max(paired_energy))

        pairs_energy.append(energy)
        leading_energy.append(energy_max)
        secondary_energy.append(energy_min)
    return np.array(pairs_energy, object), np.array(leading_energy, object), np.array(secondary_energy, object)


def removeDuplicates(pairs_list):
    """
    Used to remove identical daugher pairs and daughter pairs without identical copies
    (daugthers that are correctly paired from GetShowerPair should have two lists that are
     reveresed i.e. [0, 3] and [3, 0]).
    ----- Parameters -----
    _list   : returned list of pairs
    copy    : list to keep track of pairs already found
    ----------------------
    """
    _list = []
    copy = []
    for i in range(len(pairs_list)):
        for j in range(len(pairs_list)):
            if pairs_list[i][0] == pairs_list[j][1] and pairs_list[i][1] == pairs_list[j][0]:
                if j not in copy:
                    _list.append(pairs_list[j])
                    #print("pairs found")
                copy.append(i)
    return _list


def BeamTrackShowerAngle(beam_start, beam_end, daughter_dir):
    """
    Calculates the angle between the overall beam track direction and the 
    direction of its daughters.
    ----- Parameters -----
    beam_dists_x, ...   : component of beam distances travelled
    beam_dists          : beam distances travelled
    evt_x, ...          : direction vector compnents of daughters per event
    angle               : angle between beam track and daughter directions
    ----------------------
    """
    beam_dist_x = beam_end.x - beam_start.x
    beam_dist_y = beam_end.y - beam_start.y
    beam_dist_z = beam_end.z - beam_start.z


    beam_dist = ( beam_dist_x**2 + beam_dist_y**2 + beam_dist_z**2 )**0.5

    angles = []
    for i in range(len(daughter_dir.x)):
        evt_x = daughter_dir.x[i]
        evt_y = daughter_dir.y[i]
        evt_z = daughter_dir.z[i]

        angle = []
        if beam_dist[i] != 0:
            for j in range(len(evt_x)):
                if evt_x[j] == -999:
                    angle.append(-999)
                else:
                    angle.append( Angle(evt_x[j], evt_y[j], evt_z[j], beam_dist_x[i], beam_dist_y[i], beam_dist_z[i]) )
        else:
            # add nulls for each daughter in events with no beam
            num_daughters = len( evt_x )
            angle = [-999] * num_daughters

        angles.append(angle)

    return np.array(angles, object)


def DaughterRecoMCAngle(true_start, true_end, _dir):
    """
    Reconstrcut the angle between the reco daughter particle and the true MC particle.
    ----- Parameters -----
    true_dist_x, ...   : components of the distance of the MC particle
    angle              : recoMC angle per daughter in an event
    angles             : recoMC angles per event
    ----------------------
    """
    true_dist_x = true_end.x - true_start.x
    true_dist_y = true_end.y - true_start.y
    true_dist_z = true_end.z - true_start.z

    angles = []
    for i in range(len(true_start.x)):
        angle = []
        for j in range(len(true_start.x[i])):    
            if true_start.x[i][j] == -999:
                angle.append(-999)
            else:
                angle.append( Angle(_dir.x[i][j], 
                                    _dir.y[i][j], 
                                    _dir.z[i][j], 
                                    true_dist_x[i][j], 
                                    true_dist_y[i][j], 
                                    true_dist_z[i][j]) )
        angles.append(angle)
    return np.array(angles, object)


def CNNScore(em, track, average=False):
    """
    Calculate the em/track like CNN score per daughter or per event by taking
    the average.
    ----- Parameters -----
    cnnScore    : normalised em/shower like score
    ----------------------
    """
    cnnScore = em / (em + track)
    if average is True:
        cnnScore = np.array([np.nanmean(evt) for evt in cnnScore])

    return cnnScore


def ResidualShowerEnergy(shower_energy, true_energy):
    """
    Calculate the shower energy residual per daughter.
    ----- Parameters -----
    residual_energy    : energy residual per event per shower
    evt                : ennergy of daughters in an event
    evt_res            : energy residual per shower
    ----------------------
    """
    residual_energy = []
    for i in range(len(shower_energy)):
        evt = shower_energy[i]
        evt_res = []
        for j in range(len(evt)):
            if(evt[j] != -999):
                evt_res.append( (evt[j] - true_energy[i][j]) / true_energy[i][j] )
            else:
                evt_res.append(-999)
        residual_energy.append(evt_res)
    
    return np.array(residual_energy, object)


def InvariantMass(start_pos, _dir, shower_energy):
    """
    Calculate the invariant mass from shower pairs.
    ----- Parameters -----
    showerPair_angles       : angle between shower pairs
    showerPair_energy       : energy of the showers forming the pairs
    inv_mass                : invariant mass per event per shower pair
    evt_angle               : angle between shower pairs for an event
    showerPair_energy       : energy of the showers forming the pair for an event
    angle                   : opening angle for a pair
    energy                  : shower pair energies
    ----------------------
    """    
    showerPair_angles = ShowerPairAngle(start_pos, _dir)
    showerPair_energy, _, _ = ShowerPairEnergy(start_pos, shower_energy)
    
    inv_mass = []
    for i in range(len(showerPair_angles)):
        evt_angle = showerPair_angles[i]
        evt_energy = showerPair_energy[i]

        evt_inv_mass = []
        for j in range(len(evt_angle)):
            angle = evt_angle[j]
            energies = evt_energy[j]

            if angle != -999:
                if energies[0] != -999 and energies[1] != -999:
                    # calculate the invariant mass using particle kinematics
                    evt_inv_mass.append( np.sqrt(2 * energies[0] * energies[1] * (1 - np.cos(angle))) )
                else:
                    evt_inv_mass.append(-999)
            else:
                evt_inv_mass.append(-999)

        inv_mass.append(evt_inv_mass)

    return np.array(inv_mass, object)


def GetNumResonantParticles(start_pos):
    """
    Gets the number of pairs i.e. the number of particles which produce these daughters,
    assuming it decays only into two shower pairs.
    ----- Parameters -----
    numParticles    : number of pairs in an event
    ----------------------
    """
    pairs = GetShowerPairs(start_pos)

    numParticles = []
    for i in range(len(pairs)):
        numParticles.append(len(pairs[i]))
    return np.array(numParticles, object)


def GetResonantMomentum(start_pos, shower_energy):
    """
    Returns the monemtum of the parent per shower pair per event, assuming the
    daugthers are photons.
    ----- Parameters -----
    energies    : shower pair energies
    momenta     : resonant momenta
    evt         : energy pairs per event
    evt_mom     : resontant momenta per event
    e0, e1      : shower energies in a pair
    ----------------------
    """
    energies, _, _ = ShowerPairEnergy(start_pos, shower_energy)
    
    momenta = []
    for i in range(len(energies)):
        evt = energies[i]
        
        evt_mom = []
        for j in range(len(evt)):
            
            e0 = evt[j][0]
            e1 = evt[j][1]
            
            if evt[j] != -999:
                if e0 != -999 and e1 != -999:
                    evt_mom.append(e0 + e1)
                else:
                    evt_mom.append(-999)
            else:
                evt_mom.append(-999)

        momenta.append(evt_mom)
    return np.array(momenta, object)


def GetResonanceEnergy(start_pos, _dir, shower_energy):
    """
    Calculate the emergy of the parent particle produding the shower pairs.
    ----- Parameters -----
    mass        : Invariant mass calculated from shower pair energies
    momentum    : momentum of the parent particle assuming the daughters are massless
    res_mass    : parent mass of a pair for an event
    res_mom     : parent momentum of a pair for an event
    evt_energy  : parent energy for an event 
    energy      : parent energies per event
    ----------------------
    """
    mass = InvariantMass(start_pos, _dir, shower_energy)
    momentum = GetResonantMomentum(start_pos, _dir, shower_energy)
    
    energy = []
    for i in range(len(mass)):
        evt_energy = []
        for j in range(len(mass[i])):
            
            res_mass = mass[i][j]
            res_mom = momentum[i][j]
            
            if res_mass != -999 and res_mom != -999:
                evt_energy.append( ( res_mass**2 + res_mom**2 )**0.5 )
            else:
                evt_energy.append(-999)
        energy.append(evt_energy)
    return np.array(energy, object)


def GetSHowerStartHits(hit_radial, hit_longitudinal, geometry=[1, -1, 4]):
    """
    Will return the number of hits within a cylindrical geometry about the shower
    start point.
    ----- Parameters -----
    geometry        : cylinder bounds,
                      [radius,
                       legnth along shower direction behind the start point,
                       legnth along shower direction in front the start point]
    start_hits      : number of hits within the geometry of each shower
    evt_r           : radial value to compare for each shower per event
    evt_l           : longitudinal value to compare for each shower per event
    shower_r        : radial value to compare for a shower
    shower_l        : longitudinal value to compare for a shower
    hits            : number of hits within the geometry for a shower
    ----------------------
    """
    start_hits = []
    for i in range(len(hit_radial)):
        # per event
        evt_r = hit_radial[i]
        evt_l = hit_longitudinal[i]
        evt_hits = []
        for j in range(len(evt_r)):
            # per shower
            shower_r = evt_r[j]
            shower_l = evt_l[j]
            hits = 0
            for k in range(len(shower_r)):
                # per hit
                if shower_r[k] < geometry[0] and geometry[1] < shower_l[k] < geometry[2]:
                    # hit near the start!
                    hits += 1
            evt_hits.append(hits)
        start_hits.append(evt_hits)


#data = Master.Data("pduneana_Prod4_1GeV_2_9_21.root")
data = Master.Data("pi0Test_output.root")

r = 1 # cm
l_min = -1 # cm
l_max = 4 # cm

"""SELECTION"""

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


hit_radial = data.hit_radial()
hit_longitudinal = data.hit_longitudinal()

"""
#Calculated quantities
print("calculating cnn_score...")
cnn_score = CNNScore(data.cnn_em(), data.cnn_track())
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


print("energy residual")
energyResidual = Unwrap( (energy - mc_energy)/ mc_energy )

print("angle between beam and daughters")
beam_angle = Unwrap(BeamTrackShowerAngle(beam_start_pos, beam_end_pos, direction))

print("angle between daughters and mc particle")
mc_angle = Unwrap(DaughterRecoMCAngle(true_start_pos, true_end_pos, direction))

print("separation")
pair_separation = Unwrap(ShowerPairSeparation(start_pos))

print("pair angle")
pair_angle = Unwrap(ShowerPairAngle(start_pos, direction))

print("pair energy")
pair_energies, pair_leading, pair_second = ShowerPairEnergy(start_pos, energy)

pair_energies = Unwrap(pair_energies) / 1000
pair_leading = Unwrap(pair_leading) / 1000
pair_second = Unwrap(pair_second) / 1000

print("Invariant mass")
inv_mass = Unwrap(InvariantMass(start_pos, direction, energy))
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
inv_mass = inv_mass[inv_mass > 0] / 1000
PlotHist(inv_mass[inv_mass < 0.5], 100, "Shower pair invariant mass (GeV)", "", 3)
"""

"""
pair_angle = pair_angle[pair_angle != -999]
PlotHist(pair_angle, 100, "Angle between beam track and daughter shower (rad)")
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

"""
pair_leading = pair_leading[pair_leading > 0]
pair_second = pair_second[pair_second > 0]

#plt.figure(2, (10, 5))

#plt.subplot(122)
_, edges = PlotHist(pair_second, title="(b)", label="shower with the least energy in a pair", alpha=0.5)

#plt.subplot(121)
PlotHist(pair_leading, edges, xlabel="Shower energy (GeV)", title="(a)", label="shower with the most energy in a pair", alpha=0.5)

plt.legend()
plt.tight_layout()
"""

