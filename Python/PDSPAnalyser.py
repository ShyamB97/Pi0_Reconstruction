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


def PlotHist(data, bins=100, xlabel="", title="", sf=2):
    """
    Plot histogram of data and axes including bin width.
    ----- Parameters -----
    height      : bin height
    edges       : right edge of the bins
    ----------------------
    """
    height, edges, _ = plt.hist(data, bins)
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


def GetShowerPairs(start_x, start_y, start_z):
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
    nEvents = len(start_x) # should be the same for all
    
    # get shower pairs by comparing the distances of each daughter per event
    evt_pairs = []
    for i in range(nEvents):
        array_x = start_x[i]
        array_y = start_y[i]
        array_z = start_z[i]
    
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


def ShowerPairSeparation(start_x, start_y, start_z, actual_pairs=True):
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
    nEvents = len(start_x) # should be the same for all
    
    null_value = 1E10 # value to set for null data i.e. distances of zero
    
    if actual_pairs is False:
        # gets seperation of all possible shower pairs
        evt_pairs_dist = []
        for i in range(nEvents):
            array_x = start_x[i]
            array_y = start_y[i]
            array_z = start_z[i]
        
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
        showerPairs = GetShowerPairs(start_x, start_y, start_z)
        for i in range(len(showerPairs)):
            pair_dists = []
            if(len(showerPairs[i]) == 0):
                pair_dists.append(-999)
            else:
                for pair in showerPairs[i]:
                    x_0 = start_x[i][pair[0]]
                    x_1 = start_x[i][pair[1]]
                    y_0 = start_y[i][pair[0]]
                    y_1 = start_y[i][pair[1]]
                    z_0 = start_z[i][pair[0]]
                    z_1 = start_z[i][pair[1]]
                    
                    x = x_1 - x_0
                    y = y_1 - y_0
                    z = z_1 - z_0
                    
                    pair_dists.append( np.sqrt( x**2 + y**2 + z**2 ) )

            evt_pairs_dist.append(pair_dists)

    return np.array(evt_pairs_dist, object)


def ShowerPairAngle(start_x, start_y, start_z, dir_x, dir_y, dir_z):
    """
    Gets the angle between daughter pairs.
    ----- Parameters -----
    showerPairs     : daughters paired using GetSHowerPairs
    x_0, ...        : direction vector components of the first shower pair
    x_1, ...        : direction vector components of the second shower pair
    evt_angles      : angles of shower per event
    ----------------------
    """
    showerPairs = GetShowerPairs(start_x, start_y, start_z)    

    # get the angle between the shower pairs
    angles = []
    for i in range(len(showerPairs)):
        evt_angles = []
        if(len(showerPairs[i]) != 0):
            for pair in showerPairs[i]:
                x_0 = dir_x[i][pair[0]]
                x_1 = dir_x[i][pair[1]]
                y_0 = dir_y[i][pair[0]]
                y_1 = dir_y[i][pair[1]]
                z_0 = dir_z[i][pair[0]]
                z_1 = dir_z[i][pair[1]]

                evt_angles.append(Angle(x_0, y_0, z_0, x_1, y_1, z_1))
        angles.append(np.array(evt_angles))
    return np.array(angles, object)


def ShowerPairEnergy(start_x, start_y, start_z, shower_energy):
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
    pairs = GetShowerPairs(start_x, start_y, start_z)

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


def BeamTrackShowerAngle(start_x, start_y, start_z, end_x, end_y, end_z, dir_x, dir_y, dir_z):
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
    beam_dist_x = end_x - start_x
    beam_dist_y = end_y - start_y
    beam_dist_z = end_z - start_z

    #beam_dist = np.sqrt(beam_dist_x**2 + beam_dist_y**2 + beam_dist_z**2)
    beam_dist = ( beam_dist_x**2 + beam_dist_y**2 + beam_dist_z**2 )**0.5

    angles = []
    for i in range(len(dir_x)):
        evt_x = dir_x[i]
        evt_y = dir_y[i]
        evt_z = dir_z[i]

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


def DaughterRecoMCAngle(true_start_x, true_start_y, true_start_z, true_end_x, true_end_y, true_end_z, dir_x, dir_y, dir_z):
    """
    Reconstrcut the angle between the reco daughter particle and the true MC particle.
    ----- Parameters -----
    true_dist_x, ...   : components of the distance of the MC particle
    angle              : recoMC angle per daughter in an event
    angles             : recoMC angles per event
    ----------------------
    """
    true_dist_x = true_end_x - true_start_x
    true_dist_y = true_end_y - true_start_y
    true_dist_z = true_end_z - true_start_z

    angles = []
    for i in range(len(true_start_x)):
        angle = []
        for j in range(len(true_start_x[i])):    
            if true_start_x[i][j] == -999:
                angle.append(-999)
            else:
                angle.append( Angle(dir_x[i][j], 
                                    dir_y[i][j], 
                                    dir_z[i][j], 
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


def InvariantMass(start_x, start_y, start_z, dir_x, dir_y, dir_z, shower_energy):
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
    showerPair_angles = ShowerPairAngle(start_x, start_y, start_z, dir_x, dir_y, dir_z)
    showerPair_energy, _, _ = ShowerPairEnergy(start_x, start_y, start_z, shower_energy)
    
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


def GetNumResonantParticles(start_x, start_y, start_z):
    """
    Gets the number of pairs i.e. the number of particles which produce these daughters,
    assuming it decays only into two shower pairs.
    ----- Parameters -----
    numParticles    : number of pairs in an event
    ----------------------
    """
    pairs = GetShowerPairs(start_x, start_y, start_z)

    numParticles = []
    for i in range(len(pairs)):
        numParticles.append(len(pairs[i]))
    return np.array(numParticles, object)


def GetResonantMomentum(start_x, start_y, start_z, shower_energy):
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
    energies, _, _ = ShowerPairEnergy(start_x, start_y, start_z, shower_energy)
    
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


def GetResonanceEnergy(start_x, start_y, start_z, dir_x, dir_y, dir_z, shower_energy):
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
    mass = InvariantMass(start_x, start_y, start_z, dir_x, dir_y, dir_z, shower_energy)
    momentum = GetResonantMomentum(start_x, start_y, start_z, shower_energy)
    
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


def Selection(data, value, cut, greater=True, beamData=True):
    """
    Selection function. Applies a cut wrt the value on data. Either greater than
    or less than.
    ----- Parameters -----
    selected_data       : data post selection
    evt                 : data per event
    new_evt             : selected data per event
    beamData            : if the data is event level only i.e. not daughter data,
                          important since beam data and daughter data have different
                          struture in the per event data.
    ----------------------
    """
    selected_data = []
    for i in range(len(value)):
        evt = data[i]
        evt_value = value[i]
        new_evt = []
        for j in range(len(evt_value)):

            if greater is True and evt_value[j] > cut:
                if beamData is True:
                    selected_data.append(evt)
                    break
                else:
                    new_evt.append(evt[j])

            if greater is False and evt_value[j] < cut:
                if beamData is True:
                    selected_data.append(evt)
                    break
                else:
                    new_evt.append(evt[j])

        if len(new_evt) > 0:
            selected_data.append(new_evt)

    return np.array(selected_data, object)



data = Master.Data("pduneana_Prod4_1GeV_2_9_21.root")


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

#Calculated quantities
print("calculating cnn_score...")
cnn_score = CNNScore(data.cnn_em(), data.cnn_track())
cut = 0.6
# apply selection to all the data

print("Applying selection...")
start_pos = Vector3(Selection(start_pos.x, cnn_score, cut),
                    Selection(start_pos.y, cnn_score, cut),
                    Selection(start_pos.z, cnn_score, cut))

direction = Vector3(Selection(direction.x, cnn_score, cut),
                    Selection(direction.y, cnn_score, cut),
                    Selection(direction.z, cnn_score, cut))

energy = Selection(energy, cnn_score, cut)
mc_energy = Selection(mc_energy, cnn_score, cut)
nHits = Unwrap(Selection(nHits, cnn_score, cut))

beam_start_pos = Vector3(Selection(beam_start_pos.x, cnn_score, cut, beamData=True),
                         Selection(beam_start_pos.y, cnn_score, cut, beamData=True),
                         Selection(beam_start_pos.z, cnn_score, cut, beamData=True))

beam_end_pos = Vector3(Selection(beam_end_pos.x, cnn_score, cut, beamData=True),
                       Selection(beam_end_pos.y, cnn_score, cut, beamData=True),
                       Selection(beam_end_pos.z, cnn_score, cut, beamData=True))

true_start_pos = Vector3(Selection(true_start_pos.x, cnn_score, cut),
                         Selection(true_start_pos.y, cnn_score, cut),
                         Selection(true_start_pos.z, cnn_score, cut))

true_end_pos = Vector3(Selection(true_end_pos.x, cnn_score, cut),
                       Selection(true_end_pos.y, cnn_score, cut),
                       Selection(true_end_pos.z, cnn_score, cut))

print("calculating quantities...")

print("energy residual")
energyResidual = Unwrap( (energy - mc_energy)/ mc_energy )

print("angle between beam and daughters")
beam_angle = Unwrap(BeamTrackShowerAngle(beam_start_pos.x, beam_start_pos.y, beam_start_pos.z,
                                  beam_end_pos.x, beam_end_pos.y, beam_end_pos.z,
                                  direction.x, direction.y, direction.z))

print("angle between daughters and mc particle")
mc_angle = Unwrap(DaughterRecoMCAngle(true_start_pos.x, true_start_pos.y, true_start_pos.z,
                               true_end_pos.x, true_end_pos.y, true_end_pos.z,
                               direction.x, direction.y, direction.z))

print("separation")
pair_separation = Unwrap(ShowerPairSeparation(start_pos.x, start_pos.y, start_pos.z))

print("pair angle")
pair_angle = Unwrap(ShowerPairAngle(start_pos.x, start_pos.y, start_pos.z, direction.x, direction.y, direction.z))

print("pair energy")
pair_energies, pair_leading, pair_second = ShowerPairEnergy(start_pos.x, start_pos.y, start_pos.z, energy)

pair_energies = Unwrap(pair_energies) / 1000
pair_leading = Unwrap(pair_leading) / 1000
pair_second = Unwrap(pair_second) / 1000

print("Invariant mass")
inv_mass = Unwrap(InvariantMass(start_pos.x, start_pos.y, start_pos.z, direction.x, direction.y, direction.z, energy))



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
PlotHist(inv_mass[inv_mass < 0.5], 100, "Invaraint mass of the shower pair parent (GeV)", "", 3)
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

plt.figure(2, (10, 5))

plt.subplot(122)
_, edges = PlotHist(pair_second, xlabel="Shower with the smaller energy in a pair (GeV)", title="(b)")

plt.subplot(121)
PlotHist(pair_leading, edges, xlabel="Shower with the most energy in a pair (GeV)", title="(a)")

plt.tight_layout()
"""

