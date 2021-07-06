#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 12:09:56 2021

@author: sb16165
"""

import numpy as np
from Master import Vector3


def Angle(v_0, v_1):
    """
    Calculates the angle between two vectors (or list of vectors).
    ----- Parameters -----
    v_n         : nth Vector3
    dot         : dot product of the two vectors
    mag_n       : magnitude of vector n
    cos_angle   : cosine of the angle
    ----------------------
    """
    if type(v_0) is not Vector3 or type(v_1) is not Vector3:
        raise TypeError("Inputs must be of type Master.Vector3")

    dot = v_0.x * v_1.x + v_0.y * v_1.y + v_0.z * v_1.z

    mag_0 = v_0.Magnitude()
    mag_1 = v_1.Magnitude()

    cos_angle = dot / ( mag_0 * mag_1 )

    # account for floating point error in calculation
    if cos_angle > 1:
        cos_angle = 1
    if cos_angle < -1:
        cos_angle = -1
    
    return np.arccos(cos_angle)


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


def GetShowerPairValues(parameter, pairs):
    """
    Used to retrieve parameters of the shower pairs
    ----- Parameters -----
    paired_values       : parameters only for the showers making a pair
    evt                 : event level parameters
    evt_pairs           : shower pairs per event
    evt_paired values   : parameters of the showers forming pairs in the event
    ----------------------
    """
    paired_values = []
    for i in range(len(parameter)):
        evt = parameter[i]
        evt_pairs = pairs[i]
        
        evt_paired_values = []
        for j in range(len(evt_pairs)):
            evt_paired_values.extend(evt[evt_pairs[j]])
        paired_values.append(evt_paired_values)
    
    return np.array(paired_values, dtype=object)


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


def ShowerPairSeparation(start_pos, shower_pairs, actual_pairs=True):
    """
    Gets separation of paired daughter events, either all possible pairs or 
    the actual pairs (closest separation per daughter).
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
        for i in range(len(shower_pairs)):
            pair_dists = []
            if(len(shower_pairs[i]) == 0):
                pair_dists.append(-999)
            else:
                for pair in shower_pairs[i]:
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


def ShowerPairAngle(shower_pairs, _dir):
    """
    Gets the angle between daughter pairs.
    ----- Parameters -----
    x_0, ...        : direction vector components of the first shower pair
    x_1, ...        : direction vector components of the second shower pair
    evt_angles      : angles of shower per event
    ----------------------
    """

    # get the angle between the shower pairs
    angles = []
    for i in range(len(shower_pairs)):
        evt_angles = []
        if(len(shower_pairs[i]) != 0):
            for pair in shower_pairs[i]:
                v_0 = _dir[i][pair[0]]
                v_1 = _dir[i][pair[1]]

                evt_angles.append( Angle(v_0, v_1) )
        angles.append(np.array(evt_angles))
    return np.array(angles, object)


def ShowerPairEnergy(shower_pairs, shower_energy):
    """
    Gets the Energy of the daughter pairs.
    ----- Parameters -----
    energy          : list of energy for pair 0 and pair 1
    pairs_energy    : shower energy paired by daughters for each event
    leading_energy  : shower in pair with the higher energy
    secondary_energy: shower in pair with the lower energy
    energy_min      : shower in pair with the higher energy
    energy_min      : shower in pair with the lower energy
    ----------------------
    """

    pairs_energy = []
    leading_energy = []
    secondary_energy = []
    for i in range(len(shower_pairs)):
        evt_pairs = shower_pairs[i]
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
        evt = daughter_dir[i]

        angle = []
        if beam_dist[i] != 0:
            for j in range(len(evt.x)):
                if evt.x[j] == -999:
                    angle.append(-999)
                else:
                    angle.append( Angle(evt[j],
                                  Vector3(beam_dist_x[i], beam_dist_y[i], beam_dist_z[i]) ) )
        else:
            # add nulls for each daughter in events with no beam
            num_daughters = len( evt.x )
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
                angle.append( Angle(Vector3(_dir.x[i][j], _dir.y[i][j], _dir.z[i][j]), 
                                    Vector3(true_dist_x[i][j], true_dist_y[i][j], true_dist_z[i][j]) ) )
        angles.append(angle)
    return np.array(angles, object)


def CNNScore(em, track):
    """
    Calculate the em/track like CNN score per daughter
    ----- Parameters -----
    cnnScore    : normalised em/shower like score
    ----------------------
    """
    cnnScore = []
    for i in range(len(em)):
        cnnScore_evt = []
        for j in range(len(em[i])):
            if em[i][j] != -999:
                cnnScore_evt.append( em[i][j] / (em[i][j] + track[i][j]) )
            else:
                cnnScore_evt.append(-999)
        cnnScore.append(cnnScore_evt)

    return np.array(cnnScore, object)


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


def InvariantMass(shower_pair_angles, shower_pair_energy):
    """
    Calculate the invariant mass from shower pairs.
    ----- Parameters -----
    inv_mass                : invariant mass per event per shower pair
    evt_angle               : angle between shower pairs for an event
    showerPair_energy       : energy of the showers forming the pair for an event
    angle                   : opening angle for a pair
    energy                  : shower pair energies
    ----------------------
    """    
    inv_mass = []
    for i in range(len(shower_pair_angles)):
        evt_angle = shower_pair_angles[i]
        evt_energy = shower_pair_energy[i]

        evt_inv_mass = []
        for j in range(len(evt_angle)):
            angle = evt_angle[j]
            energies = evt_energy[j]

            if angle != -999:
                if energies[0] >= 0 and energies[1] >= 0:
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


def GetShowerStartHits(hit_radial, hit_longitudinal, geometry=[1, -1, 4]):
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
    return np.array(start_hits, object)

