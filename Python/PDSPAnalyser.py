#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 09:49:07 2021

@author: sb16165
"""

import uproot
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def PlotHist(data, bins=100, xlabel="", title="", sf=2):
    """
    plot hostogram of data and axes including bin width
    """
    height, edges, _ = plt.hist(data, bins)
    binWidth = round((edges[-1] - edges[0]) / len(edges), sf)
    plt.ylabel("Number of events (bin width=" + str(binWidth) + ")")
    plt.xlabel(xlabel)
    plt.title(title)


def PlotHist2D(data_x, data_y, bins=100, x_range=[], y_range=[], xlabel="", ylabel="", title=""):
    """
    Plots two datasets in a 2D histogram.
    """
    if len(x_range) == 2:
        data_y = data_y[data_x > x_range[0]]
        data_x = data_x[data_x > x_range[0]]
        
        data_y = data_y[data_x < x_range[1]]
        data_x = data_x[data_x < x_range[1]]

    if len(y_range) == 2:
        data_x = data_x[data_y > y_range[0]]
        data_y = data_y[data_y > y_range[0]]
        
        data_x = data_x[data_y < y_range[1]]
        data_y = data_y[data_y < y_range[1]]
        
    plt.hist2d(data_x, data_y, 100, norm=matplotlib.colors.LogNorm())
    plt.colorbar()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)


def GetData(name, unwrap=False):
    """
    Get dataset from ROOT file. will unwrap if needed
    ----- Parameters -----
    data    : recovered data from ROOT file, converted into a nested numoy array, 
              first depth is the events if the data set is nested. The following
              depths depend on the data in question
    ----------------------
    """
    data = tree_pduneana[name].arrays(library="np")[name]
    
    # consider removing
    if unwrap is True:
        return Unwrap(data)
    else:
        return data


def Unwrap(data):
    """
    Unwraps nested numpy arrays, needed because the datasets per event do not
    have the same length so is hard for plotting functions to interpret.
    Use mainly for plotting data. Note this unrwaps data at a depth of one,
    for higher depths i.e. nested arrays of nested arrays, use the function
    more than once i.e. Unravel(Unravel(...
    ----- Parameters -----
    _list   : list of data unwrapped.
    obj     : object stored at the first depth
    ----------------------
    """
    _list = []
    for obj in data:
        for item in obj:
            _list.append(item)
    return np.array(_list, object)


def Angle(x_0, y_0, z_0, x_1, y_1, z_1):
    """
    Calculates the angle between two vectors (or list of vectors)
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
    Pairs daughter events based on their start vector separation
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
    return evt_pairs


def ShowerPairSeparation(start_x, start_y, start_z, actual_pairs=False):
    """
    Gets separation of paired daughter events, either all possible pairs or 
    the true pairs (closest separation per daughter)
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
        # do stuff
        evt_pairs_dist = []
        showerPairs = GetShowerPairs(start_x, start_y, start_z)
        for i in range(len(showerPairs)):
            if(len(showerPairs[i]) == 0):
                evt_pairs_dist.append(-999)
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
                    
                    evt_pairs_dist.append( np.sqrt( x**2 + y**2 + z**2 ) )
        evt_pairs_dist = np.array( evt_pairs_dist )
    return evt_pairs_dist


def ShowerPairAngle(start_x, start_y, start_z, dir_x, dir_y, dir_z):
    """
    Gets the angle between daughter pairs
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
    Gets the Energy of the daughter pairs
    ----- Parameters -----
    pairs           : daughters paired using GetSHowerPairs
    energy          : list of energy for pair 0 and pair 1
    pairs_energy    : shower energy paired by daughters for each event
    ----------------------
    """
    pairs = GetShowerPairs(start_x, start_y, start_z)

    pairs_energy = []
    for i in range(len(pairs)):
        evt_pairs = pairs[i]
        evt_energy = shower_energy[i]
    
        energy = []
        for j in range(len(evt_pairs)):
            pair = evt_pairs[j]
            energy.append( [evt_energy[pair[0]],  evt_energy[pair[1]]] )
        pairs_energy.append(energy)
    return np.array(pairs_energy, object)


def removeDuplicates(pairs_list):
    """
    used to remove identical daugher pairs and daughter pairs without identical copies
    (daugthers that are correctly paired from GetShowerPair should have two lists that are
     reveresed i.e. [0, 3] and [3, 0])
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

    beam_dist = np.sqrt(beam_dist_x**2 + beam_dist_x**2 + beam_dist_x**2)

    angles = []
    for i in range(len(reco_daughter_allShower_dirX)):
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
    return angles


def Selection(data, value, cut, greater=True):
    """
    Selection funciton. Applies a cut wrt the value on data. Either greater than
    or less than
    ----- Parameters -----
    selected_data       : data post selection
    evt                 : data per event
    new_evt             : selected data per event
    ----------------------
    """
    selected_data = []
    for i in range(len(reco_daughter_allShower_dirX)):
        evt = reco_daughter_allShower_dirX[i]
        evt_value = value[i]
        new_evt = []
        for j in range(len(evt)):
            if evt_value[j] > cut:
                new_evt.append(evt[j])
        if len(new_evt) > 0:
            selected_data.append(new_evt)
    return np.array(selected_data, object)


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
    Calculate the shower energy residual per daughter
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
    showerPair_energy = ShowerPairEnergy(start_x, start_y, start_z, shower_energy)
    
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
        inv_mass.append(evt_inv_mass)

    return np.array(inv_mass, object)


file = uproot.open("pduneana_Prod4_1GeV_2_9_21.root")

tree_pduneana = file['pduneana/beamana']


reco_daughter_em_cnn = GetData("reco_daughter_PFP_emScore")
reco_daughter_track_cnn = GetData("reco_daughter_PFP_trackScore")
reco_daughter_nHits_collection = GetData("reco_daughter_PFP_nHits_collection")


reco_daughter_allShower_startX = GetData("reco_daughter_allShower_startX")
reco_daughter_allShower_startY = GetData("reco_daughter_allShower_startY")
reco_daughter_allShower_startZ = GetData("reco_daughter_allShower_startZ")
reco_daughter_allShower_dirX = GetData("reco_daughter_allShower_dirX")
reco_daughter_allShower_dirY = GetData("reco_daughter_allShower_dirY")
reco_daughter_allShower_dirZ = GetData("reco_daughter_allShower_dirZ")


reco_beam_startX = GetData("reco_beam_startX")
reco_beam_startY = GetData("reco_beam_startY")
reco_beam_startZ = GetData("reco_beam_startZ")
reco_beam_endX = GetData("reco_beam_endX")
reco_beam_endY = GetData("reco_beam_endY")
reco_beam_endZ = GetData("reco_beam_endZ")


reco_daughter_PFP_true_byHits_startX = GetData("reco_daughter_PFP_true_byHits_startX")
reco_daughter_PFP_true_byHits_startY = GetData("reco_daughter_PFP_true_byHits_startY")
reco_daughter_PFP_true_byHits_startZ = GetData("reco_daughter_PFP_true_byHits_startZ")
reco_daughter_PFP_true_byHits_endX = GetData("reco_daughter_PFP_true_byHits_endX")
reco_daughter_PFP_true_byHits_endY = GetData("reco_daughter_PFP_true_byHits_endY")
reco_daughter_PFP_true_byHits_endZ = GetData("reco_daughter_PFP_true_byHits_endZ")


reco_daughter_PFP_true_byHits_startE = GetData("reco_daughter_PFP_true_byHits_startE") * 1000  # GeV for some reason
reco_daughter_allShower_energy = GetData("reco_daughter_allShower_energy")


"""SELECTION"""
#cnnScore = CNNScore(reco_daughter_em_cnn, reco_daughter_track_cnn)

start_x = reco_daughter_allShower_startX
start_y = reco_daughter_allShower_startY
start_z = reco_daughter_allShower_startZ

dir_x = reco_daughter_allShower_dirX
dir_y = reco_daughter_allShower_dirY
dir_z = reco_daughter_allShower_dirZ

true_start_x = reco_daughter_PFP_true_byHits_startX
true_start_y = reco_daughter_PFP_true_byHits_startY
true_start_z = reco_daughter_PFP_true_byHits_startZ
true_end_x = reco_daughter_PFP_true_byHits_endX
true_end_y = reco_daughter_PFP_true_byHits_endY
true_end_z = reco_daughter_PFP_true_byHits_endZ

shower_energy = reco_daughter_allShower_energy

nHits = Unwrap(reco_daughter_nHits_collection)


angles = BeamTrackShowerAngle(reco_beam_startX, reco_beam_startY, reco_beam_startZ, 
                              reco_beam_endX , reco_beam_endY, reco_beam_endZ, 
                              dir_x, dir_y, dir_z)

angles = Unwrap(angles)


PlotHist2D(nHits, angles, 100, [-1, 501], [-999, np.pi], "Number of collection plane hits", "Angle between daughters and beam track")

"""
height, edges, _ = plt.hist(reco_daughter_nHits_collection, 100)
binWidth = round((edges[-1] - edges[0]) / len(edges), 2)
plt.ylabel("Number of events (bin width=" + str(binWidth) + ")")
plt.xlabel("number of collection plane hits")


height, edges, _ = plt.hist(cnnScore, 100)
binWidth = round((edges[-1] - edges[0]) / len(edges), 2)
plt.ylabel("Number of events (bin width=" + str(binWidth) + ")")
plt.xlabel("CNN score of beam daughters")
"""