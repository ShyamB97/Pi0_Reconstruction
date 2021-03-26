#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 13:50:53 2021

@author: sb16165
"""

import uproot
import numpy as np
from enum import Enum

class Conditional(Enum):
    GREATER = 1
    LESS = 2
    EQUAL = 3
    NOT_EQUAL = 4


class Vector3:
    """
    Vector data structure, currently has no operations.
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


def GetData(file, name, convert=False, unwrap=False):
    """
    Get dataset from ROOT file. will unwrap if needed.
    ----- Parameters -----
    tree_pduneana   : directory the datasets are stored in, by default this should
                      never change.
    data            : recovered data from ROOT file, converted into a nested numoy array, 
                      first depth is the events if the data set is nested. The following
                      depths depend on the data in question
    ----------------------
    """
    file = uproot.open(file)
    tree_pduneana = file['pduneana/beamana']
    data = tree_pduneana[name].arrays(library="np")[name]
    
    # converts the top level data into a list
    # needed to convert from stlVectors(c++ nested vectors) to lists
    if convert is True:
        for i in range(len(data)):
            data[i] = data[i].tolist()
    
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


# compress functions for quantities likely to be retrieved together
class Data:
    """
    Will retrieve supported datasets as needed, to keep computational times down.
    any most data applies to the shower, or MC truth particle producing the shower,
    if not it is stated in the name.
    """
    def __init__(self, filename):
        self.filename = filename
    
    # reco daughter data
    def start_pos(self):
        return Vector3(
            GetData(self.filename, "reco_daughter_allShower_startX"),
            GetData(self.filename, "reco_daughter_allShower_startY"),
            GetData(self.filename, "reco_daughter_allShower_startZ")
            )
    def direction(self):
        return Vector3(
            GetData(self.filename, "reco_daughter_allShower_dirX"),
            GetData(self.filename, "reco_daughter_allShower_dirY"),
            GetData(self.filename, "reco_daughter_allShower_dirZ")
            )
    def energy(self):
        return GetData(self.filename, "reco_daughter_allShower_energy")
    def mc_energy(self):
        return GetData(self.filename, "reco_daughter_PFP_true_byHits_startE") * 1000  # GeV for some reason
    def nHits(self):
        return GetData(self.filename, "reco_daughter_PFP_nHits_collection")
    def cnn_em(self):
        return GetData(self.filename, "reco_daughter_PFP_emScore_collection")
    def cnn_track(self):
        return GetData(self.filename, "reco_daughter_PFP_trackScore_collection")

    # mc truth daughter data
    def true_start_pos(self):
        return Vector3(
            GetData(self.filename, "reco_daughter_PFP_true_byHits_startX"),
            GetData(self.filename, "reco_daughter_PFP_true_byHits_startY"),
            GetData(self.filename, "reco_daughter_PFP_true_byHits_startZ")
            )
    def true_end_pos(self):
        return Vector3(
            GetData(self.filename, "reco_daughter_PFP_true_byHits_endX"),
            GetData(self.filename, "reco_daughter_PFP_true_byHits_endY"),
            GetData(self.filename, "reco_daughter_PFP_true_byHits_endZ")
            )

    # beam particle (parent) data
    def beam_start_pos(self):
        return Vector3(
            GetData(self.filename, "reco_beam_startX"),
            GetData(self.filename, "reco_beam_startY"),
            GetData(self.filename, "reco_beam_startZ")
            )
    def beam_end_pos(self):
        return Vector3(
            GetData(self.filename, "reco_beam_endX"),
            GetData(self.filename, "reco_beam_endY"),
            GetData(self.filename, "reco_beam_endZ")
            )
    
    def pandoraTag(self):
        return GetData(self.filename, "pandoraTag")
    
    # quantities to calculate cylinder hits
    def hit_radial(self):
        return GetData(self.filename, "hitRadial", convert=True)
    def hit_longitudinal(self):
        return GetData(self.filename, "hitLongitudinal", convert=True)
    
    # CNN score without averaging the track and EM score (i.e. not the same way as done in PDSPAnalyser.py)
    def CNNScore(self):
        return GetData(self.filename, "CNNScore")


### OLD SELECTION CLASS ###
class Selection:
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
    def __init__(self, parameter, cut, greater=True):
        self.parameter = parameter
        self.cut = cut
        self.greater = greater


    def Cut(self, data, beamData):

        selected_data = []
        for i in range(len(self.parameter)):
            evt = data[i]
            evt_value = self.parameter[i]
            new_evt = []
            for j in range(len(evt_value)):
    
                if self.greater is True and evt_value[j] > self.cut:
                    if beamData is True:
                        selected_data.append(evt)
                        break
                    else:
                        new_evt.append(evt[j])
    
                if self.greater is False and evt_value[j] < self.cut:
                    if beamData is True:
                        selected_data.append(evt)
                        break
                    else:
                        new_evt.append(evt[j])
    
            new_evt = np.array(new_evt)
            if len(new_evt) > 0:
                selected_data.append(new_evt)
        return np.array(selected_data, object)


    def CutVector(self, data, beamData):

        selected_data = Vector3(self.Cut(data.x, beamData), 
                                self.Cut(data.y, beamData), 
                                self.Cut(data.z, beamData))
        return selected_data


    def Apply(self, data, beamData=False):
        
        if type(data) is Vector3:
            selected_data = self.CutVector(data, beamData)
        else:
            selected_data = self.Cut(data, beamData)
        return selected_data


class SelectionMask:
    def __init__(self, mask=None):
        if mask is None:
            mask = []
        self.mask = mask

    def InitiliseMask(self, reference):
        for i in range(len(reference)):
            evt_mask = np.ones(reference[i].shape)
            self.mask.append(evt_mask)        
        self.mask = np.array(self.mask, object)
    
    def CutMask(self, parameter, cut, conditional):
        for i in range(len(parameter)):
            evt_parameter = parameter[i]
            evt_mask = self.mask[i]
            for j in range(len(evt_parameter)):
                
                if conditional == Conditional.EQUAL and evt_parameter[j] != cut:
                    evt_mask[j] = 0
                    
                if conditional == Conditional.NOT_EQUAL and evt_parameter[j] == cut:
                    evt_mask[j] = 0
                
                if conditional == Conditional.LESS and evt_parameter[j] > cut:
                    evt_mask[j] = 0
                
                if conditional == Conditional.GREATER and evt_parameter[j] < cut:
                    evt_mask[j] = 0

    def ApplyMask(self, data, beamData):
        selected_data = []
        for i in range(len(self.mask)):
            evt_mask = self.mask[i]
            evt_data = data[i]
            new_evt = []
            for j in range(len(evt_mask)):
                if evt_mask[j] == 1:
                    if beamData is True:
                        selected_data.append(evt_data)
                        break
                    else:
                        new_evt.append(evt_data[j])
            selected_data.append( np.array(new_evt) )

        return np.array(selected_data, object)

    def ApplyMaskVector(self, data, beamData):
        selected_data = Vector3(self.ApplyMask(data.x, beamData), 
                        self.ApplyMask(data.y, beamData), 
                        self.ApplyMask(data.z, beamData))
        return selected_data

    def Apply(self, data, beamData=False):
        if type(data) is Vector3:
            selected_data = self.ApplyMaskVector(data, beamData)
        else:
            selected_data = self.ApplyMask(data, beamData)
        return selected_data

    def ApplyMaskToSelf(self):
        self.mask = self.ApplyMask(self.mask)




