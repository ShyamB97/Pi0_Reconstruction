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


class Vector3():
    """
    Vector data structure, currently has no operations.
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __call__(self):
        return

    def __sub__(self, other):
        return Vector3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __add__(self, other):
        return Vector3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __getitem__(self, i):
        return Vector3(self.x[i], self.y[i], self.z[i])

    def Magnitude(self):
        return ( self.x**2 + self.y**2 + self.z**2 )**0.5

    def Normalise(self):
        mag = self.Magnitude()
        return Vector3(self.x / mag, self.y / mag, self.z / mag)


def GetData(file, name, convert=False):
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
    try:
        data = tree_pduneana[name].arrays(library="np")[name]
    except uproot.KeyInFileError:
        print(name + " not found in " + "file, moving on...")
        return None
    
    # converts the top level data into a list
    # needed to convert from stlVectors(c++ nested vectors) to lists
    if convert is True:
        for i in range(len(data)):
            data[i] = data[i].tolist()
    
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
        try:
            for item in obj:
                _list.append(item)
        except TypeError:
            _list.append(obj)
    return np.array(_list, object)


# compress functions for quantities likely to be retrieved together
class Data:
    """
    Will retrieve supported datasets as needed, to keep computational times down.
    Most data applies to the shower, or MC truth particle producing the shower,
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
        return GetData(self.filename, "reco_daughter_allShower_energy") / 1000 # MeV -> GeV
    def mc_energy(self):
        return GetData(self.filename, "reco_daughter_PFP_true_byHits_startE")
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
    def true_momentum(self):
        return Vector3(
            GetData(self.filename, "reco_daughter_PFP_true_byHits_pX"),
            GetData(self.filename, "reco_daughter_PFP_true_byHits_pY"),
            GetData(self.filename, "reco_daughter_PFP_true_byHits_pZ")
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


class SelectionMask:
    """
    Class which holds a mask of the data retrieved from a root file. the mask
    has values of 1 initially and you cut the mask (1 -> 0) depending on various
    selections you make on various quantities of the data. You can repeatidly cut on
    the maske and once it is finally applied to a set of data, the data will contain
    daughter and beam events which satisfy all the cuts made to the mask.
    
    This code only works on doubly nested data (nested lists can be different in size)
    """
    def __init__(self, mask=None):
        if mask is None:
            mask = []
        self.mask = mask

    def InitiliseMask(self, reference):
        """
        Create the mask, and set its shape to the reference data
        """
        for i in range(len(reference)):
            # create a list of ones with the same shape as the ith reference i.e. beam event
            evt_mask = np.ones(reference[i].shape)
            self.mask.append(evt_mask)
        self.mask = np.array(self.mask, object)
    
    def CutMask(self, parameter, cut, conditional):
        """
        goes through each element of the mask and parameter, if the nested element of the parameter
        satisfies a condition, the mask at that position is kept (remains 1), otherwise it is cut (1 -> 0).
        """
        for i in range(len(parameter)):
            evt_parameter = parameter[i]
            evt_mask = self.mask[i]
            for j in range(len(evt_parameter)):
                
                # keep events where parameter[i][j] == cut
                if conditional == Conditional.EQUAL and evt_parameter[j] != cut:
                    evt_mask[j] = 0
                
                # keep events where parameter[i][j] != cut
                if conditional == Conditional.NOT_EQUAL and evt_parameter[j] == cut:
                    evt_mask[j] = 0
                
                # keep events where parameter[i][j] < cut
                if conditional == Conditional.LESS and evt_parameter[j] >= cut:
                    evt_mask[j] = 0
                
                # keep events where parameter[i][j] > cut
                if conditional == Conditional.GREATER and evt_parameter[j] <= cut:
                    evt_mask[j] = 0

    def ApplyMask(self, data, beamData):
        """
        Remove elements of data where mask[i][j] == 0. works on beam data by
        saying if ANY daughters have a mask = 1, it is kept.
        """
        selected_data = []
        for i in range(len(self.mask)):
            evt_mask = self.mask[i]
            evt_data = data[i]
            new_evt = []
            for j in range(len(evt_mask)):
                if evt_mask[j] == 1:
                    if beamData is True:
                        # for beam data, checking for at least 1 daughter with mask value 1
                        selected_data.append(evt_data)
                        break
                    else:
                        new_evt.append(evt_data[j])
            selected_data.append( np.array(new_evt) )

        return np.array(selected_data, object)

    def ApplyMaskVector(self, data, beamData):
        """
        The ApplyMask function but to work on a Vector3 of data.
        """
        selected_data = Vector3(self.ApplyMask(data.x, beamData), 
                        self.ApplyMask(data.y, beamData), 
                        self.ApplyMask(data.z, beamData))
        return selected_data

    def Apply(self, data, beamData=False):
        """
        Function to call when actually using this class. will call ApplyMask or
        ApplyMaskVector depeding on the type.
        """
        if type(data) is Vector3:
            selected_data = self.ApplyMaskVector(data, beamData)
        else:
            selected_data = self.ApplyMask(data, beamData)
        return selected_data

    def ApplyMaskToSelf(self):
        """
        Remove all 0 elements in the mask. Done so you can cut the mask and the data in stages
        rather than at once.
        """
        self.mask = self.ApplyMask(self.mask)




