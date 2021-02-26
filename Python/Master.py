#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 13:50:53 2021

@author: sb16165
"""

import uproot
import numpy as np


class Vector3:
    """
    Vector data structure, currently has no operations.
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


def GetData(file, name, unwrap=False):
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
        return GetData(self.filename, "reco_daughter_PFP_emScore")
    def cnn_track(self):
        return GetData(self.filename, "reco_daughter_PFP_trackScore")

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
