#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 12:16:29 2021

@author: sb16165
"""

import matplotlib
import matplotlib.pyplot as plt

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