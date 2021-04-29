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


def PlotHist(data, bins=100, xlabel="", title="", label="", alpha=1, sf=2, density=True):
    """
    Plot histogram of data and axes including bin width.
    ----- Parameters -----
    height      : bin height
    edges       : right edge of the bins
    ----------------------
    """
    height, edges, _ = plt.hist(data, bins, label=label, alpha=alpha, density=density)
    binWidth = round((edges[-1] - edges[0]) / len(edges), sf)
    plt.ylabel("Number of events (bin width=" + str(binWidth) + ")")
    plt.xlabel(xlabel)
    plt.title(title)
    if label != "": plt.legend()
    plt.tight_layout()
    return height, edges


def PlotHist2D(data_x, data_y, bins=100, x_range=[], y_range=[], xlabel="", ylabel="", title="", label=""):
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
    plt.hist2d(data_x, data_y, 100, norm=matplotlib.colors.LogNorm(), label=label)
    plt.colorbar()

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if label != "": plt.legend()
    plt.tight_layout()


def PlotHistComparison(data_1, data_2, bins=100, xlabel="", title="", label_1="", label_2="", alpha=1, sf=2, density=True):
    """
    Plot two histograms on the same axes, plots larger set first based on the provided bin numbers.
    ----- Parameters -----
    height_1    : bin height
    height_2    : bin height
    edges       : right edge of the bins
    ----------------------
    """

    # data_1 should be bigger than data_2 so histograms aren't cut off
    if len(data_1) < len(data_2):
       tmp = data_1
       data_1 = data_2
       data_2 = tmp
       
       tmp = label_1
       label_1 = label_2
       label_2 = label_1

    height_1, edges, _ = plt.hist(data_1, bins, label=label_1, alpha=alpha, density=density)
    height_2, _, _ = plt.hist(data_2, edges, label=label_2, alpha=alpha, density=density)

    binWidth = round((edges[-1] - edges[0]) / len(edges), sf)
    plt.ylabel("Number of events (bin width=" + str(binWidth) + ")")
    plt.xlabel(xlabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    return height_1, height_2, edges