#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import sys
import glob
import os

## Heatmap scripts adapted from: https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html#sphx-glr-gallery-images-contours-and-fields-image-annotated-heatmap-py

def Matrix_Elements(input_list):
    """Makes list of lists of matrix elements"""
    Out = []
    for entry1 in input_list:
        Line = []
        for entry2 in input_list:
            Data = [entry1, entry2]
            Line.append(Data)
        Out.append(Line)
    return Out

def Matrix_Matcher(input_list, element_matches, elements):
    """Makes a matrix based on the linked data from an element list consisting of matches and data"""
    Elements = Matrix_Elements(input_list)
    Out = []
    for line in Elements:
        New = []
        for entry in line:
            if entry[0] == entry[1]:
                New.append(1)
            elif (entry in element_matches):
                Pos = element_matches.index(entry)
                New.append(elements[Pos])
            elif ([entry[1], entry[0]] in element_matches):
                Pos = element_matches.index([entry[1], entry[0]])
                New.append(elements[Pos])
            else:
                New.append(0)
        Out.append(New)
    return Out

def PLASMAR_Overlap_Matrix(input_overlap_file):
    """Makes a matrix of the overlap probabilities from an input overlap file"""
    Matches = []
    Probs = []
    f = open(input_overlap_file, 'r')
    String1 = f.readline()
    for line in f:
        List1 = line.split()
        Matches.append([List1[0], List1[1]])
        Probs.append(float(List1[-1]))
    f.close()
    IDs = []
    for entry in Matches:
        if (entry[0] in IDs) == False:
            IDs.append(entry[0])
        if (entry[1] in IDs) == False:
            IDs.append(entry[1])
    IDs.sort()
    New = Matrix_Matcher(IDs, Matches, Probs)
    return [IDs, New]

def Array_Maker(input_list):
    Out = np.array(input_list)
    return Out

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-45, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

def func(x, pos):
    return f"{x:.2f}".replace("0.00", "0").replace("1.00", "1")

def PLASMAR_Heatmap_Maker(input_overlap_file, output_file):
    Info = PLASMAR_Overlap_Matrix(input_overlap_file)
    Data = Array_Maker(Info[1])
    fig, ax = plt.subplots()
    im, cbar = heatmap(Data, Info[0], Info[0], ax=ax, cmap="YlGn", vmin=0, vmax=1, cbarlabel="Fraction of Models Predicting Overlap")
    texts = annotate_heatmap(im, valfmt=matplotlib.ticker.FuncFormatter(func), size=6)
    fig.tight_layout()
    plt.savefig(output_file)

List1 = glob.glob('./Reports/Overlap_Comp_Report_*.txt')
for files in List1:
    Allele = files.split('_')[-1][0:-4]
    PLASMAR_Heatmap_Maker(files, './Reports/' + Allele + '_Heatmap.png')
