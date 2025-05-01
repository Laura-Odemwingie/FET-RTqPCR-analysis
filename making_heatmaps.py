import csv
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 

import matplotlib
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# -------------- Heatmap code from Matplotlib : IGNORE ---------------
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
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.02)
    cbar = ax.figure.colorbar(im, cax=cax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "black"),
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

# -----------------------------------------------------------------------

#------------------ Merging of all agerage replicates ----------------

# ADD BELOW THE FULL PATH TO THE CSV FILE DIRECTORY - TELL IT WHERE TO FIND YOUR FILES
#   - Makes it easier to access the csv data files

csvFile_path = "~/Desktop/python_qPCR_analysis/RT_qPCR Data for Heatmap.csv"
OutputDir = "~/Desktop/python_qPCR_analysis/Heatmaps/"

#----------------------------------------------------------------------
#  DEFINE DATA STRUCTURE OR READING IN DATA?
# Define a class (data structure) for each age-region-gene expression ddCT values
# Telling my coputer that I will be creating a data set in this structure. It's like a bluprint.???????
# you will be ginen values to be stored in the classs/blueprint I have created 
class data:
  def __init__(self, age:str, region:str, gene:str, Avg_ddct_norm1:float, Avg_ddct_norm2:float):
    self.age = age
    self.region = region
    self.gene = gene
    self.Avg_ddct_norm1 = round(Avg_ddct_norm1, 3)
    self.Avg_ddct_norm2 = round(Avg_ddct_norm2, 3)


# Define global vmin and vmax values to use for range of color bar
VMIN = -2.0
VMAX = 2.0


# Open a file
dataFile = open(csvFile_path)
# Read all lines into a list called MyLines
MyLines = dataFile.readlines()
# close file
dataFile.close()

# Remove un-used header row
MyLines.pop(0)
# Go over all lines and create a `data` class instance for each. Store into a list called `rowData`
rowData = []
for line in MyLines:
    # remove newline character
    line = line.rstrip('\n')
    # Every time you see a comma, split into new column
    line = line.split(',')
    # Test to see if avg norm1 data exists
    AvgDdctNorm1 = 0.0
    if(line[3] != ''):
       AvgDdctNorm1 = float(line[3])
    # create a new `data` object for this row
    tempData = data(age=str(line[0]), region=str(line[1]), gene=str(line[2]), Avg_ddct_norm1=AvgDdctNorm1, Avg_ddct_norm2=float(line[4]))
    # store the row
    rowData.append(tempData)

Ages = []
Regions = []
Regions_noSC = []
Genes = []
for row in rowData:
    # Extract all ages into a list from rowData
    if row.age not in Ages:
      Ages.append(row.age)
    # Extract all regions into a list from rowData
    if row.region not in Regions:
      Regions.append(row.region)
      if row.region != "Spinal cord":
         Regions_noSC.append(row.region)
    # Extract all Genes into a list from rowData
    if row.gene not in Genes:
      Genes.append(row.gene)



# ------ Heatmap 1.1 : Gene x Region_noSC (ddct Norm1) 3 months ------

# Initialise heatmap data structure to correct shape
ddctNorm1_3m_gXr = np.zeros(shape=(len(Genes), len(Regions_noSC)))

for gene in range(0, len(Genes)):
   for reg in range(0, len(Regions_noSC)):
      # Find the row which corresponds to the current gene and region (excluding spinal cord)
      for row in rowData:
         if row.gene == Genes[gene] and row.region == Regions_noSC[reg] and row.age == "3m":
            # When we find the row, add it to the heatmap data structure and move onto next gene-region pair
            ddctNorm1_3m_gXr[gene][reg] = row.Avg_ddct_norm1
            break

# Create heatmap
fig, ax = plt.subplots()
# Uses matplotlib functions defined at top of file
im, cbar = heatmap(ddctNorm1_3m_gXr, Genes, Regions_noSC, ax=ax,
                   cmap="RdBu_r", cbarlabel="∆∆CT", vmin=VMIN, vmax=VMAX)

texts = annotate_heatmap(im, valfmt="{x:.3f}")
fig.tight_layout()
ax.set_title("FET expression in 3 months Ntg mouse brain", fontsize=16)
plt.savefig(OutputDir + row.age + "_" + "regions" + "_" + "Heatmap.png")

#cmap="YlGn"
# ------ Heatmap 1.2 : Gene x Region_noSC (ddct Norm1) 12 months ------

# Initialise heatmap data structure to correct shape
ddctNorm1_12m_gXr = np.zeros(shape=(len(Genes), len(Regions_noSC)))

for gene in range(0, len(Genes)):
   for reg in range(0, len(Regions_noSC)):
      # Find the row which corresponds to the current gene and region (excluding spinal cord)
      for row in rowData:
         if row.gene == Genes[gene] and row.region == Regions_noSC[reg] and row.age == "12m":
            # When we find the row, add it to the heatmap data structure and move onto next gene-region pair
            ddctNorm1_12m_gXr[gene][reg] = row.Avg_ddct_norm1
            break

# Create heatmap
fig, ax = plt.subplots()
# Uses matplotlib functions defined at top of file
im, cbar = heatmap(ddctNorm1_12m_gXr, Genes, Regions_noSC, ax=ax,
                   cmap="RdBu_r", cbarlabel="∆∆CT", vmin=VMIN, vmax=VMAX)
texts = annotate_heatmap(im, valfmt="{x:.3f}")
fig.tight_layout()
ax.set_title("FET expression in 12 months Ntg mouse brain", fontsize=16)
plt.savefig(OutputDir + row.age + "_" + "regions" + "_" + "Heatmap.png")

# ------ Heatmap 1.3 : Gene x Region_noSC (ddct Norm1) 24 months ------

# Initialise heatmap data structure to correct shape
ddctNorm1_24m_gXr = np.zeros(shape=(len(Genes), len(Regions_noSC)))

for gene in range(0, len(Genes)):
   for reg in range(0, len(Regions_noSC)):
      # Find the row which corresponds to the current gene and region (excluding spinal cord)
      for row in rowData:
         if row.gene == Genes[gene] and row.region == Regions_noSC[reg] and row.age == "24m":
            # When we find the row, add it to the heatmap data structure and move onto next gene-region pair
            ddctNorm1_12m_gXr[gene][reg] = row.Avg_ddct_norm1
            break

# Create heatmap
fig, ax = plt.subplots()
# Uses matplotlib functions defined at top of file
im, cbar = heatmap(ddctNorm1_12m_gXr, Genes, Regions_noSC, ax=ax,
                   cmap="RdBu_r", cbarlabel="∆∆CT", vmin=VMIN, vmax=VMAX)
texts = annotate_heatmap(im, valfmt="{x:.3f}")
fig.tight_layout()
ax.set_title("FET expression in 24 months Ntg mouse brain", fontsize=16)
plt.savefig(OutputDir + row.age + "_" + "regions" + "_" + "Heatmap.png")

####-------------------------------------------------------------------
####-------------------------------------------------------------------

# ------ Heatmap 2.1 : Age x Region_noSC (ddct Norm2) FUS ------

# Initialise heatmap data structure to correct shape
ddctNorm2_FUS_gXr = np.zeros(shape=(len(Ages), len(Regions)))

for age in range(0, len(Ages)):
   for reg in range(0, len(Regions)):
      # Find the row which corresponds to the current gene and region (excluding spinal cord)
      for row in rowData:
         if row.age == Ages[age] and row.region == Regions[reg] and row.gene == "FUS":
            # When we find the row, add it to the heatmap data structure and move onto next gene-region pair
            ddctNorm2_FUS_gXr[age][reg] = row.Avg_ddct_norm2
            break

# Create heatmap
fig, ax = plt.subplots()
# Uses matplotlib functions defined at top of file
im, cbar = heatmap(ddctNorm2_FUS_gXr, Ages, Regions, ax=ax,
                   cmap="RdBu_r", cbarlabel="∆∆CT", vmin=VMIN, vmax=VMAX)
texts = annotate_heatmap(im, valfmt="{x:.3f}")
fig.tight_layout()
ax.set_title("FUS expression in aging Ntg mouse brain", fontsize=16)
plt.savefig(OutputDir + row.gene + "_"+ "Ageing" + "_" + "Heatmap.png")

# ------ Heatmap 2.2 : Age x Region_noSC (ddct Norm2) EWS ------

# Initialise heatmap data structure to correct shape
ddctNorm2_EWS_gXr = np.zeros(shape=(len(Ages), len(Regions)))

for age in range(0, len(Ages)):
   for reg in range(0, len(Regions)):
      # Find the row which corresponds to the current gene and region (excluding spinal cord)
      for row in rowData:
         if row.age == Ages[age] and row.region == Regions[reg] and row.gene == "EWS":
            # When we find the row, add it to the heatmap data structure and move onto next gene-region pair
            ddctNorm2_EWS_gXr[age][reg] = row.Avg_ddct_norm2
            break

# Create heatmap
fig, ax = plt.subplots()
# Uses matplotlib functions defined at top of file
im, cbar = heatmap(ddctNorm2_EWS_gXr, Ages, Regions, ax=ax,
                   cmap="RdBu_r", cbarlabel="∆∆CT", vmin=VMIN, vmax=VMAX)
texts = annotate_heatmap(im, valfmt="{x:.3f}")
fig.tight_layout()
ax.set_title("EWS expression in aging Ntg mouse brain", fontsize=16)
plt.savefig(OutputDir + row.gene + "_"+ "Ageing" + "_" + "Heatmap.png")

# ------ Heatmap 2.3 : Age x Region_noSC (ddct Norm2) TAF15 ------

# Initialise heatmap data structure to correct shape
ddctNorm2_TAF15_gXr = np.zeros(shape=(len(Ages), len(Regions)))

for age in range(0, len(Ages)):
   for reg in range(0, len(Regions)):
      # Find the row which corresponds to the current gene and region (excluding spinal cord)
      for row in rowData:
         if row.age == Ages[age] and row.region == Regions[reg] and row.gene == "TAF15":
            # When we find the row, add it to the heatmap data structure and move onto next gene-region pair
            ddctNorm2_TAF15_gXr[age][reg] = row.Avg_ddct_norm2
            break

# Create heatmap
fig, ax = plt.subplots()
# Uses matplotlib functions defined at top of file
im, cbar = heatmap(ddctNorm2_TAF15_gXr, Ages, Regions, ax=ax,
                   cmap="RdBu_r", cbarlabel="∆∆CT", vmin=VMIN, vmax=VMAX)
texts = annotate_heatmap(im, valfmt="{x:.3f}")
fig.tight_layout()
ax.set_title("TAF15 expression in aging Ntg mouse brain", fontsize=16)
plt.savefig(OutputDir + row.gene + "_"+ "Ageing" + "_" + "Heatmap.png")