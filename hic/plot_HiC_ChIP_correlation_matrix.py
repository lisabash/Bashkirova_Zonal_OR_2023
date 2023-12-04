import numpy as np
import scipy
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import distance
import pandas as pd

# used to plot heatmap in Supp Figure 2D
# input is matrix of Hi-C contacts ordered by level of histone ChIP expression (generated with R script)

# to run:
# (make sure you are in base environment)
# python plot_matrix.py

#zonal_order_file = "baits.OEtk79_order.within_zonal_category.tsv"
#zonal_order_file = "baits.OEtk79_order.zones_intermixed.tsv"
#zonal_order_file = "baits.Z1OMPtk79_order.within_zonal_category.tsv"
#zonal_order_file = "baits.Z1OMPtk79_order.zones_intermixed.tsv"
#zonal_order_file = "baits.Z2OMPtk79_order.within_zonal_category.tsv"
#zonal_order_file = "baits.Z2OMPtk79_order.zones_intermixed.tsv"
#zonal_order_file = "baits.Z5OMPtk79_order.within_zonal_category.tsv"
#zonal_order_file = "baits.Z5OMPtk79_order.zones_intermixed.tsv"
#zonal_order_file = "baits.Z1OMPtk79_order.within_zonal_category.tsv"
#zonal_order_file = "baits.Z1Ngntk79_order.within_zonal_category.tsv"
#zonal_order_file = "baits.Z5Ngntk79_order.within_zonal_category.tsv"
#zonal_order_file = "baits.Z1Ngntk79_order.zones_intermixed.tsv"
zonal_order_file = "baits.Z5Ngntk79_order.zones_intermixed.tsv"
#matrix_file = "Z1.Ngn.Zonal_OR_to_Zonal_OR_inter.mat.csv"
matrix_file = "Z5.Ngn.Zonal_OR_to_Zonal_OR_inter.mat.csv"
#matrix_file = "OMP.Zonal_OR_to_Zonal_OR_inter.mat.csv"
#matrix_file = "P2.Zonal_OR_to_Zonal_OR_inter.mat.csv"
#matrix_file = "Z1.OMP.Zonal_OR_to_Zonal_OR_inter.mat.csv"
#matrix_file = "Z5.OMP.Zonal_OR_to_Zonal_OR_inter.mat.csv"


output_plot_file = "heatmap_custom_order.Z5.Ngn.Z5_Ngn_tk79_intermixed.vmax20.pdf"
tk79_value_column = "Z5_Ngn_tk79"
vmax = 20

def get_genomic_order(chrs, loci):
    chr_to_entries = {}
    chr_order = [str(i) for i in range(1, 20)]
    chr_order.append("X")
    chr_order.append("Y")
    for chr in chr_order:
        chr_to_entries[chr] = []
    for (i, (chr, locus)) in enumerate(zip(chrs, loci)):
        chr_to_entries[chr].append((chr, locus, i))
    entries = []
    for chr in chr_order:
        entries.extend(sorted(chr_to_entries[chr], key=lambda x: x[1]))
    order = [x[2] for x in entries]
    return order

def zone_to_color(zone):
    if zone == "1":
        return (234.0/255,29.0/255,37.0/255)
    elif (zone == "5" or zone == "4"):
        return (34.0/255,70.0/255,200.0/255)
    elif (zone == "2" or zone == "3"):
        return (99.0/255,181.0/255,51.0/255)
    elif (zone == "classI"):
        return (0.8, 0.8, 0.8)
    else:
        assert False
        return "white"


def to_chrom_and_pos(input_str):
    chr, start_str = input_str.split("_")
    start = int(float(start_str))
    return (chr, start)

def read_table():
    data = pd.read_csv(matrix_file) # test.csv
    ors = list(map(to_chrom_and_pos, data.iloc[:,0]))
    mat = np.array(data.iloc[:,1:])
    print(mat)
    assert mat.shape == (len(ors), len(ors))
    return (ors, mat)

def load_tk79val():
    df = pd.read_csv(zonal_order_file, sep="\t")
    baits = list(df.bait)
    zones = list(df.zone)
    baits = list(map(to_chrom_and_pos, baits))
    values = list(df[tk79_value_column])#Z5_OSN_tk79)
    bait_to_value = {}
    bait_to_zone = {}
    for (bait, value, zone) in zip(baits, values, zones):
        bait_to_value[bait] = value
        bait_to_zone[bait] = zone
    return baits, bait_to_value, bait_to_zone

# ordered by chrom
# custom order <--- as a list of chr_start in order
# hierarchical clustering

CHROM_ORDER = "chrom"
CUSTOM_ORDER = "custom"
CLUSTERING_ORDER = "clustering"

def make_plot(order_mode, output_fname, custom_order=None):

    (baits, matrix) = read_table()

    #bait_to_zone = get_bait_to_zone()
    #print(bait_to_zone) # TODO 
    #zones = [bait_to_zone[bait] for bait in baits]
    bait_to_original_idx = {}
    for (i, bait) in enumerate(baits):
        bait_to_original_idx[bait] = i

    plt.figure(figsize=(25, 25))

    #plt.imshow(np.log(1 + matrix), cmap="Reds")
    #plt.imshow(matrix, cmap="Reds")
    #plt.savefig(output_fname)
    #return

    tracks = []

    if order_mode == CLUSTERING_ORDER:
        assert False
        #assert matrix.shape[0] == matrix.shape[1]
        #for row in range(len(baits)):
            #print((baits[row], np.max(np.abs((matrix[row,:] - matrix[:,row])))))
        #matrix = (matrix + matrix.T)/2 # 
        #dists = np.divide(1.0, 1.0 + matrix) # TODO other ways of doing this..
        #for i in range(len(baits)):
            #dists[i,i] = 0.0
        #y = squareform(dists)
        #Z = hierarchy.linkage(y, method='single', metric='cosine', optimal_ordering=False)
        #order = hierarchy.leaves_list(Z)
    elif order_mode == CHROM_ORDER:
        chrs = [bait[0] for bait in baits]
        loci = [bait[1] for bait in baits]
        order = get_genomic_order(chrs, loci)
        # TODO add chromosome track
    elif order_mode == CUSTOM_ORDER:

        ordered_baits, bait_to_value, bait_to_zone = load_tk79val()

        #print("**************************************")
        #print(bait_to_zone[('10',79000000)])
        #print(bait_to_value[('10',79000000)])

        zones = [bait_to_zone[bait] for bait in baits]
        zz_lookup = { "zone1" : "1", "zone2to3" : "2", "zone4to5" : "4"}
        zones = list(map(lambda zz: zz_lookup[zz], zones))
        values = [bait_to_value[bait] for bait in baits]

        # order

        # get order from highest 79 value to lowest, and check it
        #order = np.flip(np.argsort(values))
        #min_so_far = np.inf
        #for j in order:
            #assert values[j] <= min_so_far
            #min_so_far = values[j]

        order = [bait_to_original_idx[bait] for bait in ordered_baits]
        values = [values[j] for j in order]
        print(values)
        max_val = np.max(values)
        tracks.append((values, lambda val: (1-(val/max_val), 1-(val/max_val), 1-(val/max_val))))

    # reorder everything
    baits = [baits[j] for j in order]
    zones = [zones[j] for j in order]
    matrix = (matrix[order,:])[:,order]
    print(matrix.shape)
    assert len(baits) == len(zones)
    assert len(baits) == matrix.shape[0]
    assert matrix.shape[0] == matrix.shape[1]

    cmap = matplotlib.cm.Reds
    cmap.set_bad('white',1.)
    plt.imshow(matrix, cmap=cmap, vmin=0.0, vmax=vmax) # TODO how to set this parameter..

    # add zone track
    #colors = map(zone_to_color, zones)
    tracks.append((zones, zone_to_color))

    # add all tracks
    TRACK_WIDTH = 20
    BOUNDARY_GAP = 10
    ax = plt.gca()
    offset = -BOUNDARY_GAP-TRACK_WIDTH
    for (track_vals, val_to_color) in tracks:
        for (j, val) in enumerate(track_vals):
            rect = matplotlib.patches.Rectangle((j, offset), 1, TRACK_WIDTH, fill=True, facecolor=val_to_color(val))
            ax.add_patch(rect)
            rect = matplotlib.patches.Rectangle((offset, j), TRACK_WIDTH, 1, fill=True, facecolor=val_to_color(val))
            ax.add_patch(rect)
        offset -= TRACK_WIDTH
    xmin, xmax = plt.gca().get_xlim()
    ymin, ymax = plt.gca().get_ylim()
    xmin = offset
    ymax = offset

    # color ORs by zone

    # if chrom order, then show lines for chromosome
    # TODO remove
    #if order_mode == CHROM_ORDER:
        #chr_idx = []
        #prev_chr = None
        #for (j, bait) in enumerate(baits):
            #chr = bait[0]
            #if chr != prev_chr:
                ##labels.append(chr)
                ##ticks.append(j)
                #chr_idx.append(j) 
            #prev_chr = chr  
        #for j in chr_idx:
            #plt.plot([xmin, xmax], [j - 0.5, j - 0.5], color="black")#, linewidth=0.7, alpha=0.5)
            #plt.plot([j - 0.5, j - 0.5], [ymin, ymax], color="black")#, linewidth=0.7, alpha=0.5)

    # tick labels
    labels = []
    ticks = []
    for (k, bait) in enumerate(baits):
        labels.append(bait[0] + "_" + str(bait[1]))
        ticks.append(k)
    plt.gca().set_yticks(ticks)
    plt.gca().set_yticklabels(labels)
    plt.gca().set_xlim(xmin, xmax)
    plt.gca().set_ylim(ymin, ymax)

    #plt.title(cell)
    plt.tight_layout()
    plt.savefig(output_fname)

#make_plot(CHROM_ORDER, "heatmap_chr_order.pdf")
#make_plot(CLUSTERING_ORDER, "heatmap_clustering_order.pdf")
make_plot(CUSTOM_ORDER, output_plot_file) # TODO put the custom order here
