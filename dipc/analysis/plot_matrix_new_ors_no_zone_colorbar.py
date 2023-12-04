import numpy as np
import scipy
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import distance

Z1_COLOR = (234.0/255, 29.0/255, 37.0/255)
Z5_COLOR = (43.0/255, 61.0/255, 150.0/255)

def read_leg_file(fname):
    row_hom_names = []
    row_ref_loci = []
    for input_leg_file_line in open(fname, "r"):
        ref_name, ref_locus, haplotype = input_leg_file_line.strip().split(",")
        ref_locus = int(ref_locus)
        hom_name = ref_name + ("(pat)" if haplotype == "0" else "(mat)")
        row_hom_names.append(hom_name)
        row_ref_loci.append(ref_locus)
    return (row_hom_names, row_ref_loci)

def read_name_file(fname):
    row_names = []
    for name_line in open(fname, "r"):
        row_names.append(name_line.strip())
    return row_names

def get_genomic_order(homs, loci):
    hom_to_entries = {}
    chr_order = [str(i) for i in range(1, 20)]
    chr_order.append("X")
    chr_order.append("Y")
    for chr in chr_order:
        for gender in ["(mat)", "(pat)"]:
            hom = "chr" + chr + gender
            hom_to_entries[hom] = []
    for (i, (hom, locus)) in enumerate(zip(homs, loci)):
        hom_to_entries[hom].append((hom, locus, i))
    entries = []
    for chr in chr_order:
        for gender in ["(mat)", "(pat)"]:
            hom = "chr" + chr + gender
            entries.extend(sorted(hom_to_entries[hom], key=lambda x: x[1]))
    order = [x[2] for x in entries]
    return order

def filter_out_chrX(homs, loci, names):
    new_homs = []
    new_loci = []
    new_names = []
    included_list = []
    for i in range(len(homs)):
        if "chrX" in homs[i]:
            continue
        included_list.append(i)
        new_homs.append(homs[i])
        new_loci.append(loci[i])
        new_names.append(names[i])
    return (new_homs, new_loci, new_names, included_list)


cells = []
z1_cells = []
z5_cells = []
with open('../good_cells.txt', 'r') as f:
    for line in f:
        (name, sex) = line.strip().split()
        cells.append((name, sex))
        if "z1" in name:
            z1_cells.append((name, sex))
        if "z5" in name:
            z5_cells.append((name, sex))

####################
# all OR vs all OR #
####################

def or_vs_or_dists(cell):
    values = []
    leg_homs, leg_loci = read_leg_file("OR_genes_zonal_anno.noGIOverlap_50kb.bed.leg")
    names = read_name_file("OR_genes_zonal_anno.noGIOverlap_50kb.bed.name")
    leg_homs, leg_loci, names, included_list = filter_out_chrX(leg_homs, leg_loci, names)
    assert len(leg_homs) == len(names)
    (cell_name, sex) = cell
    dists = np.loadtxt(
            ("pairwise_distances/" + cell_name + "." + sex +
            ".ors_vs_ors.pairwise_distances.txt"),
            dtype=float, delimiter="\t")
    dists = dists[:,included_list][included_list,:]
    assert not any(np.isnan(dists.flatten()))
    return (dists, names, leg_homs, leg_loci)

def load_olfr_name_to_zone():
    olfr_name_to_zone = {}
    with open("../bed_files/OR_genes_zonal_anno.bed", "r") as f:
        for line in f:
            (chr, start, stop, olfr_name, _, _, zone) = line.strip().split("\t")
            olfr_name_to_zone[olfr_name] = zone
    return olfr_name_to_zone

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

def make_chrom_colors():
    chrom_colors = {}
    for (i, color) in enumerate(
            ["#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#ffffff", "#000000"]):
        chrom_colors["chr" + str(i+1)] = color
    return chrom_colors

CHROM_COLORS = make_chrom_colors()

def hom_to_color(hom):
        chr = hom.split("(")[0]
        return CHROM_COLORS[chr]

def make_plot(vmax, use_chrom_order=False):
    # histogram for each cell
    i = 1
    all_z1_cells_dists = []
    all_z5_cells_dists = []
    olfr_name_to_zone = load_olfr_name_to_zone()
    for (cells, plot_title, color, aggregate_dists) in [
            (z1_cells, "z1OMP", Z1_COLOR, all_z1_cells_dists),
            (z5_cells, "z5OMP", Z5_COLOR, all_z5_cells_dists)]:
        plt.figure(figsize=(150, 75))
        i = 1
        for cell in cells:
        #for cell in [("z5OMP.16", "male")]:
            print(cell)
            dists, names, homs, loci = or_vs_or_dists(cell)
            y = squareform(dists)
            Z = hierarchy.linkage(y, method='single', metric='euclidean', optimal_ordering=False)
            if use_chrom_order:
                order = get_genomic_order(homs, loci)
            else:
                order = hierarchy.leaves_list(Z)


            dists = (dists[order,:])[:,order]
            print(dists.shape)
            names = [names[j].split("_")[-1] for j in order]
            homs = [homs[j] for j in order]
            loci = [loci[j] for j in order]
            plt.subplot(5, 10, i)
            assert dists.shape[0] == dists.shape[1]
            plt.imshow(dists, cmap="Reds_r", vmin=0.0, vmax=vmax)

            tracks = []

            # add zone track
            zones = [olfr_name_to_zone[olfr_name] for olfr_name in names]
            #tracks.append((zones, zone_to_color))

            # add chromosome track
            tracks.append((homs, hom_to_color))

            # add all tracks
            TRACK_WIDTH = 180
            IN_BETWEEN_GAP = 10
            BOUNDARY_GAP = 10
            ax = plt.gca()
            offset = -BOUNDARY_GAP-TRACK_WIDTH
            for (track_vals, val_to_color) in tracks:
                for (j, val) in enumerate(track_vals):
                    rect = matplotlib.patches.Rectangle((j, offset), 1, TRACK_WIDTH, fill=True, facecolor=val_to_color(val))
                    ax.add_patch(rect)
                    rect = matplotlib.patches.Rectangle((offset, j), TRACK_WIDTH, 1, fill=True, facecolor=val_to_color(val))
                    ax.add_patch(rect)
                offset -= (TRACK_WIDTH + IN_BETWEEN_GAP)
            xmin, xmax = plt.gca().get_xlim()
            ymin, ymax = plt.gca().get_ylim()
            #xmin = offset
            #ymax = offset
            xmin = offset + (TRACK_WIDTH + IN_BETWEEN_GAP)
            ymax = offset + (TRACK_WIDTH + IN_BETWEEN_GAP)

            # ticks and labels
            labels = []
            ticks = []
            if use_chrom_order:
                prev_chr = None
                for (j, hom) in enumerate(homs):
                    chr = hom.split("(")[0]
                    if chr != prev_chr:
                        labels.append(chr)
                        ticks.append(j)
                    prev_chr = chr  
                #for tick in ticks:
                    #plt.plot([xmin, xmax], [tick - 0.5, tick - 0.5], color="white", linewidth=0.7, alpha=0.5)
                    #plt.plot([tick - 0.5, tick - 0.5], [ymin, ymax], color="white", linewidth=0.7, alpha=0.5)
            else:
                for (k, name) in enumerate(names):
                    if k % 50 == 0:
                        labels.append(name)
                        ticks.append(k)

                # generate file with list of olfrs, in order
                with open("olfr_orders/" + cell[0] + "_or_vs_or_olfr_order.txt", "w") as f:
                    for (name, hom, locus) in zip(names, homs, loci): # already reordered
                        f.write("\t".join([name, hom, str(locus)]) + "\n")

            plt.gca().set_yticks(ticks)
            plt.gca().set_yticklabels(labels)
            plt.gca().set_xlim(xmin, xmax)
            plt.gca().set_ylim(ymin, ymax)
            i += 1
            plt.title(cell)
        plt.tight_layout()
        fname = plot_title + "_or_vs_or_matrices" + ("chrorder" if use_chrom_order else "percellorder") + ".vmax" + str(vmax) + ".nozonecolorbar.png"
        plt.savefig(fname)

make_plot(10, True)
make_plot(10, False)
