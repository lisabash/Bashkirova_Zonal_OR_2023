import sys

bed_fname1 = sys.argv[1]
bed_fname2 = sys.argv[2]

def load_bed(bed_fname):
    data = []
    with open(bed_fname, "r") as f:
        for line in f:
            tokens = line.strip().split()
            chrom, start, stop = tokens[0:3]
            data.append((chrom, int(start), int(stop)))
    return data

def middle(start, stop):
    return int((start + stop)/2)

data1 = load_bed(bed_fname1)
data2 = load_bed(bed_fname2)

for (chrom1, start1, stop1) in data1:
    mid1 = middle(start1, stop1)
    for (chrom2, start2, stop2) in data2:
        mid2 = middle(start2, stop2)
        print(chrom1 + "," + str(mid1) + "," + "." + "\t" + chrom2 + "," + str(mid2) + "," + ".")
