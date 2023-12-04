import sys
import subprocess

ref_con_fname = sys.argv[1]
cells_list_fname = sys.argv[2]
mode = sys.argv[3]

def compute_density_ratio_for_cell(cell_name, sex):

    # NOTE: the results change a little bit when you use impute versus contacts
    if mode == "impute":
        con_fname = "../merged." + cell_name + ".ca." + sex + ".impute.con.gz"
    elif mode == "contacts":
        con_fname = "../merged." + cell_name + ".ca." + sex + ".contacts.con.gz"
    else:
        raise ValueError("mode must be impute or contacts")

    # total counts for 10Mb x 10Mb regions
    output = subprocess.check_output(
            ['../dip-c-master/dip-c', 'ard', '-d10000000', '-c', ref_con_fname, '-n', con_fname])
    big_total = float(sum(map(int, output.split())))
    
    # total counts for 100kb x 100kb regions
    output = subprocess.check_output(
            ['../dip-c-master/dip-c', 'ard', '-d100000', '-c', ref_con_fname, '-n', con_fname])
    little_total = float(sum(map(int, output.split())))
    
    numerator = little_total / (100000 * 100000)
    denominator = (big_total - little_total) / ((10000000 * 10000000) - (100000 * 100000))
    ratio = numerator / denominator
    return ratio

with open(cells_list_fname, "r") as f:
    for line in f:
        (cell_name, sex) = line.strip().split("\t")
        ratio = compute_density_ratio_for_cell(cell_name, sex)
        print(cell_name + "\t" + sex + "\t" + str(ratio))
    
