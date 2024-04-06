import math
import matplotlib.pyplot as plt
import numpy as np
import os
from graphviz import Graph
from math import cos, pi, sin
from operator import itemgetter

inp = "../../Conferences and Papers/2023 CIBCB/WSDA/BitOutDone/"
outp = "./BitProcessed/ForJames/"
finame = "best.lint"
samps = 30
precision = 6
col_width = 6 + precision


def get_data(dir_path: str):
    fits = []
    networks = []
    SDAs = []
    edges = []
    weights = []
    hists = []
    fi_str = finame
    cnt = 0
    with open(dir_path + fi_str) as f:
        lines = f.readlines()
        next_graph = False
        cont_graph = False
        next_SDA = False
        cont_SDA = False
        network = []
        SDA = []
        for line in lines:
            if line.__contains__(str("-fitness")):
                cnt += 1
                d = line.split(" ")
                fits.append(float(d[0]))
                if cont_graph:
                    networks.append(network)
                    network = []
                cont_graph = False
                pass
            elif line.__contains__(str("Self-Driving Automata")):
                next_SDA = True
                pass
            elif line.__contains__("Graph"):
                cont_SDA = False
                SDAs.append(SDA)
                SDA = []
                pass
            elif line.__contains__("Edges: "):
                edges.append(int(line.rstrip("\n").split(" ")[1]))
                pass
            elif line.__contains__("Tot Weight: "):
                weights.append(int(line.rstrip("\n").split(" ")[2]))
                pass
            elif line.__contains__("W Hist: "):
                hists.append([int(v) for v in line.rstrip("\n").split(" ")[3:8]])
                next_graph = True
                pass
            elif next_SDA:
                next_SDA = False
                SDA.append(line)
                cont_SDA = True
                pass
            elif cont_SDA:
                SDA.append(line)
                pass
            elif next_graph:
                next_graph = False
                network.append(line)
                cont_graph = True
                pass
            elif cont_graph:
                network.append(line)
                pass
            pass
        networks.append(network)
        pass

    samps = cnt
    p0_list = get_patient_zero(dir_path)
    # run number, fitness, profileS, dnaS, edge count
    if len(fits) != samps:
        print("ERROR in fits: " + dir_path)
        pass
    if len(networks) != samps:
        print("ERROR in networks: " + dir_path)
        pass
    if len(SDAs) != samps:
        print("ERROR in SDAs: " + dir_path)
        pass
    if len(p0_list) != samps:
        print("ERROR in p0_list: " + dir_path)
        pass
    if len(edges) != samps:
        print("ERROR in edges: " + dir_path)
        pass
    if len(weights) != samps:
        print("ERROR in weights: " + dir_path)
        pass
    if len(hists) != samps:
        print("ERROR in hists: " + dir_path)
        pass

    data = [[i + 1, fits[i], SDAs[i], networks[i], edges[i], weights[i], hists[i], p0_list[i]] for i in range(samps)]
    data.sort(key=itemgetter(1))  # Ascending
    if "ED" in dir_path:
        data.reverse()
        pass
    return data


def get_patient_zero(dir_path: str):
    p0_list = []
    with open(dir_path + "patient0.dat") as f:
        lines = f.readlines()
        prev = None
        for line in lines:
            if line.__contains__("mating event"):
                li = line.split(" ")
                prev = int(li[4])
                pass
            if line == '\n':
                p0_list.append(prev)
                pass
            pass
        pass
    return p0_list


def main():
    print("START")
    folder_names = os.listdir(inp)
    mode_itms = [["ED"], ["PM", "P1"], ["PM", "P7"], ["PM", "DUB"]]
    mode_info = ["ED", "PM1", "PM7", "DUBPM"]
    sizes = [128, 256]
    mode_dirs = [[[] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    muts = [1, 2, 3]
    states = [12, 16, 20]
    nets_to_grab_per_exp = [[[1, 2, 3], [3, 8, 9]], [[3, 5, 6], [1, 2, 3]], [[1, 6, 9], [1, 6, 9]],
                            [[1, 6, 8], [1, 6, 8]]]

    exp_lbls = []
    exp_dat = []
    exp_descriptions = []
    exp_idx = 1
    for s in states:
        for m in muts:
            exp_dat.append([str(s) + "S", str(m) + "M"])
            exp_lbls.append(str(exp_idx) + "(" + str(s) + ", " + str(m) + ")")
            exp_descriptions.append("Exp " + str(exp_idx) + ": " + str(s) + " States, " + str(m) + " Muts")
            exp_idx += 1
            pass
        pass

    for midx, itms in enumerate(mode_itms):
        for sidx, size in enumerate(sizes):
            for dat in exp_dat:
                for fld in folder_names:
                    if all(itm in fld for itm in itms) and str(size) in fld and all(d in fld for d in dat):
                        mode_dirs[midx][sidx].append(fld)
                        pass
                    elif all(itm in fld for itm in itms) and all(str(si) not in fld for si in sizes) and \
                            all(d in fld for d in dat):
                        mode_dirs[midx][sidx].append(fld)
                        pass
                    pass
                pass
            pass
        pass

    # mode_data[mode][size][exp][run][0] = run num
    # mode_data[mode][size][exp][run][1] = run's fit (sorted based on this)
    # mode_data[mode][size][exp][run][2][:] = run's SDA
    # mode_data[mode][size][exp][run][3][:] = run's network
    # mode_data[mode][size][exp][run][4] = run's network's edges
    # mode_data[mode][size][exp][run][5] = run's network's weight
    # mode_data[mode][size][exp][run][6][:] = run's network's weight hist
    # mode_data[mode][size][exp][run][7] = run's network's best p0
    mode_data = [[[] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    for midx, itms in enumerate(mode_itms):
        for sidx, size in enumerate(sizes):
            for fld in mode_dirs[midx][sidx]:
                mode_data[midx][sidx].append(get_data(inp + fld + "/"))
                pass
            pass
        pass

    for midx, itms in enumerate(mode_itms):
        for sidx, size in enumerate(sizes):
            file_name = ""
            if midx < 3:
                file_name = str(size) + mode_info[midx]
                pass
            elif midx == 3 and sidx == 0:
                file_name = mode_info[midx]
                pass
            if midx < 3 or sidx == 0:
                for exp in nets_to_grab_per_exp[midx][sidx]:
                    with open(outp + file_name + "_BestNet" + str(exp) + ".dat", "w") as f:
                        f.write(exp_descriptions[exp - 1] + "\n")
                        f.write("Best Fitness: " + str(mode_data[midx][sidx][exp - 1][0][1]) + "\n")
                        f.write("Representation: directed linked lists, repetition shows integer weight, and blank line "
                                "before next experiment\n")
                        f.write("Note: if node 0 and 1 share an edge, these appear in both lists\n")
                        f.write("Patient Zero: " + str(mode_data[midx][sidx][exp - 1][0][7]) + "\n")
                        f.write("Network Starts on Next Line:\n")
                        if midx < 3:
                            f.write(str(size) + " " + str(mode_data[midx][sidx][exp - 1][0][4]) + " " +
                                str(mode_data[midx][sidx][exp - 1][0][5]) + "\n")
                        else:
                            f.write("200 " + str(mode_data[midx][sidx][exp - 1][0][4]) + " " +
                                str(mode_data[midx][sidx][exp - 1][0][5]) + "\n")
                            pass
                        f.writelines(mode_data[midx][sidx][exp - 1][0][3])
                    pass
                pass
            pass
        pass
    print("END")
    pass


main()
