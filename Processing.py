import math
import os
from operator import itemgetter

import matplotlib.pyplot as plt
from graphviz import Graph
import numpy as np

inp = "../../Conferences and Papers/2023 CIBCB/WSDA/BitOutDone/"
outp = "./BitProcessed/"
james_inp = "./BitOutput/"
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


def get_base_data(file_path: str):
    data = []
    with open(file_path) as f:
        vals = []
        lines = f.readlines()
        for line in lines:
            if line == "\n":
                data.append(vals)
                vals = []
                pass
            elif not line.__contains__("EE"):
                vals.append(float(line))
                pass
            pass
        data.append(vals)
        pass
    return data


def get_james_data():
    data = []
    # ED
    mode_data = []
    exp_dat = get_base_data(james_inp + "ED128results.dat")
    mode_data.append(exp_dat)
    exp_dat = get_base_data(james_inp + "ED256results.dat")
    mode_data.append(exp_dat)
    data.append(mode_data)

    # PM 1
    mode_data = []
    exp_dat = get_base_data(james_inp + "PM128Prof1results.dat")
    mode_data.append(exp_dat)
    exp_dat = get_base_data(james_inp + "PM256Prof1results.dat")
    mode_data.append(exp_dat)
    data.append(mode_data)

    # PM 7
    mode_data = []
    exp_dat = get_base_data(james_inp + "PM128Prof7results.dat")
    mode_data.append(exp_dat)
    exp_dat = get_base_data(james_inp + "PM256Prof7results.dat")
    mode_data.append(exp_dat)
    data.append(mode_data)

    # DUB
    mode_data = []
    exp_dat = get_base_data(james_inp + "PMdublinResults.dat")
    mode_data.append(exp_dat)
    data.append(mode_data)
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


def writeStat(data: [], out):
    mean = float(np.mean(data))
    mean = round(mean, precision)
    std = float(np.std(data, ddof=0))
    std = round(std, precision)  # Population standard deviation
    diff = 1.96 * std / math.sqrt(30)  # 95% CI
    diff = round(diff, precision)
    if lower_better:
        maxima = float(min(data))
        pass
    else:
        maxima = float(max(data))
        pass
    maxima = round(maxima, precision)
    out.write(str(mean).ljust(col_width))
    out.write(str(std).ljust(col_width))
    out.write(str(diff).ljust(col_width))
    out.write(str(maxima).ljust(col_width))
    return mean, maxima


def make_table(many_data: [], exp_info: [], fname: str):
    with open(fname, "w") as f:
        f.write("EXP".ljust(col_width))
        f.write("Parameters".ljust(2 * col_width))
        f.write("Best Run".ljust(col_width))
        f.write("Mean".ljust(col_width))
        f.write("SD".ljust(col_width))
        f.write("95% CI".ljust(col_width))
        f.write("Best".ljust(col_width))
        f.write("\n")
        for di, data in enumerate(many_data):
            f.write(str("EXP" + str(di)).ljust(col_width))
            f.write(exp_info[di].ljust(2 * col_width))
            f.write(str("Run " + str(data[0][0])).ljust(col_width))
            d = [data[i][1] for i in range(samps)]
            writeStat(d, f)
            f.write("\n")
            pass
        pass
    pass


def write_best(data: [], info: str, seed: str):
    with open(seed, "w") as f:
        f.write("Run " + str(data[0][0]) + "\n")
        f.write(str(data[0][1]) + "\n")
        f.write(info + "\n")
        profs = data[0][2]
        for prof in profs:
            pstr = ""
            idx = 0
            while idx < prof[0]:
                pstr += "\t"
                idx += 1
                continue
            for p in prof[2]:
                pstr += str(p) + "\t"
                pass
            f.write(pstr)
            f.write("\n")
            pass
        dnas = data[0][3]
        for dna in dnas:
            f.write(dna)
            f.write("\n")
            pass
        pass
    pass


def box_plot(bp, skip: int, edge_color):
    base_colors = ["#808080", "#696969"]
    colors = ['#0000FF', '#00FF00', '#FFFF00']
    for whisker in bp['whiskers']:
        whisker.set(color='#8B008B', linewidth=1)
        pass

    for cap in bp['caps']:
        cap.set(color='#8B008B', linewidth=1)
        pass

    for median in bp['medians']:
        median.set(color='red', linewidth=1)
        pass

    for flier in bp['fliers']:
        flier.set(marker='.', color='#e7298a', alpha=0.5, markersize=3)

    for idx, patch in enumerate(bp['boxes']):
        if idx < skip:
            patch.set(facecolor=base_colors[idx % len(base_colors)])
        else:
            patch.set(facecolor=colors[idx % len(colors)])
        pass
    pass


def inf_data(data):
    spreads = []

    for run in data:
        spread = 0
        for prof in run[2]:
            spread += prof[4]
            pass
        spreads.append(spread)
        pass
    return np.mean(spreads), spreads[0]


def calc(data):
    mean = float(np.mean(data))
    std = float(np.std(data, ddof=0))
    diff = 1.96 * std / math.sqrt(30)  # 95% CI
    print(str(mean) + "+-" + str(diff))
    pass


def high_low_deg(el: [], verts: int):
    deg = [int(0) for _ in range(verts)]
    low_deg = []
    high_deg = []
    for ed in el:
        deg[ed[0]] += 1
        deg[ed[1]] += 1
        pass
    max = 0
    for idx, deg in enumerate(deg):
        if deg > max:
            max = deg
            pass
        if deg > 20:
            high_deg.append(idx)
            pass
        elif deg > 10:
            low_deg.append(idx)
            pass
        pass
    print("Max: " + str(max))
    return low_deg, high_deg


def edge_list(lines):
    el = []
    lines.__delitem__(-1)
    for from_node, line in enumerate(lines):
        line = line.rstrip()
        line = line.split("\t")
        for to_node in line:
            if to_node != '':
                if [from_node, int(to_node)] not in el:
                    if [int(to_node), from_node] not in el:
                        el.append([from_node, int(to_node)])
                        pass
                    pass
                pass
            pass
        pass

    edge_lists = []
    for d in el:
        if d not in edge_lists:
            edge_lists.append(d)
            pass
        pass
    return edge_lists


def make_graph(el: [], low_deg: [], high_deg: [], out_file: str, verts: int):
    g = Graph(engine='sfdp')
    e_cout = 0

    g.graph_attr.update(dpi='1000', size="10,10", outputorder='edgesfirst', overlap='false', splines='true')
    g.node_attr.update(color='black', shape='point', width='0.02', height='0.02')
    g.edge_attr.update(color='black', penwidth='0.2')
    for n in range(verts):
        if n == 0:
            if n in low_deg:
                g.node(str(n), label=str(n), color='red', width='0.03', height='0.03')
                pass
            elif n in high_deg:
                g.node(str(n), label=str(n), color='red', width='0.04', height='0.04')
                pass
            else:
                g.node(str(n), label=str(n), color='red')
                pass
        elif n in low_deg:
            g.node(str(n), label=str(n), width='0.03', height='0.03')
        elif n in high_deg:
            g.node(str(n), label=str(n), width='0.04', height='0.04')
        else:
            g.node(str(n), label=str(n))
        pass

    for idx, d in enumerate(el):
        if d[0] < d[1]:
            g.edge(str(d[0]), str(d[1]), color='black')
            pass
        e_cout += 1
        pass
    print("Edges: " + str(e_cout))
    g.render(filename=out_file, directory=outp, cleanup=True, format='png')

    # g.save(filename=out_file, directory=outp, format='png')
    pass


def get_settings(dir: str):
    pass


def make_from_dirs(dirs: [], splits: [], titles: []):
    exp_lbls = []
    data = []

    split_dirs = []
    for spl in splits:
        dir_list = []
        for dir in dirs:
            if spl in dir:
                dir_list.append(dir)
                pass
            pass
        split_dirs.append(dir_list)
        pass

    for split in split_dirs:
        lbl_list = []
        split_data = []
        for dir in split:
            lbl_list.append(dir[:10])
            split_data.append(get_data(inp + dir + "/"))
            pass
        exp_lbls.append(lbl_list)
        data.append(split_data)
        pass

    # mode_stats[mode][size][exp] = [run's fitness vals]
    just_fits = []
    for split in data:
        split_fits = []
        for exp in split:
            exp_fits = []
            for run in exp:
                exp_fits.append(run[1])
                pass
            split_fits.append(exp_fits)
            pass
        just_fits.append(split_fits)
        pass

    for split_idx in range(len(splits)):
        plt.rc('xtick', labelsize=8)
        plt.rc('ytick', labelsize=8)

        f = plt.figure()
        f.set_dpi(900)
        f.set_figheight(8)
        f.set_figwidth(10)

        plot = f.add_subplot(111)

        bp = plot.boxplot(just_fits[split_idx], patch_artist=True)
        box_plot(bp, "#FFFFFF")

        # plot.set_xticks(xpos[idx])
        plot.set_xticklabels(exp_lbls[split_idx], rotation=90)

        f.suptitle(titles[split_idx], fontsize=12)
        plot.set_xlabel("Experiment", fontsize=10)
        plot.set_ylabel("Fitness", fontsize=10)
        # for x in lxpos:
        #     plot.axvline(x=x, color='black', linestyle='--', linewidth=0.75)
        #     pass
        plot.grid(visible="True", axis="y", which='major', color="darkgray", linewidth=0.75)
        f.tight_layout()
        f.savefig(outp + titles[split_idx] + ".png", dpi=900)
        plt.close()
        pass
    pass


def main():
    print("START")
    folder_names = os.listdir(inp)
    # mode_itms = [["ED"]]
    mode_itms = [["ED"], ["PM", "P1"], ["PM", "P7"], ["PM", "DUB"]]
    # sizes = [128, 256, 384, 512]
    sizes = [128]
    mode_dirs = [[[] for _ in range(len(sizes))] for _ in range(len(mode_itms))]  # ED[size], PM1[size], PM7[size], DUB[size]
    muts = [1, 2, 3]
    # states = [12]
    states = [12, 16, 20]

    # plt.style.use("seaborn-dark")
    # make_from_dirs(folder_names, ["ED", "PM"], ["Epidemic Length Tests", "Profile 1 Matching Tests"])

    exp_lbls = ["Base1", "Base2"]
    exp_dat = []
    exp_idx = 1
    for s in states:
        for m in muts:
            exp_dat.append([str(s) + "S", str(m) + "M"])
            exp_lbls.append(str(exp_idx) + "(" + str(s) + ", " + str(m) + ")")
            exp_idx += 1
            pass
        pass

    for midx, itms in enumerate(mode_itms):
        for sidx, size in enumerate(sizes):
            for dat in exp_dat:
                for fld in folder_names:
                    if all(itm in fld for itm in itms) and str(size) not in fld and all(d in fld for d in dat):
                        mode_dirs[midx][sidx].append(fld)
                        pass
                    if all(itm in fld for itm in itms) and str(size) in fld and all(d in fld for d in dat):
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

    # base_stats[mode][size][exp] = [run's fitness vals]
    base_stats = get_james_data()

    # mode_stats[mode][size][exp] = [run's fitness vals]
    mode_stats = [[[] for _ in range(len(sizes))] for _ in range(len(mode_itms))]
    for midx, itms in enumerate(mode_itms):
        for sidx, size in enumerate(sizes):
            mode_stats[midx][sidx].append(base_stats[midx][sidx][0])
            mode_stats[midx][sidx].append(base_stats[midx][sidx][1])
            for exp in mode_data[midx][sidx]:
                exp_fits= []
                for run in exp:
                    exp_fits.append(run[1])
                    pass
                mode_stats[midx][sidx].append(exp_fits)
                pass
            pass
        pass

    # titles = ["Epidemic Length"]
    titles = ["Epidemic Length", "Epidemic Profile Matching P1", "Epidemic Profile Matching P7", "Epidemic Profile Matching Dublin"]
    names = ["EL_boxplot", "PM1_boxplot", "PM7_boxplot", "DUB_boxplot"]
    # xsp = [[i for i in range(len(all_data[0]))], [i for i in range(len(all_data[1]))]]
    # xpos = [xsp[0], xsp[1], xsp[0], xsp[1], xsp[0], xsp[1], xsp[0], xsp[1]]
    ylb = ["Fitness", "Fitness", "Fitness", "Fitness"]
    xlb = ["Experiment (Num. States, Max. Mutations)",
           "Experiment (Num. States, Max. Mutations)",
           "Experiment (Num. States, Max. Mutations)",
           "Experiment (Num. States, Max. Mutations)"]

    # lxpos = []
    # for i in range(2, len(all_data[0]) - 3, 3):
    #     lxpos.append(i + 0.5)
    #     pass

    for idx in range(len(titles)):
        for sidx, size in enumerate(sizes):
            plt.rc('xtick', labelsize=6)
            plt.rc('ytick', labelsize=6)

            f = plt.figure()
            f.set_figheight(5)
            f.set_figwidth(8)
            plot = f.add_subplot(111)

            bp = plot.boxplot(mode_stats[idx][sidx], patch_artist=True)
            box_plot(bp, 2, "#FFFFFF")

            # plot.set_xticks(xpos[idx])
            plot.set_xticklabels(exp_lbls, rotation=90)

            if not titles[idx].__contains__("Dublin"):
                f.suptitle(titles[idx] + " w " + str(size) + " Nodes", fontsize=12)
            else:
                f.suptitle(titles[idx], fontsize=12)
                pass
            plot.set_xlabel(xlb[idx], fontsize=10)
            plot.set_ylabel(ylb[idx], fontsize=10)
            # for x in lxpos:
            #     plot.axvline(x=x, color='black', linestyle='--', linewidth=0.75)
            #     pass
            plot.grid(visible="True", axis="y", which='major', color="darkgray", linewidth=0.75)
            f.tight_layout()
            f.savefig(outp + names[idx] + ".png", dpi=300)
            plt.close()
            pass
        pass


    # plt.rc('xtick', labelsize=6)
    # plt.rc('ytick', labelsize=6)
    # for fidx, many_data in enumerate(all_data):
    #     make_table(many_data, x_lbl[fidx], outp + mode[fidx] + "_table.dat")
    #     for didx, data in enumerate(many_data):
    #         profs = data[0][2]
    #         mean_inf, best_inf = inf_data(data)
    #         f = plt.figure()
    #         f.set_dpi(600)
    #         f.set_figheight(4)
    #         f.set_figwidth(8)
    #
    #         plot = f.add_subplot(111)
    #
    #         xs = []
    #         ys = []
    #         infs = 0
    #         for prof in profs:
    #             xs.append(prof[2])
    #             ys.append([i for i in range(prof[0], prof[1] + 1)])
    #             infs += prof[4]
    #             pass
    #
    #         for idx in range(len(profs)):
    #             if idx == 0:
    #                 lbl = profs[idx][3]
    #             else:
    #                 lbl = profs[idx][3]
    #             plot.plot(ys[idx], xs[idx], label=lbl)
    #             pass
    #
    #         y_max = plot.get_ylim()[1]
    #         x_max = plot.get_xlim()[1]
    #         plot.set_ylim(0, y_max)
    #         plot.set_xlim(0, x_max)
    #
    #         f.suptitle(titles[fidx] + " Experiment " + str(didx + 1) + " Best Epidemic Curve", fontsize=12)
    #         sp_info = "Total Infections: " + str(infs) + '    '
    #         sp_info += "Mean Infections: " + str(format(mean_inf, '.2f'))
    #         plot.text(plot.get_xlim()[1] / 2, plot.get_ylim()[1], sp_info, horizontalalignment='center',
    #                   verticalalignment='bottom', bbox=dict(fc='white', ec='none', pad=0), fontsize=10)
    #         plot.set_xlabel("Day", fontsize=10)
    #         plot.set_ylabel("New Infections", fontsize=10)
    #         plot.minorticks_on()
    #         plot.grid(visible="True", axis="x", which='minor', color="darkgray", linewidth=0.5)
    #         plot.grid(visible="True", axis="x", which='major', color="black", linewidth=0.5)
    #         plot.grid(visible="True", axis="y", which='major', color="black", linewidth=0.5)
    #         plt.legend(title="Variant ID (Parent)", borderaxespad=0, bbox_to_anchor=(1, 0.5),
    #                    loc="center left", fontsize=6, title_fontsize=6, ncol=2)
    #         f.tight_layout()
    #         f.savefig(outp + mode[fidx] + "_EXP" + str(didx + 1) + "_profile.png", dpi=600)
    #         plt.close()
    #         write_best(data, x_lbl[fidx][didx], outp + mode[fidx] + "_EXP" + str(didx + 1) + "_best.dat")
    #         pass
    #     pass
    #
    # plt.rc('xtick', labelsize=6)
    # plt.rc('ytick', labelsize=6)
    # for fidx, many_data in enumerate(all_data):
    #     for didx, data in enumerate(many_data):
    #         profs = data[0][2]
    #         mean_inf, best_inf = inf_data(data)
    #         f = plt.figure()
    #         f.set_dpi(600)
    #         f.set_figheight(4)
    #         f.set_figwidth(8)
    #
    #         plot = f.add_subplot(111)
    #
    #         # xs = []
    #         # ys = []
    #         # infs = 0
    #         # for prof in profs:
    #         #     xs.append(prof[2])
    #         #     ys.append([i for i in range(prof[0], prof[1] + 1)])
    #         #     infs += prof[4]
    #         #     pass
    #
    #         ys = data[0][6]
    #         cont = True
    #         idx = 127
    #         while cont:
    #             if ys[idx] == 0:
    #                 ys = ys[0:-1]
    #                 pass
    #             else:
    #                 cont = False
    #                 pass
    #             idx -= 1
    #             pass
    #
    #         xs = [i for i in range(0, len(ys))]
    #         x_lbl = [str(i) for i in xs]
    #         plot.bar(xs, ys, label=x_lbl)
    #
    #         y_max = plot.get_ylim()[1]
    #         x_max = plot.get_xlim()[1]
    #         plot.set_ylim(0, y_max)
    #         plot.set_xlim(0, x_max)
    #
    #         f.suptitle(titles[fidx] + " Experiment " + str(didx + 1) + " Best Epidemic Histogram", fontsize=12)
    #         sp_info = "Total Infections: " + str(best_inf) + '    '
    #         sp_info += "Mean Infections: " + str(format(mean_inf, '.2f'))
    #         plot.text(plot.get_xlim()[1] / 2, plot.get_ylim()[1], sp_info, horizontalalignment='center',
    #                   verticalalignment='bottom', bbox=dict(fc='white', ec='none', pad=0), fontsize=10)
    #         plot.set_xlabel("Severity of Infection", fontsize=10)
    #         plot.set_ylabel("Number of Infections", fontsize=10)
    #         # plot.minorticks_on()
    #         # plot.grid(visible="True", axis="x", which='minor', color="darkgray", linewidth=0.5)
    #         # plot.grid(visible="True", axis="x", which='major', color="black", linewidth=0.5)
    #         # plot.grid(visible="True", axis="y", which='major', color="black", linewidth=0.5)
    #         # plt.legend(title="Variant ID (Parent)", borderaxespad=0, bbox_to_anchor=(1, 0.5),
    #         #            loc="center left", fontsize=10, title_fontsize=10)
    #         f.tight_layout()
    #         f.savefig(outp + mode[fidx] + "_EXP" + str(didx + 1) + "_hist.png", dpi=600)
    #         plt.close()
    #         # write_best(data, x_lbl[fidx][didx], outp + mode[fidx] + "_EXP" + str(didx) + "_best.dat")
    #         # print("From File: " + str(data[0][4]))
    #         # el = edge_list(data[0][7])
    #         # low_deg, high_deg = high_low_deg(el, 256)
    #         # make_graph(el, low_deg, high_deg, mode[fidx] + "_EXP" + str(didx + 1) + "_graph", 256)
    #         pass
    #     pass

    print("END")
    pass


main()
