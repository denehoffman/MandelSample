#!/usr/bin/env python3
"""MandelSample
Usage: mandelsample [options] <data> <acc> <gen> 

Options:
    --help                                  Show this screen.
    --t-branch-data <t_branch_data>         Specify branch with Mandelstam t for Data [default: t]
    --t-branch-acc <t_branch_acc>           Specify branch with Mandelstam t for Accepted Monte Carlo [default: t]
    --t-branch-gen <t_branch_gen>           Specify branch with Mandelstam t for Generated Monte Carlo [default: t_p2p4]
    --run-branch-acc <run_branch_acc>       Specify branch with run number for Accepted Monte Carlo [default: run]
    --event-branch-acc <event_branch_acc>   Specify branch with event number for Accepted Monte Carlo [default: event]
    --run-branch-gen <run_branch_gen>       Specify branch with run number for Generated Monte Carlo [default: RunNumber]
    --event-branch-gen <event_branch_gen>   Specify branch with event number for Generated Monte Carlo [default: EventNumber]
    -l <min_t>                              Set lower bound of t-distribution [default: 0.1]
    -h <max_t>                              Set upper bound of t-distribution [default: 2.0]
    -n <n_bins>                             Set number of bins to use [default: 25]
"""


import numpy as np
from pathlib import Path
import uproot
from docopt import docopt
import matplotlib.pyplot as plt
import ROOT
from tqdm import tqdm

def get_dist(data_path, acc_path, gen_path, t_branch_data, t_branch_acc, t_branch_gen, min_t, max_t, n_bins, run_branch_gen, event_branch_gen):
    with uproot.open(data_path) as data_file, uproot.open(acc_path) as acc_file, uproot.open(gen_path) as gen_file:
        data_tree = data_file[data_file.keys()[0]]
        acc_tree = acc_file[acc_file.keys()[0]]
        gen_tree = gen_file[gen_file.keys()[0]]
        data_t = data_tree.arrays([t_branch_data], f"({t_branch_data} > {min_t}) & ({t_branch_data} < {max_t})", library='np')[t_branch_data]
        acc_t = acc_tree.arrays([t_branch_acc], f"({t_branch_acc} > {min_t}) & ({t_branch_acc} < {max_t})", library='np')[t_branch_acc]
        gen_df = gen_tree.arrays([t_branch_gen, run_branch_gen, event_branch_gen], f"({t_branch_gen} > {min_t}) & ({t_branch_gen} < {max_t})", library='np')
        gen_t = gen_df[t_branch_gen]
        gen_run = gen_df[run_branch_gen]
        gen_event = gen_df[event_branch_gen]
        data_hist, bins = np.histogram(data_t, range=(min_t, max_t), bins=n_bins)
        acc_hist, _ = np.histogram(acc_t, range=(min_t, max_t), bins=n_bins)
        gen_hist, _ = np.histogram(gen_t, range=(min_t, max_t), bins=n_bins)
        ratio = data_hist / acc_hist
        normalized_ratio = ratio / np.abs(sum(ratio) * np.diff(bins)[0])
        #fig, ax = plt.subplots(figsize=(15, 7))
        #fig = plt.figure(constrained_layout=True)
        #fig, ((ax0, ax1), (ax2)) plt.subplot_mosaic()
        #axs = fig.subplot_mosaic([['TopLeft', 'Right'],['BottomLeft', 'Right']],
        #                  gridspec_kw={'width_ratios':[1, 2]})
        #axs['Right'].hist(bins[:-1], bins, weights=normalized_ratio, color='tab:red')
        #axs['Right'].set_title('Normalized Ratio')
        #axs['TopLeft'].hist(data_t, density=True, range=(min_t, max_t), bins=n_bins)
        #axs['TopLeft'].set_title('Data Hist')
        #axs['BottomLeft'].hist(acc_t, density=True, range=(min_t, max_t), bins=n_bins)
        #axs['BottomLeft'].set_title('Acc Hist')
        #ax3.set_title("normalized_ratio")     
        #ax.stairs(normalized_ratio, bins, label='Data/Accepted MC', color='C0')
        #ax.legend(loc="upper left")
        #ax.set_ylabel("Normalized Ratio")
        #ax2 = ax.twinx()
        #ax2.stairs(data_hist, bins, label='Data', color='C1')
        #ax2.stairs(acc_hist, bins, label='Accepted MC', color='C2')
        #ax2.legend(loc="upper right")
        #bin_width = round(np.abs(max_t - min_t) / n_bins, 2)
        #ax2.set_ylabel(f"Counts / {bin_width} GeV${{}}^2$")
        #plt.xlabel("Mandelstam t")
        #plt.tight_layout()
        #fig.tight_layout()
        #plt.show()
        Cg = normalized_ratio  / np.amax(normalized_ratio);
        print("Evaluating probabilities...")
        accepted_ids = set()
        for index, t_val in tqdm(enumerate(gen_t), total=len(gen_t)):
            t_bin = np.digitize(t_val, bins) - 1
            prob = Cg[t_bin]
            if np.random.uniform(0, 1) <= prob:
                accepted_ids.add(f"{gen_run[index]}: {gen_event[index]}")
        return accepted_ids

def resample(accepted_ids, gen_path, acc_path, run_branch_gen, event_branch_gen, run_branch_acc, event_branch_acc):
    gen_file = ROOT.TFile.Open(str(gen_path))
    gen_tree = gen_file.GetListOfKeys().At(0).ReadObj()
    gen_out_file = ROOT.TFile(str(gen_path.parent / (gen_path.stem + "_sampled.root")), "RECREATE")
    gen_out_tree = gen_tree.CloneTree(0)
    print("Resampling Generated MC...")
    for entry in tqdm(range(gen_tree.GetEntries())):
        gen_tree.GetEntry(entry)
        run_number = getattr(gen_tree, run_branch_gen)
        event_number = getattr(gen_tree, event_branch_gen)
        event_id = f"{run_number}: {event_number}"
        if event_id in accepted_ids:
            gen_out_tree.Fill()
    gen_out_tree.Write()
    gen_out_file.Close()
    gen_file.Close()

    # now that we have the resampled generated MC, we need to remove the events that were rejected in the accepted MC

    acc_file = ROOT.TFile.Open(str(acc_path))
    acc_tree = acc_file.GetListOfKeys().At(0).ReadObj()
    acc_out_file = ROOT.TFile(str(acc_path.parent / (acc_path.stem + "_sampled.root")), "RECREATE")
    acc_out_tree = acc_tree.CloneTree(0)
    print("Resampling Accepted MC...")
    print(len(accepted_ids))
    print(acc_tree.GetEntries())
    for entry in tqdm(range(acc_tree.GetEntries())):
        acc_tree.GetEntry(entry)
        run_number = getattr(acc_tree, run_branch_acc)
        event_number = getattr(acc_tree, event_branch_acc)
        event_id = f"{run_number}: {event_number}"
        if event_id in accepted_ids:
            acc_out_tree.Fill()
    acc_out_tree.Write()
    acc_out_file.Close()
    acc_file.Close()


def main():
    args = docopt(__doc__)
    accepted_ids = get_dist(Path(args['<data>']),
                            Path(args['<acc>']),
                            Path(args['<gen>']),
                            args['--t-branch-data'],
                            args['--t-branch-acc'],
                            args['--t-branch-gen'],
                            float(args['-l']),
                            float(args['-h']),
                            int(args['-n']),
                            args['--run-branch-gen'],
                            args['--event-branch-gen'])
    resample(accepted_ids,
             Path(args['<gen>']),
             Path(args['<acc>']),
             args['--run-branch-gen'],
             args['--event-branch-gen'],
             args['--run-branch-acc'],
             args['--event-branch-acc'])

if __name__ == '__main__':
    main()
