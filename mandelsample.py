#!/usr/bin/env python3
"""MandelSample
Usage: mandelsample [options] <data> <acc> <gen> 

Options:
    --help                         Show this screen.
    --t-branch <t_branch>          Specify branch with Mandelstam t [default: t]
    --run-branch <run_branch>      Specify branch with run number [default: RunNumber]
    --event-branch <event_branch>  Specify branch with event number [default: EventNumber]
    -l <min_t>                     Set lower bound of t-distribution [default: -1.0]
    -h <max_t>                     Set upper bound of t-distribution [default: 0.0]
    -n <n_bins>                    Set number of bins to use [default: 50]
"""


import numpy as np
from pathlib import Path
import uproot
from docopt import docopt
import matplotlib.pyplot as plt
import ROOT
from tqdm import tqdm

def get_dist(data_path, acc_path, gen_path, t_branch, min_t, max_t, n_bins, run_branch, event_branch):
    with uproot.open(data_path) as data_file, uproot.open(acc_path) as acc_file, uproot.open(gen_path) as gen_file:
        data_tree = data_file[data_file.keys()[0]]
        acc_tree = acc_file[acc_file.keys()[0]]
        gen_tree = gen_file[gen_file.keys()[0]]
        data_t = data_tree[t_branch].array(library='np')
        acc_t = acc_tree[t_branch].array(library='np')
        gen_t = gen_tree[t_branch].array(library='np')
        gen_run = gen_tree[run_branch].array(library='np')
        gen_event = gen_tree[event_branch].array(library='np')
        data_hist, bins = np.histogram(data_t, range=(min_t, max_t), bins=n_bins)
        acc_hist, _ = np.histogram(acc_t, range=(min_t, max_t), bins=n_bins)
        gen_hist, _ = np.histogram(gen_t, range=(min_t, max_t), bins=n_bins)
        ratio = data_hist / acc_hist
        normalized_ratio = ratio / np.abs(sum(ratio) * np.diff(bins)[0])
        fig, ax = plt.subplots(figsize=(15, 7))
        ax.stairs(normalized_ratio, bins, label='Data/Accepted MC', color='C0')
        ax.legend(loc="upper left")
        ax.set_ylabel("Normalized Ratio")
        ax2 = ax.twinx()
        ax2.stairs(data_hist, bins, label='Data', color='C1')
        ax2.stairs(acc_hist, bins, label='Accepted MC', color='C2')
        ax2.legend(loc="upper right")
        bin_width = round(np.abs(max_t - min_t) / n_bins, 2)
        ax2.set_ylabel(f"Counts / {bin_width} GeV${{}}^2$")
        plt.xlabel("Mandelstam t")
        plt.tight_layout()
        plt.show()
        const = 1.0
        for ibin, value in enumerate(normalized_ratio):
            gen_value = gen_hist[ibin]
            proposal = gen_value / value
            if proposal > const:
                const = proposal
        Cg = normalized_ratio * const
        f = gen_hist
        print("Evaluating probabilities...")
        accepted_ids = set()
        for index, t_val in tqdm(enumerate(gen_t)):
            t_bin = np.digitize(t_val, bins) - 1
            prob = f[t_bin] / Cg[t_bin]
            if np.random.uniform(0, 1) <= prob:
                accepted_ids.add(f"{gen_run[index]}: {gen_event[index]}")
        return accepted_ids

def resample(accepted_ids, gen_path, acc_path, run_branch, event_branch):
    gen_file = ROOT.TFile.Open(str(gen_path))
    gen_tree = gen_file.GetListOfKeys().At(0).ReadObj()
    gen_out_file = ROOT.TFile(str(gen_path.parent / (gen_path.stem + "_sampled.root")), "RECREATE")
    gen_out_tree = gen_tree.CloneTree(0)
    print("Resampling Generated MC...")
    for entry in tqdm(range(gen_tree.GetEntries())):
        gen_tree.GetEntry(entry)
        run_number = getattr(gen_tree, run_branch)
        event_number = getattr(gen_tree, event_branch)
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
    for entry in tqdm(range(acc_tree.GetEntries())):
        acc_tree.GetEntry(entry)
        run_number = getattr(acc_tree, run_branch)
        event_number = getattr(acc_tree, event_branch)
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
                            args['--t-branch'],
                            float(args['-l']),
                            float(args['-h']),
                            int(args['-n']),
                            args['--run-branch'],
                            args['--event-branch'])
    resample(accepted_ids,
             Path(args['<gen>']),
             Path(args['<acc>']),
             args['--run-branch'],
             args['--event-branch'])
