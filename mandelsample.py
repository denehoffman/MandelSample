#!/usr/bin/env python3
"""MandelSample
Usage: mandelsample [options] <data> <acc> <gen> 

Options:
    --t-data <t_data>   Specify branch with Mandelstam t for Data [default: t]
    --t-acc <t_acc>     Specify branch with Mandelstam t for Accepted MC (defaults to <t_data>)
    --t-gen <t_gen>     Specify branch with Mandelstam t for Generated MC (defaults to <t_data>)
    -a <branches>       Specify list of branches to identify Accepted MC [default: RunNumber,EventNumber]
    -g <branches>       Specify list of branches to identify Generated MC [default: RunNumber,EventNumber]
    -l <min_t>          Set lower bound of t-distribution [default: 0.1]
    -h <max_t>          Set upper bound of t-distribution [default: 2.0]
    -n <n_bins>         Set number of bins to use [default: 25]
    --help              Show this screen.
"""


import numpy as np
from pathlib import Path
import uproot
from docopt import docopt
import matplotlib.pyplot as plt
import ROOT
from tqdm import tqdm

def main():
    args = docopt(__doc__)
    data_path = Path(args['<data>'])
    acc_path = Path(args['<acc>'])
    gen_path = Path(args['<gen>'])
    t_low = args['-l']
    t_high = args['-h']
    n_bins = args['-n']
    with uproot.open(data_path) as data_file, uproot.open(acc_path) as acc_file, uproot.open(gen_path) as gen_file:
        data_tree = data_file[data_file.keys()[0]]
        acc_tree = acc_file[acc_file.keys()[0]]
        gen_tree = gen_file[gen_file.keys()[0]]
        t_filter = lambda branch_name: f"({branch_name} > {t_low}) & ({branch_name} < {t_high})"
        data_t = data_tree.arrays([args['--t-data']], t_filter(args['--t-data']), library='np')[args['--t-data']]
        acc_t = acc_tree.arrays([args['--t-acc']], t_filter(args['--t-acc']), library='np')[args['--t-acc']]
        gen_df = gen_tree.arrays([args['--t-gen'], *args['-g'].split(',')], t_filter(args['--t-gen']), library='np')
        gen_t = gen_df[args['--t-gen']]
        data_hist, bins = np.histogram(data_t, range=(t_low, t_high), bins=n_bins)
        acc_hist, _ = np.histogram(acc_t, range=(t_low, t_high), bins=n_bins)
        gen_hist, _ = np.histogram(gen_t, range=(t_low, t_high), bins=n_bins)
        ratio = data_hist / acc_hist
        normalized_ratio = ratio / np.abs(sum(ratio) * np.diff(bins)[0])
        Cg = normalized_ratio  / np.amax(normalized_ratio);
        accepted_ids = set()
        for index, t_val in tqdm(enumerate(gen_t), total=len(gen_t)):
            t_bin = np.digitize(t_val, bins) - 1
            prob = Cg[t_bin]
            if np.random.uniform(0, 1) <= prob:
                accepted_ids.add("_".join([str(gen_df[branch][index]) for branch in args['-g'].split(',')]))

    gen_file = ROOT.TFile.Open(str(gen_path))
    gen_out_file = ROOT.TFile(str(gen_path.parent / (gen_path.stem + "_sampled.root")), "RECREATE")
    try:
        gen_tree = gen_file.GetListOfKeys().At(0).ReadObj()
        gen_out_tree = gen_tree.CloneTree(0)
        print("Resampling Generated MC...")
        for entry in tqdm(range(gen_tree.GetEntries())):
            gen_tree.GetEntry(entry)
            event_id = "_".join([str(getattr(gen_tree, branch)) for branch in args['-g'].split(',')])
            if event_id in accepted_ids:
                gen_out_tree.Fill()
        gen_out_tree.Write()
    finally:
        gen_out_file.Close()
        gen_file.Close()

    # now that we have the resampled generated MC, we need to remove the events that were rejected in the accepted MC

    acc_file = ROOT.TFile.Open(str(acc_path))
    acc_out_file = ROOT.TFile(str(acc_path.parent / (acc_path.stem + "_sampled.root")), "RECREATE")
    try:
        acc_tree = acc_file.GetListOfKeys().At(0).ReadObj()
        acc_out_tree = acc_tree.CloneTree(0)
        print("Resampling Accepted MC...")
        print(len(accepted_ids))
        print(acc_tree.GetEntries())
        for entry in tqdm(range(acc_tree.GetEntries())):
            acc_tree.GetEntry(entry)
            event_id = "_".join([str(getattr(acc_tree, branch)) for branch in args['-a'].split(',')])
            if event_id in accepted_ids:
                acc_out_tree.Fill()
        acc_out_tree.Write()
    finally:
        acc_out_file.Close()
        acc_file.Close()

if __name__ == '__main__':
    main()
