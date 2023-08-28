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

def get_dist(data_path, acc_path, gen_path, t_branch_data, t_branch_acc, t_branch_gen, min_t, max_t, n_bins, gen_branches):
    with uproot.open(data_path) as data_file, uproot.open(acc_path) as acc_file, uproot.open(gen_path) as gen_file:
        data_tree = data_file[data_file.keys()[0]]
        acc_tree = acc_file[acc_file.keys()[0]]
        gen_tree = gen_file[gen_file.keys()[0]]
        t_filter = lambda branch_name: f"({branch_name} > {min_t}) & ({branch_name} < {max_t})"
        data_t = data_tree.arrays([t_branch_data], t_filter(t_branch_data), library='np')[t_branch_data]
        acc_t = acc_tree.arrays([t_branch_acc], t_filter(t_branch_acc), library='np')[t_branch_acc]
        gen_df = gen_tree.arrays([t_branch_gen, *gen_branches], t_filter(t_branch_gen), library='np')
        gen_t = gen_df[t_branch_gen]
        data_hist, bins = np.histogram(data_t, range=(min_t, max_t), bins=n_bins)
        acc_hist, _ = np.histogram(acc_t, range=(min_t, max_t), bins=n_bins)
        gen_hist, _ = np.histogram(gen_t, range=(min_t, max_t), bins=n_bins)
        ratio = data_hist / acc_hist
        normalized_ratio = ratio / np.abs(sum(ratio) * np.diff(bins)[0])
        Cg = normalized_ratio  / np.amax(normalized_ratio);
        accepted_ids = set()
        for index, t_val in tqdm(enumerate(gen_t), total=len(gen_t)):
            t_bin = np.digitize(t_val, bins) - 1
            prob = Cg[t_bin]
            if np.random.uniform(0, 1) <= prob:
                accepted_ids.add("_".join([str(gen_df[branch][index]) for branch in gen_branches]))
        return accepted_ids

def resample(accepted_ids, gen_path, acc_path, gen_branches, acc_branches):
    gen_file = ROOT.TFile.Open(str(gen_path))
    gen_tree = gen_file.GetListOfKeys().At(0).ReadObj()
    gen_out_file = ROOT.TFile(str(gen_path.parent / (gen_path.stem + "_sampled.root")), "RECREATE")
    gen_out_tree = gen_tree.CloneTree(0)
    print("Resampling Generated MC...")
    for entry in tqdm(range(gen_tree.GetEntries())):
        gen_tree.GetEntry(entry)
        event_id = "_".join([str(getattr(gen_tree, branch)) for branch in gen_branches])
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
        event_id = "_".join([str(getattr(acc_tree, branch)) for branch in acc_branches])
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
                            args['--t-data'],
                            args['--t-acc'],
                            args['--t-gen'],
                            float(args['-l']),
                            float(args['-h']),
                            int(args['-n']),
                            args['-g'].split(','))
    resample(accepted_ids,
             Path(args['<gen>']),
             Path(args['<acc>']),
             args['-g'].split(','),
             args['-a'].split(','))

if __name__ == '__main__':
    main()
