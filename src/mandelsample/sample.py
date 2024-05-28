#!/usr/bin/env python3
"""MandelSample
Usage: mandelsample [options] <data> <acc> <gen> 

Options:
    --t-data <t_data>    Specify branch with Mandelstam t for Data [default: t]
    --t-acc <t_acc>      Specify branch with Mandelstam t for Accepted MC (defaults to <t_data>)
    --t-gen <t_gen>      Specify branch with Mandelstam t for Generated MC (defaults to <t_data>)
    -a <branches>        Specify list of branches to identify Accepted MC [default: RunNumber,EventNumber]
    -g <branches>        Specify list of branches to identify Generated MC [default: RunNumber,EventNumber]
    -l <min_t>           Set lower bound of t-distribution [default: 0.1]
    -h <max_t>           Set upper bound of t-distribution [default: 2.0]
    -n <n_bins>          Set number of bins to use [default: 25]
    -b --binning <type>  Specify binning type ('sqrt' or 'rice') [default: sqrt]
    --help               Show this screen.

Notes:
    1. The options '-a' and '-g' are used to construct unique identifiers for each event. The default
       will only work if the combination of RunNumber and EventNumber is unique for every event, and
       maps to the same combination of numbers in both Accepted and Generated MC. 

       If these alone do not uniquely identify an event, you can add others with those arguments. The
       branch names can even differ between the TTrees as long as the underlying data is the same.
       For instance, in the current versions of MCWrapper, the EventNumber is reset for every batch,
       and the batches are all merged before the user has access to the individual files. This means
       there will be ~N events with the same RunNumber and EventNumber where N is the number of batches.
    
       To circumvent this, you could additionally specify a branch such as the thrown beam energy, 
       which should be approximately unique to that set of ~N events in this case. The numeric value
       of the branches are converted to strings separated by '_'. For example:

       $ mandelsample data.root acc.root gen.root -a RunNumber,EventNumber,E_ThrownBeam -g RunNumber,EventNumber,E_Beam
    
       This would create unique IDs like "37040_1403_7.834902" (or something similar).
    2. The values for -l <min_t> and -h <max_t> are used with respect to "-t", so you'll probably
       want to use positive numbers with <min_t> < <max_t>.
"""


import sys
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
    t_low = float(args['-l'])
    t_high = float(args['-h'])
    if t_low >= t_high:
        print(f"Error! Bounds on -t are unphysical ({t_low} ≮ {t_high})")
        sys.exit(1)
    
    custom_n_bins = args['-n']
    binning_type = args['--binning']
    
    t_data = args['--t-data']
    t_acc = args['--t-acc']
    if not t_acc:
        t_acc = t_data  # default to t_data if not given
    t_gen = args['--t-gen']
    if not t_gen:
        t_gen = t_data  # default to t_data if not given
    gen_branches = args['-g'].split(',')
    acc_branches = args['-a'].split(',')
    if len(gen_branches) != len(acc_branches):
        print("Error! You must specify the same number of branches in both Accepted and Generated MC for unique identification of events")
        sys.exit(1)
    
    # Open all three files with uproot
    with uproot.open(data_path) as data_file, uproot.open(acc_path) as acc_file, uproot.open(gen_path) as gen_file:
        # Grab first TTree in each file
        data_tree = data_file[data_file.keys()[0]]
        acc_tree = acc_file[acc_file.keys()[0]]
        gen_tree = gen_file[gen_file.keys()[0]]

        # Get desired branches from each file, selecting t in specified bounds
        print("Loading datasets... \n")
        t_filter = lambda branch_name: f"({branch_name} > {t_low}) & ({branch_name} < {t_high})"
        data_t = data_tree.arrays([t_data], t_filter(t_data), library='np')[t_data]
        acc_t = acc_tree.arrays([t_acc], t_filter(t_acc), library='np')[t_acc]
        gen_df = gen_tree.arrays([t_gen, *gen_branches], t_filter(t_gen), library='np')
        gen_t = gen_df[t_gen]

        # Dynamically calculate and store the number of bins using info from the data
        if custom_n_bins:
            n_bins_data = int(custom_n_bins)
        elif binning_type == 'rice':
            n_bins_data = int(np.ceil(2 * np.cbrt(len(data_t))))
        else:  # Default to sqrt if no valid binning type is provided
            n_bins_data = int(np.ceil(np.sqrt(len(data_t))))
        
        # Make some histograms of t for each dataset
        data_hist, bins = np.histogram(data_t, range=(t_low, t_high), bins=n_bins_data)
        acc_hist, _ = np.histogram(acc_t, range=(t_low, t_high), bins=n_bins_data)

        # Get probability distribution
        pdf = data_hist / acc_hist

        # Account for instances when inf is found in the pdf when there are 0 events for acc_hist 
        if np.isinf(pdf).any():
            print("Infinite values found in PDF, setting them to zero \n")
            pdf[np.isinf(pdf)] = 0

        # Scale pdf so that max(pdf) = 1.0 (this speeds up sampling)
        scaled_pdf = pdf / np.amax(pdf);

        # Create a set of unique "IDs" for the sampled events
        accepted_ids = set()

        # Sample generated MC
        print("Calculating sampled events...")
        for index, t_val in tqdm(enumerate(gen_t), total=len(gen_t)):
            # shortcut for getting the bin where this t would be placed
            t_bin = np.digitize(t_val, bins) - 1
            # sample the point if the scaled_efficiency of the bin exceeds a uniform random number
            if np.random.uniform(0, 1) <= scaled_pdf[t_bin]:
                # IDs are formated as an underscore-separated list of values
                # For example, the default might combine RunNumber and EventNumber
                # to formulate IDs like 37040_13, 37040_14, etc.
                accepted_ids.add("_".join([str(gen_df[branch][index]) for branch in gen_branches]))

    # We will write the sampled files in ROOT because it's faster and we can copy the whole
    # tree structure dynamically with CloneTree(0)
    gen_file = ROOT.TFile.Open(str(gen_path))
    gen_out_file = ROOT.TFile(str(gen_path.parent / (gen_path.stem + "_" + binning_type + "_sampled.root")), "RECREATE")
    try:
        # Get the first tree in the file
        gen_tree = gen_file.GetListOfKeys().At(0).ReadObj()
        gen_out_tree = gen_tree.CloneTree(0)
        print("Resampling Generated MC...")
        for entry in tqdm(range(gen_tree.GetEntries())):
            gen_tree.GetEntry(entry)
            event_id = "_".join([str(getattr(gen_tree, branch)) for branch in gen_branches])
            # if the ID is in the set of IDs from the random sampling, we fill the event in the output tree
            if event_id in accepted_ids:
                gen_out_tree.Fill()
        gen_out_tree.Write()
    finally:
        gen_out_file.Close()
        gen_file.Close()

    # Now that we have the resampled generated MC, we need to remove the rejected events in the accepted MC
    acc_file = ROOT.TFile.Open(str(acc_path))
    acc_out_file = ROOT.TFile(str(acc_path.parent / (acc_path.stem  + "_" + binning_type + "_sampled.root")), "RECREATE")
    try:
        # Get the first tree in the file
        acc_tree = acc_file.GetListOfKeys().At(0).ReadObj()
        acc_out_tree = acc_tree.CloneTree(0)
        print("Resampling Accepted MC...")
        print(len(accepted_ids))
        print(acc_tree.GetEntries())
        for entry in tqdm(range(acc_tree.GetEntries())):
            acc_tree.GetEntry(entry)
            event_id = "_".join([str(getattr(acc_tree, branch)) for branch in acc_branches])
            # if the ID is in the set of IDs from the random sampling, we fill the event in the output tree
            if event_id in accepted_ids:
                acc_out_tree.Fill()
        acc_out_tree.Write()
    finally:
        acc_out_file.Close()
        acc_file.Close()

if __name__ == '__main__':
    main()
