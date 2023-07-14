# MandelSample
Resample some MC according to the ratio of data to reconstructed MC. Credit to Lawrence Ng for the original idea!

## Requirements
See the [requirements file](requirements.txt) for the Python requirements. This package also requires a [ROOT/PyROOT installation](https://root.cern/).

## Installation
1. Clone the repository: `git clone git@github.com:denehoffman/MandelSample.git`
2. Install the package `pip install MandelSample/.` (or `cd MandelSample && pip install .`)
3. You now have access to the `mandelsample` command (assuming your ROOT installation works with the Python version in which you installed this package)

## Usage
```shell
$ mandelsample --help
MandelSample
Usage: mandelsample [options] <data> <acc> <gen>

Options:
    --help                         Show this screen.
    --t-branch <t_branch>          Specify branch with Mandelstam t [default: t]
    --run-branch <run_branch>      Specify branch with run number [default: RunNumber]
    --event-branch <event_branch>  Specify branch with event number [default: EventNumber]
    -l <min_t>                     Set lower bound of t-distribution [default: -1.0]
    -h <max_t>                     Set upper bound of t-distribution [default: 0.0]
    -n <n_bins>                    Set number of bins to use [default: 50]
```
