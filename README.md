# MandelSample
Resample some MC according to the ratio of data to reconstructed MC. Credit to Lawrence Ng for the original idea!

# Theory
When generating Monte Carlo (MC), we want the physics to match as closely as possible to our data. However, it is often the case that the data don't have some nice smooth exponential distribution in Mandelstam `t`. That's where MandelSample comes in. The idea is that we create a probability density function (PDF) based on the `t` from the data and use this PDF to resample the Monte Carlo. This can be thought of mathematically as follows:

```math
G^* = G \frac{D}{A}
```
where $G^*$ is the resampled generated MC and $G$, $A$, and $D$ are the original generated MC, accepted MC, and data distributions, respectively. Since we are sampling from the generated MC already, we only need to calculate a PDF proportional to $\frac{D}{A}$. This fraction scales inversely with the size of the accepted MC, and that means that if we take it as given, it will usually create very small probabilities in each bin in `t`. Therefore, we divide by the maximum across all bins to ensure that the sampling PDF's largest bin is equal to $1$ and the other bins are scaled appropriately. Another way to think about this is that we want to sample the largest subset of the MC as possible, and if there is a value of `t` which maximizes the PDF, we always want to sample events with that value (or similar/within a bin). The fraction $\frac{D}{A}$ effectively gives the proportional difference between the data and Monte Carlo.

For simplicity, suppose we reduce our problem to just two bins in `t`. In the data, we have $20$ events with $t\in[0.1, 0.2]$ and $40$ events with $t\in[0.2, 0.3]$. In the accepted MC, these values are $30$ and $10$ respectively, and in the generated, they are $300$ and $100$. The fractions are then $\frac{20}{30} = 0.66$ and $\frac{40}{10} = 4.0$. Dividing by the max gives $0.165$ and $1.0$. Now when we resample, we get generated bins with approximately $300 * 0.165 \sim 50$ and $100 * 1.0 \sim 100$, and this makes the resampled accepted MC have bins with $5$ and $10$ events, which is the same proportion as the data. This is an idealized example, of course, but it should demonstrate the general idea. Note that while we've significantly reduced the total amount of MC ($400$ generated and $40$ accepted), we've done it in a maximal way, since all of the events in the second bin, which was the largest in the data, were accepted.


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
    --t-data <t_data>   Specify branch with Mandelstam t for Data [default: t]
    --t-acc <t_acc>     Specify branch with Mandelstam t for Accepted MC (defaults to <t_data>)
    --t-gen <t_gen>     Specify branch with Mandelstam t for Generated MC (defaults to <t_data>)
    -a <branches>       Specify list of branches to identify Accepted MC [default: RunNumber,EventNumber]
    -g <branches>       Specify list of branches to identify Generated MC [default: RunNumber,EventNumber]
    -l <min_t>          Set lower bound of t-distribution [default: 0.1]
    -h <max_t>          Set upper bound of t-distribution [default: 2.0]
    -n <n_bins>         Set number of bins to use [default: 25]
    --help              Show this screen.

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
```
