edgeCase
========

## ./kmer-scan.py

Performs calculation of density of a given kmer in a rolling window along each
read.  
Description of command-line options is available by running
`python3 kmer-scan.py -h`.

#### By default, performs two passes:

1. Calculates densities in terminal windows (either heads or tails, depending on
whether --head-test xor --tail-test is specified with the size of the terminal
window). These densities are then fit to a Gaussian mixture model to separate
background level of densities and a distribution of significant densities
(the clustering can be stored by supplying --output-gmm).
2. For each read, tests if the density of the terminal window belongs to the
distribution of significant densities (with a p-value of 1e-5, which can be
modified by passing --pmax), and if so, calculates and outputs densities in the
rolling window along the read.

This allows to select telomeric candidates based on how they relate to the reads
in the entire sample, regardless of absolute density or true "repeatedness" of
the kmer. The advantage of this approach is that it does not require a priori
knowledge of how conserved the repeat is; however, if the entire sample has low
density of the target kmer, it will yield false positives.

#### If --cutoff is specified instead, performs only one pass:

1. Calculates densities in terminal windows (either heads or tails, depending on
whether --head-test xor --tail-test is specified with the size of the terminal
window). For any given read, if this density exceeds the specified hard cutoff,
calculates and outputs densities in the rolling window along the read.

This allows to select telomeric candidates based on an a priori expected
density. The advantage of this approach is that we are still able to select
candidate reads regardless of how much like a repeat their ends look; however,
we are biasing the search by providing a hard cutoff.

#### If neither --cutoff nor --head-test nor --tail-test are specified, performs only one pass:

1. Calculates and outputs densities in the rolling window along each read.

This is useful for cross-validation of results from other tools.


## ./plot-metric.py

Plots heatmap-like visualizations of the metric (e.g. density) along candidate
reads.  
Description of command-line options is available by running
`python3 plot-metric.py -h`.

The input file must be in the output format of kmer-scan.py, tab-separated:

```
read_name    value    value    value   value    value    ...
read_name    value    value    value   value    value    ...
read_name    value    value    value   value    value    ...
read_name    value    value    value   value    value    ...
```
