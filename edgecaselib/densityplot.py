from sys import stdout
from numpy import linspace, array, mean, nan, concatenate, fromstring, full, vstack
from pandas import read_csv, DataFrame, concat
from matplotlib.pyplot import switch_backend, subplots
from tqdm import tqdm
from os import path


def binned(A, bins, func=mean):
    """Return array data compressed into bins (smoothed by func)"""
    coords = linspace(0, len(A), bins+1).astype(int)
    return array([
        func(A[start:end])
        for start, end in zip(coords, coords[1:])
    ])


def get_binned_density_dataframe(raw_densities, chrom, bin_size):
    """Subset to densities of reads on `chrom`, bin densities, convert to dataframe"""
    indexer = (raw_densities["chrom"] == chrom)
    densities_subset = raw_densities[indexer].reset_index(drop=True)
    leftmost_pos = (densities_subset["pos"] - densities_subset["clip_5prime"]).min()
    # convert density strings into numpy array:
    density_arrays = []
    entry_iterator = tqdm(
        densities_subset.itertuples(), total=densities_subset.shape[0],
        desc="Binning densities ({})".format(chrom), unit="entry"
    )
    # align density data to reference and bin it:
    for entry in entry_iterator:
        unaligned_density = fromstring(entry.density, sep=",", dtype="float32")
        padder = full(entry.pos - entry.clip_5prime - leftmost_pos, nan)
        aligned_density = concatenate([padder, unaligned_density])
        binned_density = binned(
            aligned_density, bins=aligned_density.shape[0]/bin_size
        )
        density_arrays.append(binned_density)
    # pad densities on the right so they are all the same length:
    max_density_length = max(
        binned_density.shape[0] for binned_density in density_arrays
    )
    for i, binned_density in enumerate(density_arrays):
        padder = full(max_density_length - binned_density.shape[0], nan)
        density_arrays[i] = concatenate([binned_density, padder])
    naked_binned_density_dataframe = DataFrame(
        data=vstack(density_arrays),
        columns=[leftmost_pos + j * bin_size for j in range(max_density_length)]
    )
    return concat(
        [densities_subset.iloc[:,:-1], naked_binned_density_dataframe], axis=1
    )


def load_densities(dat, bin_size):
    """Load densities from dat file, split into dataframes per chromosome"""
    raw_densities = read_csv(dat, compression="gzip", sep="\t", escapechar="#")
    return {
        chrom: get_binned_density_dataframe(raw_densities, chrom, bin_size)
        for chrom in raw_densities["chrom"].drop_duplicates()
    }


def plot_densities(densities, figsize, palette, hide_names, bin_size, xtick_density, title="", file=stdout.buffer):
    """Plot binned densities as a heatmap"""
    raise NotImplementedError
    if file is not None:
        switch_backend("Agg")
    width, height = tuple(map(int, figsize.split("x")))
    figure, ax = subplots(figsize=(width, height))
    # adjust xticks according to bin size and alignment:
    xticks = [
        int(tick.get_text()) for tick in ax.get_xticklabels()
    ]
    if xtick_density != 1:
        each = int(len(xticks)*xtick_density)
        xticks = [
            xtick if i%each==0 else ""
            for i, xtick in enumerate(xticks)
        ]
    ax.set_xticklabels(
        [tick * bin_size for tick in xticks],
        rotation=60
    )
    # drop yticks (useful for long names / big sets):
    if hide_names:
        ax.set(yticks=[])
    ax.set(title=title, xlabel="position (bp)")
    figure.savefig(file)


def main(dat, bin_size=100, figsize="16x9", palette="viridis", hide_names=False, title=None, xtick_density=.05, file=stdout.buffer, **kwargs):
    """Dispatch data to subroutines"""
    densities = load_densities(dat, bin_size=bin_size)
    if title is None:
        title = path.split(dat)[-1]
    plot_densities(
        densities, figsize, palette, hide_names, bin_size,
        xtick_density, title, file
    )
    return 0
