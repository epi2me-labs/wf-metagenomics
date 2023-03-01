"""Calculate different alpha diversity metrics."""

# Each different taxa within a rank is consider
# as an approach of one OTU.

import numpy as np

# General
#

BRACKEN_TO_RANK = {
    'S': 'species',
    'G': 'genus',
    'F': 'family',
    'O': 'order',
    'C': 'class',
    'P': 'phylum'
}


def extract_rank(d, rank_dict, rank):
    """Extract count for a rank level. Return {tx:count}.

    :param d: nested dictionary coming from bracken json:
        {"Sample": {"Taxon name":{"rank":str, "count":int, "children": dict},}}.
    :param rank_dict: dictionary {taxa:counts}.
    :param rank: String indicating the taxonomic rank at
        which perform the counts extraction.
    :return: dictionary {taxa:counts}.
    :rtype: dict.
    """
    if BRACKEN_TO_RANK.get(rank):
        parsed_rank = BRACKEN_TO_RANK[rank]
    else:
        parsed_rank = rank

    for tx, v in d.items():
        if d[tx]['rank'] == parsed_rank:
            rank_dict[tx] = d[tx]['count']
        else:
            extract_rank(d[tx]['children'], rank_dict, rank=parsed_rank)
    return rank_dict


# Rarefy
#

def random_choice(arr, sample_size):
    """Random sample n elements from the total of inviduals.

    :param arr: vector with number of counts per taxa. Type: numpy array.
    :param sample_size: number of elements (n) to be randomly sampled. Type: int.
    :return: vector with rarefied counts.
    :rtype: numpy array.
    """
    arr = np.int64(arr)
    # To avoid working with strings, index the positions
    arr_idx = np.arange(len(arr))
    # Extend numpy vector with their corresponding counts
    long_arr = np.repeat(arr_idx, arr)
    # Take sample_size elements and count frequencies
    rarefied_vector = np.bincount(
        np.random.choice(long_arr, sample_size, replace=False))
    # If the last element of arr has not been subseted,
    # the bin count has one element less length.
    # Add a 0 counts for it
    if len(arr) > len(rarefied_vector):
        rarefied_vector = np.append(rarefied_vector, np.repeat(
            0, len(arr)-len(rarefied_vector)))
    return rarefied_vector


def global_rarefaction(df, min_counts=0.95):
    """Rarefy all samples to the same size.

    :param df: pandas dataframe with counts per taxon and sample.
    :param min_counts: percentage of sequences to take for
        the minimum number of counts in all the samples.
    :type min_counts: float. Value from 0 to 1.
    :return: dataframe with rarefied counts.
        All samples have the same number of counts.
    :rtype: pandas dataFrame.

    From all the samples, take the one with
    the lower number of counts and rarefy
    to the 95% of it.
    This is the number of sequences to take for all the samples
    """
    size = int(min_counts*min(df.sum()))
    rarefied_df = df.apply(lambda x: random_choice(
        np.array(x), size), axis=0)
    return rarefied_df


def rarefaction_curve(counts, n_repeats=10, n_points=15):
    """Count richness (S) per set of subsample reads.

    :param counts: vector with number of counts per taxa. Type: pandas Series.
    :param n_repeats: number of repeats to estimate richness. Type: int.
    :param n_points: number of points to divide the interval
        to estimate richness. Type: int.
    :return: dictionary {number of sampled reads (str): richness (int)}.
    :rtype: dict.
    """
    arr = np.array(counts)
    total_counts = arr.sum()
    richness = {}
    # Take n counts from the total number of counts.
    # Do it in intervals saved in the lists steps.
    # Divide each interval in 15 points.
    # do not step more than available counts.
    steps = np.linspace(0, total_counts, num=n_points, dtype=int)
    for s in steps:  # iterate over number of counts to be sampled
        # Repeat n times and return mean S
        spc = np.zeros(n_repeats, dtype=int)
        for i in range(n_repeats):
            spc[i] = round(observed_richness(
                random_choice(np.array(counts), s)))
        richness[int(s)] = round(np.mean(spc))
    return richness

# Alpha-diversity functions
#


def observed_richness(counts):
    """Calculate species richness.

    i.e. the total number of unique species.

    :param counts: vector with number of counts per taxa.
    :type counts: numpy array or pandas Series.
    :return: number of different taxa.
    :rtype: int.
    """
    species = np.count_nonzero(counts)
    return species


def effective_nspecies(shannon_index):
    """Calculate effective number of Species.

    Effective number of species = exp(H)

    :param shannon_index: shannon diversity index. Type: float.
    :return: effective number of species.
    :rtype: float.
    """
    return np.exp(shannon_index)


def shannon_entropy(counts):
    """Compute the Shannon entropy of a vector.

    H = -sum(pi*ln(pi))

    :param counts: vector with number of counts per taxa. Type: pandas Series.
    :return: Shannon entropy of a vector.
    :rtype: float.
    """
    # pi = unique sequences counts* /
    #       total counts in the sample
    #       for each i species
    # * we don't have it with kraken, is an approach
    # log base = math.e
    p = counts/counts.sum()
    logp = np.where(p > 0, np.log(p), 0)
    shannon_index = abs(np.dot(p, logp))
    return shannon_index


def simpson(counts):
    """Calculate Simpson Diversity Index.

    D = 1 - sum(pi**2) for each i-species.

    :param counts: vector with number of counts per taxa. Type: pandas Series.
    :return: simpson diversity index.
    :rtype: float.
    """
    p = (counts/counts.sum())**2
    simpson_index = 1 - p.sum()
    return simpson_index


def pielou(shannon_index, species):
    """Calculate Pielou Evenness.

    J = H/log(S)

    :param shannon_index: shannon diversity index. Type: float.
    :param species: richness (number of species). Type: float.
    :return: Pielou evennes index.
    :rtype: float.
    """
    if species <= 0:
        return float("nan")
    elif species == 1:
        return 0
    else:
        return shannon_index/np.log(species)


# Diversity metrics
#


def alpha_diversity(df):
    """Calculate alpha diversity metrics.

    :param df: pandas dataframe with counts per taxa and sample.
    :return: dataframe with some diversity metrics.
    :rtype: pandas DataFrame.
    """
    # df expects samples in columns, taxa in rows
    # Read table
    diversity = df.apply(sum, axis=0, raw=False).to_frame()
    diversity.columns = ['Total counts']
    # Richness: n of Species in the ecosystem
    diversity['S'] = df.apply(observed_richness, axis=0)
    diversity['H'] = df.apply(shannon_entropy, axis=0)
    diversity['Effective N Species'] = diversity[
        'H'].apply(effective_nspecies)
    diversity['Simpson\'s Index (D)'] = df.apply(simpson, axis=0)
    diversity['Inverse of D'] = diversity[
        'Simpson\'s Index (D)'].apply(lambda x: 1/x if x != 0 else float("nan"))
    # Evenness: the extent to which species are evenly distributed
    diversity['Pielou evenness'] = diversity.apply(lambda x: pielou(x.H, x.S), axis=1)
    return diversity
