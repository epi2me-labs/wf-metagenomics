"""Calculate different alpha diversity metrics."""

# Each different taxa within a rank is consider
# as an approach of one OTU.

import numpy as np
import scipy as scp

# Rarefy
#


def random_choice(arr, sample_size, seed=4):
    """Random sample n elements from the total of inviduals.

    :param arr: vector with number of counts per taxa. Type: numpy array.
    :param sample_size: number of elements (n) to be randomly sampled. Type: int.
    :param seed (int): seed value.
    :return: vector with rarefied counts.
    :rtype: numpy array.
    """
    np.random.seed(seed)
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


def global_rarefaction(df, min_counts=0.95, pc_reads=0.5, seed=4):
    """Rarefy all samples to the same size.

    :param df: pandas dataframe with counts per taxon and sample.
    :param min_counts (float): percentage of sequences to take for
        the minimum number of counts in the samples.
    :param pc_reads (float): percentage of sequences relative to the
        median number of reads of all the samples.
    :param seed (int): seed value.
    :return: dataframe with rarefied counts.
        All samples have the same number of counts.
    :rtype: pandas dataFrame.

    From all the samples, take the one with
    the lower number of counts and rarefy
    to the 95% of it.
    This is the number of sequences to take for all the samples
    """
    # Get median number of reads of all the samples
    median_nreads = df.sum().median()
    # Filter those samples with less than x% reads than average
    threshold = int(pc_reads * median_nreads)
    samples = df.columns[df.sum() >= threshold]
    df_to_rarefy = df[samples]
    rarefaction_size = int(min_counts * min(df_to_rarefy.sum()))
    rarefied_df = df_to_rarefy.apply(lambda x: random_choice(
        x, rarefaction_size, seed=seed), axis=0)
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
    if np.isnan(shannon_index):
        return float("nan")
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
    if counts.sum() == 0:
        return float("nan")
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
    if counts.sum() == 0:
        return float("nan")
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


def berger_parkers(counts):
    """Calculate Berger–Parker index.

    BP = max(pi)

    :param counts: vector with number of counts per taxa. Type: pandas Series.
    :return: berger parkers index.
    :rtype: float.
    """
    if counts.sum() == 0:
        return float("nan")
    return max(counts/counts.sum())


def fishers_alpha(counts):
    """Calculate Fisher's slope constant.

    S = alpha * ln(1 + N/alpha)
    :param counts: vector with number of counts per taxa. Type: pandas Series.
    :rtype: float.
    """
    s = len(counts)
    n = counts.sum()
    if n == 0:
        return float("nan")

    def eq(alpha):
        return alpha * np.log(1 + n/alpha) - s
    alpha_fisher = scp.optimize.fsolve(eq, 1)[0]
    return alpha_fisher


# Diversity metrics

def alpha_diversity(df):
    """Calculate alpha diversity metrics.

    :param df: pandas dataframe with counts per taxa and sample.
    :return: dataframe with some diversity metrics.
    :rtype: pandas DataFrame.
    """
    # df expects samples in columns, taxa in rows
    # Read table
    diversity = df.apply(sum, axis=0, raw=False).to_frame()
    indices_dict = {
        "S": "Richness", "H": "Shannon diversity index",
        "ENS": "Effective number of species", "D": "Simpson\'s index",
        "Inverse of D": "Inverse Simpson\'s index", "J": "Pielou\'s evenness",
        "F": "Fisher\'s alpha", "BP": "Berger Parker index"
        }
    diversity.columns = ['Total counts']
    # Richness: n of Species in the ecosystem
    diversity['S'] = df.apply(observed_richness, axis=0)
    diversity['H'] = df.apply(shannon_entropy, axis=0)
    diversity['ENS'] = diversity[
        'H'].apply(effective_nspecies)
    diversity['D'] = df.apply(simpson, axis=0)
    diversity['Inverse of D'] = diversity[
        'D'].apply(lambda x: 1/x if x != 0 else float("nan"))
    diversity['BP'] = df.apply(berger_parkers, axis=0)
    diversity['F'] = df.apply(fishers_alpha, axis=0)
    # Evenness: the extent to which species are evenly distributed
    diversity['J'] = diversity.apply(lambda x: pielou(x.H, x.S), axis=1)
    diversity.rename(columns=indices_dict, inplace=True)
    return diversity.T
