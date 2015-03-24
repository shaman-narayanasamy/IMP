"""VizBin Python wrapper

Usage:
    vizbin.py INPUT OUTPUT DIMENSIONS ... [options]

Description:
    INPUT       The input fasta file.
    OUTPUT      The output prefix path ("dimension" is added to the file).
    DIMENSIONS  The dimensions reductions (could have several).

Options:
    -h --help           Show this screen.
    --version           Show version.
    -l LENGTH           The length of the kmers to compute [default: 5].
    -s SIZE             Size of the thread pool [default: 1].
    -t THETA            Theta [default: 0.5].
    -p PERP             Perp [default: 30]
    -f MINLEN           Filter sequences < MINLEN [default: 1000].
    -o FILTERED_FASTA   Output the filtered fasta file
"""

import itertools
import numpy as np
from collections import defaultdict
import sys
import math
import multiprocessing
from multiprocessing import Manager
from sklearn.decomposition import PCA
from docopt import docopt


def get_product(alphabet, len_word):
    for word in itertools.product(alphabet, repeat=len_word):
        yield ''.join(word)


def isheader(line):
    return line[0] == '>'


def read_fasta(fasta_path, minlen, filtered_fasta):

    # function to read the fasta file that yield the sequences
    def read(rhandle, whandle, fn):
        ident = None
        prev_head = None
        for head, group in itertools.groupby(rhandle, key=isheader):
            if not head:
                seq = format_sequence(group)
                if len(seq) >= minlen:
                    fn(whandle, prev_head, seq)
                    yield seq
                else:
                    sequences_excluded.append(ident)
            else:
                prev_head = format_sequence(group)
                # ident = format_sequence(group)[1:]

        if whandle is not None:
            print("[x] Wrote filtered fasta file to: %s" % filtered_fasta)
        print("[x] %s sequences excluded." % len(sequences_excluded))

    # write function is write_filtered fasta option is given
    def write_fa(whandle, prev_head, seq):
        whandle.write('%s\n' % format_sequence(prev_head))
        whandle.write('%s\n' % ''.join([i.strip() for i in seq]))

    # null function
    def null(whandle, prev_head, seq):
        pass

    # main
    format_sequence = lambda x: ''.join([i.strip() for i in x])
    sequences_excluded = []

    if filtered_fasta is not None:
        with open(filtered_fasta, 'w') as whandle:
            with open(fasta_path, 'r') as rhandle:
                for s in read(rhandle, whandle, write_fa):
                    yield s
    else:
        with open(fasta_path, 'r') as rhandle:
            for s in read(rhandle, whandle, null):
                yield s

    # if filtered_fasta is not None:
    #     with open(filtered_fasta, 'w') as whandle:
    #         format_sequence = lambda x: ''.join([i.strip() for i in x])
    #         sequences_excluded = []
    #         with open(fasta_path, 'r') as rhandle:
    #             ident = None
    #             prev_head = None
    #             for head, group in itertools.groupby(rhandle, key=isheader):
    #                 if head:
    #                     prev_head = group
    #                 if not head:
    #                     seq = format_sequence(group)
    #                     if len(seq) >= minlen:
    #                         whandle.write('%s\n' % format_sequence(prev_head))
    #                         whandle.write('%s\n' % ''.join([i.strip() for i in seq]))
    #                         yield seq
    #                     else:
    #                         sequences_excluded.append(ident)
    #                 else:
    #                     ident = format_sequence(group)[1:]
    #         print("[x] Wrote filtered fasta file to: %s" % filtered_fasta)
    #         print("[x] %s sequences excluded." % len(sequences_excluded))
    # else:
    #     format_sequence = lambda x: ''.join([i.strip() for i in x])
    #     sequences_excluded = []
    #     with open(fasta_path, 'r') as rhandle:
    #         ident = None
    #         for head, group in itertools.groupby(rhandle, key=isheader):
    #             if not head:
    #                 seq = format_sequence(group)
    #                 if len(seq) >= minlen:
    #                     yield seq
    #                 else:
    #                     sequences_excluded.append(ident)
    #             else:
    #                 ident = format_sequence(group)[1:]
    #     print("[x] %s sequences excluded." % len(sequences_excluded))


def window(seq, n=2):
    it = iter(seq)
    result = tuple(itertools.islice(it, n))
    if len(result) == n:
        yield ''.join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield ''.join(result)


def get_reverse_complement(word):
    comp = {
        'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'
    }
    return ''.join(reversed([comp[letter] for letter in word]))


def get_words_without_complement(words):
    result = []
    mapper = {}
    for word in words:
        rev = get_reverse_complement(word)
        if rev in result:
            mapper[word] = rev
        elif word in result:
            mapper[word] = word
        else:
            result.append(word)
            mapper[word] = word
    return result, mapper


def update_progress(progress, message="Completion"):
    barLength = 100  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength * progress))
    text = "\r{0}: [{1}] {2:.0f}% {3}".format(message, "#"*block + "-"*(barLength-block), progress * 100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


def setup_progress_bar(toolbar_width):
    sys.stdout.write("[%s]" % (" " * toolbar_width))
    sys.stdout.flush()
    sys.stdout.write("\b" * (toolbar_width+1))  # return to start of line, after '['


def row_frequency(index, frequencies, mapper, sequence, lword, words, words_without_complement, excluded_kmers):
    nb_kmers = 0
    # apply a sliding window and add +1 to the word found in the counter
    for w in window(sequence, LEN):
        if w not in words:
            excluded_kmers.append(w)
        else:
            words[w] += 1
        nb_kmers += 1

    # merge reverse complement frequencies
    for k, v in words.items():
        frequencies[mapper[k]] += v

    # calculate frequencies of all words for this sequence
    for k, v in frequencies.items():
        frequencies[k] = v / float(nb_kmers + len(words_without_complement))
        # frequencies[k] = v - 1 # if count, do v - 1

    # put frequencies in the matrix
    return [frequencies[word] for word in words_without_complement]


def calculate_frequencies_matrix(fasta_file, alphabet, l, minlen, n_threads, nb_max_sequences, verbose, filtered_fasta):
    dictionnary = get_product(alphabet, l)
    words = defaultdict(int)
    # initialize counts to 1
    for word in dictionnary:
        words[word] = 0
    # get alphabet without the reverse complements and sort
    words_without_complement, mapper = get_words_without_complement(words.keys())
    words_without_complement = sorted(words_without_complement)

    # read fasta sequence by sequence
    t1, t2 = itertools.tee(read_fasta(fasta_file, minlen, filtered_fasta))
    # get total of sequences
    total = len(list(t2))
    if nb_max_sequences > 0:
        total = nb_max_sequences
    print("[x] %s sequences kept." % total)
    # initialize matrix
    matrix = np.zeros((total, len(words_without_complement)), dtype=float)
    # initialize frequencies to 1
    freq = dict([(w, 1) for w in words_without_complement])
    # launch calculations in a pool
    pool = multiprocessing.Pool(n_threads)
    manager = Manager()
    # store the excluded kmers
    excluded_kmers = manager.list()
    results = []
    for idx, sequence in enumerate(t1):
        if idx >= total:
            break
        if verbose:
            update_progress((idx + 1) / float(total), 'Preparing')
        r = pool.apply_async(row_frequency, [
            idx, freq.copy(), mapper, sequence,
            l, words, words_without_complement, excluded_kmers])
        results.append(r)

    for idx, result in enumerate(results):
        if verbose:
            update_progress((idx + 1) / float(total), 'Waiting to finish')
        matrix[idx] = result.get()

    sys.stdout.write("\n")
    return words_without_complement, matrix, excluded_kmers


def apply_geomean(matrix, verbose):
    dim, l = matrix.shape
    l = float(l)

    def apply_transform(array):
        geomean = 1
        for a in array:
            if verbose:
                pass
                # update_progress((idx + 1) / float(total), 'Applying')
            geomean *= math.pow(a, 1 / l)
        return [math.log(a / geomean) for a in array]

    return np.apply_along_axis(apply_transform, 1, matrix)


def apply_PCA(matrix, nDim):
    pca = PCA(n_components=nDim)
    return pca, pca.fit(matrix).transform(matrix)


def make_matrix(fasta_file, alphabet, kmer_len, nb_dimension_reduction, minlen, n_threads=1, exclude_complement=True,
                nb_max_sequences=-1, verbose=True, filtered_fasta=None):
    if not isinstance(nb_dimension_reduction, (list, tuple)):
        nb_dimension_reduction = [nb_dimension_reduction]
    print("[x] Calculating frequencies matrix")
    headers, matrix, excluded_kmers = calculate_frequencies_matrix(
        fasta_file, alphabet, kmer_len, minlen, n_threads, nb_max_sequences, verbose, filtered_fasta
    )
    print("[x]", len(excluded_kmers), "excluded kmers.")
    print('[x] Applying geometric mean')
    matrix = apply_geomean(matrix, verbose)
    results = []
    for nDim in nb_dimension_reduction:
        print("[x] Running PCA")
        if(matrix.shape[1] <= nDim):
            print("[x] Out-dimensionality sames as original -> Returning unprocessed original matrix: %s" % (matrix.shape[1]))
        else:
            print("[x] Reducing dimensions %s -> %s." % (matrix.shape[1], nDim))
            info, newmatrix = apply_PCA(matrix, nDim)
            print('[x] Explained amount of variance - components:', np.sum(info.explained_variance_ratio_))
            results.append((nDim, newmatrix))
    return matrix, results


if __name__ == "__main__":
    arguments = docopt(__doc__, version='VizBin wrapper 1.0')
    POOL = int(arguments['-s'])
    LEN = int(arguments['-l'])
    FATASFILE = arguments['INPUT']
    OUTPUT = arguments['OUTPUT']
    DIMS = [int(a) for a in arguments['DIMENSIONS']]
    THETA = arguments['-t']
    PERP = arguments['-p']
    MINLEN = int(arguments['-f'])
    FILTERED_FASTA = arguments['-o']
    ALPHABET = 'ATGC'
    TOTAL_SEQUENCES = -1

    original_matrix, results = make_matrix(
        FATASFILE, ALPHABET, LEN, DIMS, MINLEN,
        n_threads=8, exclude_complement=True, nb_max_sequences=-1, verbose=True, filtered_fasta=FILTERED_FASTA)
    for dim, matrix in results:
        with open('%s.%s.dat' % (OUTPUT, dim), 'w') as whandle:
            whandle.write('%s\n' % matrix.shape[0])
            whandle.write('%s\n' % dim)
            whandle.write('%s\n' % THETA)
            whandle.write('%s\n' % PERP)
            size = matrix.shape[0]
            print('DATA file: (%s.%s.dat)' % (OUTPUT, dim))
            for idx, line in enumerate(matrix):
                update_progress((idx + 1) / float(size), 'Writing to file')
                whandle.write('%s\n' % '\t'.join([str(s) for s in line]))

    # for item in row:
    #     print item,
