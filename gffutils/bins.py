# Implementation of the UCSC genome binning strategy -- heavily commented and
# with tests to help understand what's going on.
#
# Ryan Dale 2013
#
# With help from implementations in kent src and Brent Pedersen's cruzdb,
# specifically:
#
# http://genome-source.cse.ucsc.edu/gitweb/?p=kent.git\
#                                    ;a=blob;f=src/lib/binRange.c
# and
#
# https://github.com/brentp/cruzdb/blob/master/cruzdb/__init__.py
#
# For reference, Fig 7 in http://genome.cshlp.org/content/12/6/996.abstract
# looks like this:
# ----------------------------------------------------------------------------
# |                                  1                                       |
# ----------------------------------------------------------------------------
# |       2      |         3         |         4         |         5         |
# ----------------------------------------------------------------------------
# |6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 |
# ----------------------------------------------------------------------------
#
#      AAAAAAAAAAAAAAAAAAAAAAA          BBBBBBBBBBBBBBB             CC
#
# The smallest bin feature "A" fits within is 1, but it also overlaps 2-3, 7-9.
# The smallest bin feature "B" fits within is 4, but it also overlaps 1, 14-17.
# The smallest bin feature "C" fits within is 20, but it also overlaps 1, 5.


# kent src: "How much to shift to get to next larger bin."
# So each level splits by 2**3 = 8 fold (octree?)
NEXT_SHIFT = 3

# kent src: "How much to shift to get to finest bin."
# 2 ** FIRST_SHIFT is the size of the smallest bin.
FIRST_SHIFT = 17

# These offsets are the bin numbers at the beginning of each level.  In Fig 7,
# these would be 6, 2, and 1.
#
# Bins with higher numbers (e.g. 585 to 4681) are the smaller-sized bins.  Run
# print_bin_sizes() for more info on this.
OFFSETS = [
    4096 + 512 + 64 + 8 + 1,  # bins 4681-585
    512 + 64 + 8 + 1,         # bins 585-73
    64 + 8 + 1,               # bins 73-9
    8 + 1,                    # bins 9-1
    1                         # bin  0
]

# for BED (0-based, half-open) or GFF (1-based, closed intervals)
COORD_OFFSETS = {'bed': 0, 'gff': 1}

MAX_CHROM_SIZE = 2 ** 29

def bins(start, stop, fmt='gff', one=True):
    """
    Uses the definition of a "genomic bin" described in Fig 7 of
    http://genome.cshlp.org/content/12/6/996.abstract.

    Parameters
    ----------
    one : boolean
        If `one=True` (default), then only return the smallest bin that
        completely contains these coordinates (useful for assigning a single
        bin).

        If `one=False`, then return the set of *all* bins that overlap these
        coordinates (useful for looking for features that could intersect)

    fmt : 'gff' | 'bed'
        This specifies 1-based start coords (gff) or 0-based start coords (bed)
    """

    # For very large coordinates, return 1 which is "somewhere on the
    # chromosome".
    if start >= MAX_CHROM_SIZE or stop >= MAX_CHROM_SIZE:
        if one:
            return 1
        else:
            return set([1])

    # Jump to highest resolution bin that will fit these coords (depending on
    # whether we have a BED or GFF-style coordinate).
    #
    # Some GFF files include negative coords, which will throw off this
    # calculation.  If negative coords, then set the bin to the largest
    # possible.
    if start < 0:
        if one:
            return 1
        else:
            return set([1])

    if stop < 0:
        if one:
            return 1
        else:
            return set([1])

    start = (start - COORD_OFFSETS[fmt]) >> FIRST_SHIFT
    stop = (stop) >> FIRST_SHIFT

    # We always at least fit within the chrom, which is bin 1.
    bins = set([1])

    for offset in OFFSETS:
        # Since we're going from smallest to largest bins, the first one where
        # the feature's start and stop positions are both within the same bin
        # is the smallest one these coords fit within.
        if one:
            if start == stop:
                # Note that at this point, because of the bit-shifting, `start`
                # is the number of bins (at this current level).  So we need to
                # add it to `offset` to get the actual bin ID.
                return offset + start

        # See the Fig 7 reproduction above to see why range().
        bins.update(list(range(offset + start, offset + stop + 1)))

        # Move to the next level (8x larger bin size; i.e., 2**NEXT_SHIFT
        # larger bin size)
        start >>= NEXT_SHIFT
        stop >>= NEXT_SHIFT

    return bins


def print_bin_sizes():
    """
    Useful for debugging: how large is each bin, and what are the bin IDs?
    """
    for i, offset in enumerate(OFFSETS):
        binstart = offset
        try:
            binstop = OFFSETS[i + 1]
        except IndexError:
            binstop = binstart

        bin_size = 2 ** (FIRST_SHIFT + (i * NEXT_SHIFT))
        actual_size = bin_size

        # nice formatting
        bin_size, suffix = bin_size / 1024, 'Kb'
        if bin_size >= 1024:
            bin_size, suffix = bin_size / 1024, 'Mb'
        if bin_size >= 1024:
            bin_size, suffix = bin_size / 1024, 'Gb'
        size = '(%s %s)' % (bin_size, suffix)
        actual_size = '%s bp' % (actual_size)

        print('level: {i:1};  bins {binstart:<4} to {binstop:<4}; '
              'size: {actual_size:<12} {size:<6}'.format(**locals()))


def test():
    # These should obviously fit inside the first bin
    assert bins(0, 1, fmt='bed') == 4681
    assert bins(1, 1, fmt='gff') == 4681

    # All bins that overlap:
    assert bins(0, 1, fmt='bed', one=False) == set([1, 73, 9, 585, 4681])

    # Or, more generally,
    assert bins(0, 1, fmt='bed', one=False) == set(OFFSETS)

    # Since the smallest bin size is 2**17, this rather large feature should
    # fit completely inside the first bin, too
    assert bins(2 ** 17 - 1, 2 ** 17 - 1, fmt='bed') == 4681

    # or, more generally,
    assert bins(
        2 ** FIRST_SHIFT - 1,
        2 ** FIRST_SHIFT - 1,
        fmt='bed') == OFFSETS[0]

    # At exactly 2**17 in BED coords it moves over to the next bin
    assert bins(2 ** 17, 2 ** 17, fmt='bed') == 4682

    # or
    assert bins(
        2 ** FIRST_SHIFT,
        2 ** FIRST_SHIFT,
        fmt='bed') == OFFSETS[0] + 1

    # The other bins it overlaps should be similar to before, with the
    # exception of 4682
    assert bins(2 ** 17, 2 ** 17, fmt='bed', one=False) \
        == set([1, 9, 73, 585, 4682])

    assert bins(2 ** 17, 2 ** 17, fmt='bed', one=False)\
        .difference(bins(2 ** 17 - 1, 2 ** 17 - 1, fmt='bed', one=False)) \
        == set([4682])

    # Spanning the first and second smallest bins should cause it to bump up
    # a level
    assert bins(2 ** 17 - 1, 2 ** 17, fmt='bed') == 585

    # or
    assert bins(
        2 ** FIRST_SHIFT - 1,
        2 ** FIRST_SHIFT,
        fmt='bed') == OFFSETS[1]

    # Make it as big as it can get in this bin:
    assert bins(
        2 ** FIRST_SHIFT - 1,
        2 ** (FIRST_SHIFT + NEXT_SHIFT) - 1,
        fmt='bed') == OFFSETS[1]

    # Then make it big enough to jump another level
    assert bins(
        2 ** FIRST_SHIFT - 1,
        2 ** (FIRST_SHIFT + NEXT_SHIFT),
        fmt='bed') == OFFSETS[2]

    # When start or stop or both are negative, bin should be the largest
    # possible.
    assert bins(-1, 1000, one=True) == 1
    assert bins(-1, 1000, one=False) == set([1])
    assert bins(-1, -1000, one=True) == 1
    assert bins(-1, -1000, one=False) == set([1])
    assert bins(1, -1000, one=True) == 1
    assert bins(1, -1000, one=False) == set([1])

    # Juuuust fits inside the max chrom size
    assert bins(536870910, 536870911, one=True) == 8776

    # Too big; falls back to 1.
    assert bins(536870911, 536870912, one=True) == 1

if __name__ == "__main__":
    print_bin_sizes()
    test()
