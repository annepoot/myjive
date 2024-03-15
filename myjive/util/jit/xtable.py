from numba import njit

##########################
# numba helper functions #
##########################


@njit
def set_block_jit(data, irows, jcols, block):
    for i, irow in enumerate(irows):
        for j, jcol in enumerate(jcols):
            data[irow, jcol] = block[i, j]


@njit
def add_block_jit(data, irows, jcols, block):
    for i, irow in enumerate(irows):
        for j, jcol in enumerate(jcols):
            data[irow, jcol] += block[i, j]
