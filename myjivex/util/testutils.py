import matplotlib.pyplot as plt
import numpy as np


def showmat(A, mask_zeros=False, title=None, **kwargs):
    if hasattr(A, "toarray"):
        A = A.toarray()

    if mask_zeros:
        A = np.ma.masked_where(np.isclose(A, 0), A)

    plt.figure()
    plt.imshow(A, **kwargs)

    if title is not None:
        plt.title(title)
    plt.show()
