import numpy as np
from numba import *


def batch(iterable, n=1):
    """
    Generate batches of `n` elements from `iterable`.

    Parameters
    ----------
    iterable : iterable
        The iterable to be batched.
    n : int, optional
        The batch size. Default is 1.

    Yields
    ------
    batch : list
        A list of `n` elements from `iterable`, except possibly the last batch
        which may have fewer than `n` elements.
    """
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx : min(ndx + n, l)]

# Define the function for phi(Iext)
@vectorize([int32(int32), float32(float32), float64(float64)])
@njit
def phi(x):
    """
    Applies the rectified linear activation function to the input.

    Parameters
    ----------
    x : int, float, or ndarray
        Input to the activation function.

    Returns
    -------
    int, float, or ndarray
        Result of applying the activation function to the input.
    """
    if x > 0:
        return x
    return 0


def gaussian(N, amp=1, sig=1.0):
    """
    Generates a 1D Gaussian distribution.

    Parameters
    ----------
    N : int
        Number of points in the distribution.
    amp : float, optional
        Amplitude (height) of the Gaussian curve. Default is 1.
    sig : float, optional
        Standard deviation of the Gaussian curve. Default is 1.0.

    Returns
    -------
    numpy.ndarray
        Array of N points corresponding to a Gaussian distribution with the specified amplitude and standard deviation.
    """
    x = np.linspace(-(N - 1) / 2.0, (N - 1) / 2.0, N)
    gauss = amp * np.exp(-0.5 * np.square(x) / np.square(sig))
    if not N % 2:
        return np.roll(gauss, 1 - int(N / 2))
    return np.roll(gauss, -int(N / 2))


def create_square_kernel(N: int, sigma: int) -> np.ndarray:
    """
    Create a 1D square kernel of size N with a flat center of size 2*sigma + 1.

    Args:
    - N: int, the size of the kernel in pixels
    - sigma: int, the size of the flat center of the kernel in pixels

    Returns:
    - kernel: numpy.ndarray, a 1D square kernel of size N with a flat center of size 2*sigma + 1
    """
    # Kernel center
    kc = int(N // 2)
    # Create kernel
    kernel = np.zeros(N)
    # Center kernel at middle point
    kernel[kc - int(sigma / 2) : kc + int(sigma / 2) + 1] = 1
    return kernel
