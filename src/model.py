import numpy as np
from .utils import phi, create_square_kernel, gaussian
from numba import *


# Define the derivative function
@vectorize(
    [
        float32(float32, float32, float32, float32, float32),
        float64(float64, float64, float64, float64, float64),
    ]
)
@njit
def dRdt(t, R, Iext, Istim, Isyn):
    """
    Computes the derivative of the neural activity R with respect to time.

    Parameters
    ----------
    t : int or float
        Time.
    R : int, float, or ndarray
        Neural activity.
    Iext : int or float
        External input.
    Istim : int or float
        Stimulus input.
    Isyn : int, float, or ndarray
        Synaptic input.

    Returns
    -------
    int, float, or ndarray
        The derivative of R with respect to time.
    """
    return -R + phi(Iext + Istim + Isyn)


def set_one_ring_connectivity(theta, J0, J1):
    """
    Sets the connectivity matrix J for a neural network with one ring of neurons.

    Parameters
    ----------
    theta : ndarray
        The angular positions of the neurons in the ring, in radians.
    J0 : int or float
        The base strength of the connections between neurons.
    J1 : int or float
        The strength modulation of the connections due to the difference in angular positions.

    Returns
    -------
    ndarray
        The connectivity matrix J, where J[i,j] is the strength of the connection from neuron i to neuron j.
    """
    J = (J0 + J1 * np.cos(-np.subtract.outer(theta, theta))) / len(theta)
    np.fill_diagonal(J, 0)
    return J


def set_external_ring_coupling(N: int, Amp: float, sigma: int) -> np.ndarray:
    """
    Create an external ring coupling matrix of size N x N, where each column
    is shifted by a distance proportional to its distance from the center of the matrix.

    Args:
    - N: int, the size of the coupling matrix in pixels
    - sigma: int, the size of the flat center of the coupling matrix in pixels

    Returns:
    - J: numpy.ndarray, an N x N external ring coupling matrix
    """
    # Initial kernel center
    kc = int(N // 2)
    # Create kernel
    kernel = create_square_kernel(N, sigma) * Amp
    # Initialize coupling matrix
    J = np.empty((N, N))
    # Create external ring coupling matrix
    for k in range(N):
        J[:, k] = np.roll(kernel, k - kc)
    return J


def create_ring(
    state="SU",
    N=100,
    trials=1,
    steps=1000,
    s_pos=0,
    s_on=0,
    s_off=0,
    A_stim=0,
    sig_stim=1,
):  
    """
    Generates a ring network of N nodes with specified connectivity and external driving current, and 
    applies external stimulus to the network.

    Parameters
    ----------
    state : str, optional
        Specifies the type of network to generate. Possible values are "SU" for synchronous state and 
        "SB" for splay state. (default: "SU")
    N : int, optional
        Number of nodes in the ring network. (default: 100)
    trials : int, optional
        Number of trials for which to generate the network. (default: 1)
    steps : int, optional
        Number of time steps for which to simulate the network. (default: 1000)
    s_pos : int, optional
        Position of the stimulus. (default: 0)
    s_on : ndarray, optional
        Array of stimulus onset times. (default: np.zeros(0))
    s_off : ndarray, optional
        Array of stimulus offset times. (default: np.zeros(0))
    A_stim : float, optional
        Amplitude of the stimulus. (default: 0)
    sig_stim : float, optional
        Standard deviation of the stimulus. (default: 1)

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        J : NxN array of connectivity strength.
        Iext : TxNxN array of external driving current.
        Istim : TxNxN array of external stimulus current.
        theta : Array of node angular positions.
        R : TxNxN array of network state.
    """
    assert state in ["SU", "SB"]

    # Get coupling values for state
    if state == "SU":
        J0, J1 = -30, -8
    else:
        J0, J1 = -25, 11

    # Nodes angular position
    theta = 2 * np.pi / N * np.arange(0, N, 1, dtype=int)
    # Connectivity strength
    J = set_one_ring_connectivity(theta, J0, J1)

    # Compute external driving current
    Iext = 0.1 * (1 - J.sum(1))
    Iext = Iext * (1 + np.random.uniform(-0.5, 0.5, size=(trials, steps, N)))
    Iext[Iext < 0] = 0

    # External stimulus
    Sstim = np.zeros((steps, N))
    Istim = np.zeros((steps, N))
    pos = 0
    for onset, offset in zip(s_on, s_off):
        g = gaussian(N=N, amp=A_stim, sig=sig_stim)
        Sstim[onset:offset, :] = True
        Istim[onset:offset, :] = (
            Sstim[onset:offset, :] * np.roll(g, s_pos[pos])[None, :]
        )
        pos = pos + 1

    R = np.empty((trials, steps, N))

    return J, Iext, Istim, theta, R
