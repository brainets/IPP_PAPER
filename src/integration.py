import numpy as np

# Define the EULER solver function
def EULER(derivative, R, t, h, Iext, Istim, Isyn):
    """
    Computes the next value of the neural activity R using the Euler method.

    Parameters
    ----------
    derivative : function
        The derivative function that computes the rate of change of R with respect to time.
    R : int, float, or ndarray
        The current value of R.
    t : int or float
        The current time.
    h : int or float
        The time step.
    Iext : int or float
        External input.
    Istim : int or float
        Stimulus input.
    Isyn : int, float, or ndarray
        Synaptic input.

    Returns
    -------
    int, float, or ndarray
        The next value of R.
    """
    return R + h * derivative(t + h, R, Iext, Istim, Isyn)


# Define the RK4 solver function
def RK4(derivative, R, t, h, Iext, Istim, Isyn):
    """
    Computes the next value of the neural activity R using the fourth-order Runge-Kutta method.

    Parameters
    ----------
    derivative : function
        The derivative function that computes the rate of change of R with respect to time.
    R : int, float, or ndarray
        The current value of R.
    t : int or float
        The current time.
    h : int or float
        The time step.
    Iext : int or float
        External input.
    Istim : int or float
        Stimulus input.
    Isyn : int, float, or ndarray
        Synaptic input.

    Returns
    -------
    int, float, or ndarray
        The next value of R.
    """
    k1 = derivative(t, R, Iext, Istim, Isyn)
    k2 = derivative(t + h / 2, R + k1 / 2, Iext, Istim, Isyn)
    k3 = derivative(t + h / 2, R + k2 / 2, Iext, Istim, Isyn)
    k4 = derivative(t + h, R + k3, Iext, Istim, Isyn)
    return R + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def config_simulation(integration, tsim, h):
    """
    Configure a simulation of a dynamical system with a specified integration method.

    Args:
    - integration: str, the name of the integration method to use; either "euler" or "rk4"
    - tsim: float, the total simulation time in seconds
    - h: float, the time step size in seconds

    Returns:
    - fint: callable, the integration function corresponding to the specified integration method
    - steps: int, the number of time steps in the simulation
    - times: numpy.ndarray, an array of times corresponding to each time step in the simulation
    """
    assert integration in ["euler", "rk4"]

    # Get integration function
    fint = dict(euler=EULER, rk4=RK4)[integration]

    tinit = 0.0
    # Number of time steps
    steps = int((tsim - tinit) / h)  # Number of steps
    times = np.linspace(tinit, tsim, steps)

    return fint, steps, times
