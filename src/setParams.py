#########################################################################################
# Parameters to use in the rate based and spikin neuron network simulations
#########################################################################################
import numpy    as np  

#########################################################################################
# Return the parameters for the rate based network model
#########################################################################################
def get_params_rate_model(gba='weak-gba'):
    params = {}

    # Time constants
    params['tau_ex'] = 20.0 # ms
    params['tau_in'] = 10.0 # ms
    # Scale for excitatory weights
    params['eta']    = 0.68
    # Intrinsic time constants
    params['beta_e'] = 0.066
    params['beta_i'] = 0.351
    # Local/Long-range circuit connection weights
    params['weights'] = {'wee':24.3, 'wie':12.2, 'wii':12.5, 'wei':19.7, 'muie':25.3, 'muee':33.7} # in pA/Hz
    # External inputs to excitatory and inhibitory populations
    params['I']       = {'exc': 10.0, 'inh': 35.0} # in Hz
    # Input to V1
    params['Iext']    = 22.05*1.9

    if gba == 'strong-gba':
        params['weights']['wei']  = 25.2
        params['weights']['muee'] = 51.5
        params['Iext']            = 11.54*1.9

    return params
