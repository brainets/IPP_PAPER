# import NEST & NEST rasterplot
import os
import nest
import nest.raster_plot
import xarray as xr
import sys
import scipy.io
import matplotlib.pyplot as plt
import numpy             as np

########################################################################################
# Reading connectivity data
########################################################################################
data = np.load(os.path.expanduser('~/storage1/projects/IPP_WANG/src/interareal/markov2014.npy'),
               allow_pickle=True).item()

# Graph parameters
Nareas     = 29             # Number of areas
Npop       = 2             # Number of local populations (one excitatory and one inhibitory)
Npop_total = Npop * Nareas # Number of total populations
# FLN matrix
#M      = (flnMat > 0).astype(int) # Binarized
flnMat = data['FLN']
M      = (flnMat > 0).astype(int) # Binarized
delayMat = data['Distances'] / 3.5
print(delayMat.max())
# Hierarchy values
h = np.squeeze(data['Hierarchy'].T)
area_names = ['V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc',
              '7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B',
              'F2','STPi','PROm','F7','8B','STPr','24c']

def simulate(simtime = 1000.0, dt = 0.2, params=None, max_cond = True,
             sigma=0., seed = 0, tON=[2000.], tOFF=[2250.], s_pos=0,
             fixation=True, lr_delays=False):

    tON = np.asarray(tON)
    tOFF = np.asarray(tOFF)
    ##########################################################################
    # Configure NEST kernel and parameters
    ##########################################################################
    lnt  = 1

    nest.ResetKernel()
    nest.SetKernelStatus({
        'resolution': dt,        # Set simulation resolution
        'print_time': False,     # Enable printing of simulation progress (-> terminal)
        'local_num_threads': lnt,  # Use two threads to build & simulate the network
        'rng_seed' : seed
    })

    nest.SetDefaults('threshold_lin_rate_ipn',  {'theta': 0.0,
                                                 'alpha': 1e6,
                                                 'mu': 0.0,
                                                 'lambda': 1.0,
                                                 'sigma': sigma,
                                                 'rectify_output': max_cond,
                                                 'linear_summation': True})

    ##########################################################################
    # Computing connectivity parameters parameters
    ##########################################################################
    tau_ex = params['tau_ex']
    tau_in = params['tau_in']

    eta    = params['eta']

    beta_e = params['beta_e']
    beta_i = params['beta_i']

    # Local weights
    local_ee = beta_e * params['weights']['wee']  * (1.0 + eta * h)
    local_ie = beta_i * params['weights']['wie']  * (1.0 + eta * h)
    local_ei = beta_e * params['weights']['wei']  * (-1)
    local_ii = beta_i * params['weights']['wii']  * (-1)
    # Long range weights
    longr_ee = beta_e * params['weights']['muee'] * (1.0 + eta * h)
    longr_ie = beta_i * params['weights']['muie'] * (1.0 + eta * h)
    # Long range connectivity matrices
    lr_ee    = np.apply_along_axis(np.multiply,0, flnMat, longr_ee)
    lr_ie    = np.apply_along_axis(np.multiply,0, flnMat, longr_ie)

    wEe      = np.diag(-1.0 + local_ee) + lr_ee
    wEi      = local_ei * np.eye(Nareas)
    wIe      = np.diag(local_ie) + lr_ie
    wIi      = (-1.0 + local_ii) * np.eye(Nareas)

    A1 = np.concatenate((wEe/tau_ex, wEi/tau_ex), axis=1) 
    A2 = np.concatenate((wIe/tau_in, wIi/tau_in), axis=1)
    A  = np.concatenate((A1, A2),   axis=0) 

    B  = A.copy()
    B[0:Nareas,:] = B[0:Nareas,:] * tau_ex
    B[Nareas:,:]  = B[Nareas:,:]  * tau_in
    BG          = np.ones([2*Nareas, 1]) * params['I']['exc']
    BG[Nareas:] = params['I']['inh']
    currs       = -np.matmul(B, BG)

    bgExc       = currs[:Nareas,0]
    bgInh       = currs[Nareas:,0]

    #########################################################################################
    # Creating local populations
    #########################################################################################
    pop = [None] * Nareas

    # Creating local populations and setting parameters
    for i in range(Nareas):
        # Creating populations 
        pop[i] = nest.Create('threshold_lin_rate_ipn', Npop)
        # Excitatory population
        nest.SetStatus(pop[i][0], {'tau': tau_ex, 'g': 1.0, 'rate': params['I']['exc']})
        # Inhibitory population
        nest.SetStatus(pop[i][1], {'tau': tau_in, 'g': 1.0, 'rate': params['I']['inh']})
        if max_cond:
            nest.SetStatus(pop[i][0], {'rectify_rate': 10.0})
            nest.SetStatus(pop[i][1], {'rectify_rate': 35.0})

    #########################################################################################
    # Connecting multimeters NEST
    #########################################################################################
    multi = nest.Create("multimeter")
    #  nest.SetStatus(multi, {"record_from":["rate"], 'withgid':True, 'interval':.2})
    nest.SetStatus(multi, {"record_from":["rate"], 'interval':.2})
    #  nest.Connect(multi, list(chain.from_iterable(pop)))
    for i in range(Nareas):
        nest.Connect(multi, pop[i][0])
        nest.Connect(multi, pop[i][1])

    #########################################################################################
    # Connecting external inputs
    #########################################################################################

    Iext_e = nest.Create('threshold_lin_rate_ipn', 1) # Excitatory external input
    Iext_i = nest.Create('threshold_lin_rate_ipn', 1) # Inhibitory external input

    nest.SetStatus(Iext_e, {'rate': 1.0, 'sigma':0.0, 'tau': 10.0, 'mu': 1.0})                                                            
    nest.SetStatus(Iext_i, {'rate': 1.0, 'sigma':0.0, 'tau': 10.0, 'mu': 1.0})                                                             

    conn = {'rule': 'one_to_one'}                       

    for n in range(Nareas):
        syn_e  = {'weight': bgExc[n], 'synapse_model': 'rate_connection_instantaneous'} 
        nest.Connect(Iext_e, pop[n][0], conn, syn_e)
        syn_i  = {'weight': bgInh[n], 'synapse_model': 'rate_connection_instantaneous'} 
        nest.Connect(Iext_i, pop[n][1], conn, syn_i)

    if not fixation:
        Iext = []
        conn = {'rule': 'one_to_one'}  
        syn_e  = {'weight': 1.0, 'delay': 0.1,  'synapse_model': 'rate_connection_delayed'}
        for ton, toff in zip(tON, tOFF):
            Iext += [nest.Create('step_rate_generator')]
            nest.SetStatus(Iext[-1], {'amplitude_times': [ton, toff], 'amplitude_values': [params['Iext'] , 0.]})
            nest.Connect(Iext[-1], pop[0][s_pos], conn, syn_e)

    #########################################################################################
    # Creating local populations connections
    #########################################################################################
    conn = {'rule': 'one_to_one'}

    for n in range(Nareas):
        # E->E
        syn_ee = {'weight': local_ee[n], 'synapse_model': 'rate_connection_instantaneous'}
        nest.Connect(pop[n][0], pop[n][0], conn, syn_ee)
        # E->I
        syn_ie = {'weight': local_ie[n], 'synapse_model': 'rate_connection_instantaneous'}
        nest.Connect(pop[n][0], pop[n][1], conn, syn_ie)
        # I->E
        syn_ei = {'weight': local_ei, 'synapse_model': 'rate_connection_instantaneous'}
        nest.Connect(pop[n][1], pop[n][0], conn, syn_ei)
        # I->I
        syn_ii = {'weight': local_ii, 'synapse_model': 'rate_connection_instantaneous'}
        nest.Connect(pop[n][1], pop[n][1], conn, syn_ii)

    #########################################################################################
    # Creating long-range populations connections
    #########################################################################################
    
    pos, pre = np.where(M == 1)
    for i,j in zip(pre,pos):
        if not lr_delays:
            syn_ee = {'weight': lr_ee[j,i], 'synapse_model': 'rate_connection_instantaneous'}
            syn_ie = {'weight': lr_ie[j,i], 'synapse_model': 'rate_connection_instantaneous'}
        else:
            d_e = nest.random.normal(mean=delayMat[j,i], std=delayMat[j,i]*0.1)
            d_i = nest.random.normal(mean=delayMat[j,i], std=delayMat[j,i]*0.1)
            syn_ee = {'weight': lr_ee[j,i], 'delay': d_e,
                      'synapse_model': 'rate_connection_delayed'}
            syn_ie = {'weight': lr_ie[j,i], 'delay': d_i,
                      'synapse_model': 'rate_connection_delayed'}
        nest.Connect(pop[i][0], pop[j][0], conn, syn_ee)
        # Long range E->I
        nest.Connect(pop[i][0], pop[j][1], conn, syn_ie)

    nest.Simulate(simtime)

    # Retrieving recorded data
    index = nest.GetStatus(multi)[0]['events']['senders']
    times = nest.GetStatus(multi)[0]['events']['times']
    rates = nest.GetStatus(multi)[0]['events']['rate']

    r     = np.zeros([Npop_total, int(rates.shape[0]/Npop_total)])

    # Getting rate for each population
    for i in range(1,Npop_total+1):
        idx    = index == i
        tidx   = times[idx]
        ridx   = rates[idx]
        r[i-1,:] = ridx

    # Finding maximum frequencies
    count    = 0
    max_freq = np.zeros(Nareas)
    for i in range(Npop_total):
        if i%2 == 0:
            idx = (tidx>1750)*(tidx<5000)
            max_freq[count] = r[i,idx].max()
            count += 1

    # Convert to data array
    r = xr.DataArray(r[::2, :], dims=("roi", "times"),
                     coords=(area_names, tidx))
    r.attrs["max_freq"] = max_freq

    return r
