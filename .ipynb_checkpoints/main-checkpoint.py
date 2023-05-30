import argparse
import numpy as np
import xarray as xr
import nest

import src.meanfieldcircuit as meanfield
import src.setParams as setParams

from frites.utils import parallel_func

from tqdm import tqdm

###############################################################################
# Argument parsing
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument("PROTOCOL",
                    help="Which protocol to simulate",
                    type=int)
args = parser.parse_args()

protocol = args.PROTOCOL

assert protocol in [0, 1]


def simulate_trials(gba="strong-gba", ntrials=100, simtime=1000., dt=.1,
                    sigma=.1, onsets=[], offsets=[], n_jobs=1, fixation=False,
                    lr_delays=False, verbose=False):
    params = setParams.get_params_rate_model(gba=gba)
    
    def _trial(t):
        r_s = meanfield.simulate(simtime = simtime, dt = dt, tON=onsets,
                                 tOFF=offsets, params=params, fixation=fixation,
                                 max_cond = True, seed = t + 1,
                                 sigma=sigma, lr_delays=lr_delays)
        
        return r_s
    
    parallel, p_fun = parallel_func(
        _trial, n_jobs=n_jobs, verbose=verbose,
        total=ntrials)

    rates = parallel(p_fun(t) for t in range(ntrials))
    
    rates = xr.concat(rates, "trials")
    
    return rates


if protocol == 0:
    simtime = 7000.
    onsets = [1000., 3500.]
    offsets = [1500., 3800.]
elif protocol == 1:
    simtime = 4000.
    onsets = [1500.]
    offsets = [2000.]
elif protocol == 2:
    simtime = 7000.
    onsets = [4500.]
    offsets = [5000.]
sigma = 1.

rates_w_f = simulate_trials(gba="weak-gba", ntrials=1000, simtime=simtime, dt=.1,
                            sigma=sigma, onsets=onsets, offsets=offsets,
                            n_jobs=30, fixation=True, verbose=False)

rates_w_t = simulate_trials(gba="weak-gba", ntrials=1000, simtime=simtime, dt=.1,
                            sigma=sigma, onsets=onsets, offsets=offsets,
                            n_jobs=30, fixation=False, verbose=False)

rates_s_f = simulate_trials(gba="strong-gba", ntrials=1000, simtime=simtime, dt=.1,
                            sigma=sigma, onsets=onsets, offsets=offsets,
                        n_jobs=30, fixation=True, verbose=False)

rates_s_t = simulate_trials(gba="strong-gba", ntrials=1000, simtime=simtime, dt=.1,
                            sigma=sigma, onsets=onsets, offsets=offsets,
                        n_jobs=30, fixation=False, verbose=False)

trials = np.array([0] * 1000 + [1] * 1000)

rates_w = xr.concat([rates_w_f, rates_w_t], "trials")
rates_s = xr.concat([rates_s_f, rates_s_t], "trials")

rates = xr.concat([rates_w, rates_s], "gba")
rates = rates.assign_coords({"gba": ["weak", "strong"]})
rates = rates.assign_coords({"trials": trials})

rates.to_netcdf(f"data/protocol{protocol}.nc")
