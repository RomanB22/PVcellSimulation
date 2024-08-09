import numpy as np

def poisson_generator(rate, t_start=0.0, t_stop=1000.0, seed=None):
    """
    Returns a SpikeTrain whose spikes are a realization of a Poisson process
    with the given rate (Hz) and stopping time t_stop (milliseconds).

    Note: t_start is always 0.0, thus all realizations are as if
    they spiked at t=0.0, though this spike is not included in the SpikeList.

    Inputs:
    -------
        rate    - the rate of the discharge (in Hz)
        t_start - the beginning of the SpikeTrain (in ms)
        t_stop  - the end of the SpikeTrain (in ms)
        array   - if True, a np array of sorted spikes is returned,
                    rather than a SpikeTrain object.

    Examples:
    --------
        >> gen.poisson_generator(50, 0, 1000)
        >> gen.poisson_generator(20, 5000, 10000, array=True)

    See also:
    --------
        inh_poisson_generator, inh_gamma_generator, inh_adaptingmarkov_generator
    """

    rng = np.random.RandomState(seed)

    # number = int((t_stop-t_start)/1000.0*2.0*rate)

    # less wasteful than double length method above
    n = (t_stop - t_start) / 1000.0 * rate
    number = np.ceil(n + 3 * np.sqrt(n))
    if number < 100:
        number = min(5 + np.ceil(2 * n), 100)

    if number > 0:
        isi = rng.exponential(1.0 / rate, int(number)) * 1000.0
        if number > 1:
            spikes = np.add.accumulate(isi)
        else:
            spikes = isi
    else:
        spikes = np.array([])

    spikes += t_start
    i = np.searchsorted(spikes, t_stop)

    extra_spikes = []
    if i == len(spikes):
        # ISI buf overrun

        t_last = spikes[-1] + rng.exponential(1.0 / rate, 1)[0] * 1000.0

        while (t_last < t_stop):
            extra_spikes.append(t_last)
            t_last += rng.exponential(1.0 / rate, 1)[0] * 1000.0

        spikes = np.concatenate((spikes, extra_spikes))

    else:
        spikes = np.resize(spikes, (i,))

    return spikes


def inh_poisson_generator(rate, t, t_stop, seed=None):
    """
    Returns a SpikeTrain whose spikes are a realization of an inhomogeneous
    poisson process (dynamic rate). The implementation uses the thinning
    method, as presented in the references.

    Inputs:
    -------
        rate   - an array of the rates (Hz) where rate[i] is active on interval
                    [t[i],t[i+1]]
        t      - an array specifying the time bins (in milliseconds) at which to
                    specify the rate
        t_stop - length of time to simulate process (in ms)
        array  - if True, a np array of sorted spikes is returned,
                    rather than a SpikeList object.

    Note:
    -----
        t_start=t[0]

    References:
    -----------

    Eilif Muller, Lars Buesing, Johannes Schemmel, and Karlheinz Meier
    Spike-Frequency Adapting Neural Ensembles: Beyond Mean Adaptation and Renewal Theories
    Neural Comput. 2007 19: 2958-3010.

    Devroye, L. (1986). Non-uniform random variate generation. New York: Springer-Verlag.

    Examples:
    --------
        >> time = arange(0,1000)
        >> stgen.inh_poisson_generator(time,sin(time), 1000)

    See also:
    --------
        poisson_generator, inh_gamma_generator, inh_adaptingmarkov_generator
    """

    rng = np.random.RandomState(seed)

    if np.shape(t) != np.shape(rate):
        raise ValueError('shape mismatch: t,rate must be of the same shape')

    # get max rate and generate poisson process to be thinned
    rmax = np.max(rate)
    ps = poisson_generator(rate=rmax, t_start=t[0], t_stop=t_stop, seed=seed)

    # return empty if no spikes
    if len(ps) == 0:
        np.array([])

    # gen uniform rand on 0,1 for each spike
    rn = np.array(rng.uniform(0, 1, len(ps)))

    # instantaneous rate for each spike
    idx = np.searchsorted(t, ps) - 1
    # spike_rate = rate[idx]
    spike_rate = np.array([rate[i] for i in idx])

    # thin and return spikes
    spike_train = ps[rn < spike_rate / rmax]

    return list(spike_train)