import numpy as np
from tqdm import tqdm
from MCMC import Sampler, MetropolisHastings
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
from scipy.stats import norm

class MultipleChains(object):
    """
    Contains several Monte Carlo Markov Chains*,
    sampling using 'sampler' (MetropolisHastings only for now)
    the posterior 'posterior'
    starting from all the 'initial_positions'
    and going on for 'number_steps' for each.
    
    * or Cholesky sampler, included in the same framework for simplicity
    """

    def __init__(self, sampler_class, posterior, initial_positions, number_steps, *args, **kwargs):

        self.sampler_class = sampler_class
        self.number_steps = number_steps
        self.number_chains = len(initial_positions)
        self.initial_positions = initial_positions
        self.args = args
        self.posterior = posterior
        self.calculate_chains(**kwargs)
    
    def __str__(self):
        return(f'{self.number_chains} {self.sampler_class.name} chains.')

    def calculate_chains(self, **kwargs):
        
        parallel = kwargs.pop('parallel', True)

        global _initialize_samplers  # need this for parallelization

        def _initialize_samplers(pos):
            return(self.sampler_class(*self.args, self.posterior, pos, self.number_steps, **kwargs))

        if parallel:

            pool = Pool(cpu_count()-1)

            self.samplers = list(pool.map(_initialize_samplers, self.initial_positions))

            pool.close()  # no more tasks
            pool.join()  # wrap up current tasks

        else:
            # non-parallel version, for bugfixing:
            self.samplers = list(map(_initialize_samplers, self.initial_positions))

    def extend_chains(self, number_steps, **kwargs):
        parallel = kwargs.pop('parallel', True)
        self.number_steps += number_steps

        global _extend_sampler  # need this for parallelization

        def _extend_sampler(sampler):
            sampler.calculate_chain(number_steps, **kwargs)
            return(sampler)

        if parallel:

            pool = Pool(cpu_count()-1)

            self.samplers = list(pool.map(_extend_sampler, self.samplers))

            pool.close()  # no more tasks
            pool.join()  # wrap up current tasks

        else:
            # non-parallel version, for bugfixing:
            self.samplers = list(map(_extend_sampler, self.samplers))
    

    @property
    def all_chains(self):
        return np.concatenate([s.chain for s in self.samplers], axis=0)

    def trim_chains(self, trim_number):
        if trim_number > 0:
            print(f'Trimming at {trim_number}')
            for s in self.samplers:
                s.trim_chain(trim_number)
        elif trim_number == 0:
            print('No trimming')
            
    @property
    def means(self):
        return np.array([sampler.mean for sampler in self.samplers])

    @property    
    def average_mean(self):
        return(np.average(self.means, axis=0))
        
    @property
    def means_covariance(self):
        means = self.means
        means_mean = np.average(means, axis=0)

        deviation = means - means_mean[np.newaxis, ...]
        cov = np.sum(deviation[:, :, np.newaxis] *
                     deviation[:, np.newaxis, :], axis=0)
        return (cov / (self.number_chains - 1))

    @property
    def covariances(self):
        return np.array([sampler.covariance for sampler in self.samplers])

    @property
    def average_covariance(self):
        return(np.average(self.covariances, axis=0))

    def R_estimator(self):
        first_term = (self.number_steps - 1) / self.number_steps
        second_term = self.means_covariance / self.average_covariance

        return (first_term + second_term)
        
    def autocorrelation_times(self, *args):
        return([sampler.autocorrelation_time(*args) for sampler in self.samplers])

    @property
    def optimal_trimming(self):
        over_thrs = []
        too_low = 0
        for sampler in self.samplers:
            steps, trace = sampler.steps_trace(every=1)
            N = sampler.effective_steps
            late_trace = trace[N // 2:]
            std_trace = np.std(late_trace)
            mean_trace = np.average(late_trace)
            
            # adaptive threshold, the number is arbitrary
            # but ~5 seems to be a good choice
            number_sigmas = norm.isf(1 / sampler.effective_steps) * 5

            # so that it does not depend on the normalization of the posterior
            thr = mean_trace + number_sigmas * std_trace
            
            # get index of last occurrence of "tr > thr"
            over_thr = N - next((i for i, tr in enumerate(reversed(trace)) if tr > thr), N) - 1
            
            if(over_thr > N // 2):
                print('threshold is too high!')
                over_thr = N // 2 - 1
            if(over_thr <= 0):
                too_low += 1
                over_thr = 0

            over_thrs.append(over_thr)
        
        if too_low == self.number_chains:
            print('All thresholds too low!')
        
        # another rather arbitrary number here
        # ~2 seems good

        return(2 * max(over_thrs))

    def traces_plot(self, **kwargs):
        for s in self.samplers:
            s.trace_plot(**kwargs)
        plt.legend()
        plt.xlabel('Step number')
        plt.ylabel('$-\\log $ posterior')

def errors_sampler(sampler, mvn, trimming_index=40, num=200, max_num=int(1e6)):
    errors_cov = []
    errors_mean = []

    N = np.geomspace(sampler.number_steps, max_num, dtype=int, num=num)
    diffs = np.ediff1d(N, to_end=0) # discrete gradient of N

    for i, diff in tqdm(enumerate(diffs)):
        errors_cov.append(np.abs(sampler.average_covariance - mvn.cov))
        errors_mean.append(np.abs(sampler.average_mean - mvn.mean))

        if i == trimming_index:
            sampler.trim_chains(sampler.optimal_trimming)

        sampler.extend_chains(diff)

    return N, np.array(errors_mean), np.array(errors_cov)

def plot_errors(N, errors_mean, errors_cov, multiple_chains, trimming_index=None):
    
    plt.figure(dpi=100, figsize=(8, 6))
    for i in range(2):
        for j in range(2):
            if i <= j:
                plt.loglog(N, errors_cov[:, i, j], label = f'{i}{j} component, cov')

    for i in range(2):
        plt.loglog(N, errors_mean[:, i], label=f'{i} component, mean')

    min_error = min(errors_cov.min(), errors_mean.min())
    max_error = max(errors_cov.max(), errors_mean.max())

    leg = 'Reference $N^{-1/2}$ lines'
    for C in np.geomspace(max_error * np.sqrt(N[0]), min_error * np.sqrt(N[-1]), num=10):
        plt.plot(N, N**(-1/2) * C, ls=':', c='black', alpha=.6, label=leg)
        leg=None

    if trimming_index:
        plt.axvline(N[trimming_index], ls='--', alpha = .4, label='Trimmed beginning of chains')

    plt.legend()
    plt.title('Errors for ' + str(multiple_chains))
    plt.xlabel('Number of MCMC steps')
    plt.ylabel('Absolute error')

def plot_autocorrelations(multiple_chains, max_tau=100):
    
    # get the autocorrelations at variable delays, from 1 to 100:
    t, a = multiple_chains.samplers[0].autocorrelation_array(max_tau)

    fig, ax = plt.subplots(figsize=(8,6))

    multiple_chains.samplers[0].autocorrelation_plot(max_tau, ax=ax, c='purple')
    times = 1 + 2 * np.cumsum(a, axis=0)
    ax.tick_params(axis='y', colors='purple')
    ax.grid('on')

    # duplicate y axis, in order to show both the autocorrelation and 
    # its cumulative sum on the same axis
    ax2 = ax.twinx()

    linear = np.array(t) / 4
    end_integration = min(np.argmax(linear[:, np.newaxis] > times, axis=0))

    ax2.plot(t, times, c='green', label='Running autocorrelation time')
    ax2.plot(list(t), linear, ls=':', label='Linear function', c='black')
    ax2.tick_params(axis='y', colors='green')
    ax2.axvline(end_integration, ls='-.', c='orange', label='End of integration')

    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels+labels2, handles+handles2))
    ax.legend(by_label.values(), by_label.keys())
    ax.set_title(f'{multiple_chains.samplers[0].name} sampling autocorrelation')
