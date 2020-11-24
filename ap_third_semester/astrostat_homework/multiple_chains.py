import numpy as np
from tqdm import tqdm
from MCMC import Sampler, MetropolisHastings
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt

class MultipleChains(object):
    """
    Contains several Monte Carlo Markov Chains, 
    sampling using 'sampler' (MetropolisHastings only for now)
    the posterior 'posterior'
    starting from all the 'initial_positions'
    and going on for 'number_steps' for each.
    """

    def __init__(self, sampler_class, posterior, initial_positions, number_steps, *args, **kwargs):

        self.sampler_class = sampler_class
        self.number_steps = number_steps
        self.number_chains = len(initial_positions)
        self.initial_positions = initial_positions
        self.args = args
        self.posterior = posterior
        self.calculate_chains(**kwargs)

    def calculate_chains(self, **kwargs):
        
        if 'parallel' in kwargs.keys():
            parallel = kwargs['parallel']
        else:
            parallel = True

        global _func  # need this for parallelization

        def _func(pos):
            return(self.sampler_class(*self.args, self.posterior, pos, self.number_steps))

        if parallel:

            pool = Pool(cpu_count()-1)

            self.samplers = list(pool.map(_func, self.initial_positions))

            pool.close()  # no more tasks
            pool.join()  # wrap up current tasks

        else:
            # non-parallel version, for bugfixing:
            self.samplers = list(map(_func, self.initial_positions))

    @property
    def all_chains(self):
        return np.concatenate([s.chain for s in self.samplers], axis=0)

    def trim_chains(self, trim_number):
        for s in self.samplers:
            s.trim_chain(trim_number)

    @property
    def means(self):
        return np.array([sampler.mean for sampler in self.samplers])

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
        return(np.sum(self.covariances, axis=0) / self.number_chains)

    def R_estimator(self):
        first_term = (self.number_steps - 1) / self.number_steps
        second_term = self.means_covariance / self.average_covariance

        return (first_term + second_term)

    def traces_plot(self):
        for s in self.samplers:
            s.trace_plot()
        plt.legend()


if __name__ == "__main__":
    mean_1 = np.array([4, 2])
    mean_2 = np.array([-1, 0])
    covariance = np.array([
        [1.44, -.702],
        [-.702, .81]
    ])

    from scipy.stats import multivariate_normal

    def my_MVN(x):
        a = multivariate_normal(mean=mean_1, cov=covariance).pdf(x)
        b = multivariate_normal(mean=mean_2, cov=covariance).pdf(x)
        return (a + 2 * b)

    def gaussian_proposal(theta=None):
        return (np.random.normal(scale=1, size=2))

    initial_positions = np.stack((np.arange(-10, 10), np.arange(-10, 10))).T

    trim_amount = 1000
    chain_length = 10000

    mc = MultipleChains(MetropolisHastings, my_MVN, initial_positions,
                       chain_length + trim_amount, gaussian_proposal)
    mc.trim_chains(trim_amount)

    # dim = 5

    # c = np.random.random(size=(dim, dim))
    # cov = c.T @ c

    # def new_MVN(x):
    #     return (multivariate_normal(mean= np.zeros(dim),cov=cov).pdf(x))
    # def big_gaussian_proposal(theta=None):
    #     return (np.random.normal(scale=1, size=dim))

    # big_initial_positions = np.random.normal(size=(7, dim))

    # mc = MultipleChains(MetropolisHastings, new_MVN, big_initial_positions, chain_length + trim_amount, big_gaussian_proposal)
    # mc.trim_chains(trim_amount)
