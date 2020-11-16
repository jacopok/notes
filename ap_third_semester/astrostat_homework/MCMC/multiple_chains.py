import numpy as np
from tqdm import tqdm
from MCMC import Sampler, MetropolisHastings
from multiprocessing import Pool, cpu_count

class MultipleChain(object):
    """
    Contains several Monte Carlo Markov Chains, 
    sampling using 'sampler' (MetropolisHastings only for now)
    which 
    the posterior 'posterior'
    starting from all the 'initial_positions'
    and going on for 'number_steps' for each.
    """

    def __init__(self, sampler_class, posterior, initial_positions, number_steps, *args):

        self.sampler_class = sampler_class
        self.number_steps = number_steps
        self.number_chains = len(initial_positions)
        
        global _func # need this for parallelization
        
        def _func(pos):
            return(sampler_class(*args, posterior, pos, number_steps))
        
        pool = Pool(cpu_count()-1)
        
        self.samplers = list(pool.map(_func, initial_positions))
        
        pool.close()  # no more tasks
        pool.join()  # wrap up current tasks
        
        # non-parallel version, for bugfixing:
        # self.samplers = list(map(_func, initial_positions))
 
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
        
if __name__ == "__main__":
    mean_1 = np.array([4, 2])
    mean_2 = np.array([-1, 0])
    covariance = np.array([
        [1.44, -.702],
        [-.702, .81]
    ])

    from scipy.stats import multivariate_normal

    # def my_MVN(x):
    #     a = multivariate_normal(mean=mean_1, cov=covariance).pdf(x)
    #     b = multivariate_normal(mean=mean_2, cov=covariance).pdf(x)
    #     return (a + 2 * b)

    # def gaussian_proposal(theta=None):
    #     return (np.random.normal(scale=1, size=2))
            
    # initial_positions = np.stack((np.arange(-10, 10), np.arange(-10, 10))).T
    
    # trim_amount = 1000
    # chain_length = 10000
    
    # mc = MultipleChain(MetropolisHastings, my_MVN, initial_positions, chain_length + trim_amount, gaussian_proposal)
    # mc.trim_chains(trim_amount)
    
    dim = 2
    
    c = np.random.random(size=(dim, dim))
    cov = c.T @ c
    
    def my_MVN(x):
        return (multivariate_normal(mean= np.zeros(dim),cov=cov).pdf(x))
    def gaussian_proposal(theta=None):
        return (np.random.normal(scale=1, size=dim))
    
    initial_positions = np.random.normal(size=(7, dim))
    
    mc = MultipleChain(MetropolisHastings, my_MVN, initial_positions, chain_length + trim_amount, gaussian_proposal)
    mc.trim_chains(trim_amount)
    