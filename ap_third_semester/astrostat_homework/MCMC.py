import numpy as np
import matplotlib.pyplot as plt
from time import time
from tqdm import tqdm

SEED = 42
np.random.seed(SEED)


class Sampler(object):
    """
    General Markov Chain Monte Carlo sampler. 
    To be subclassed exclusively.
    """

    __slots__ = 'chain'

    def __init__(self, posterior, initial_position, number_steps=None, **kwargs):
        self.posterior = posterior
        self.initial_position = initial_position
        self.trim = 0
        if not number_steps:
            number_steps = 100
        self.calculate_chain(number_steps)

    def calculate_chain(self, number_steps):
        self.number_steps = number_steps

        t1 = time()
        theta = self.initial_position
        theta_arr = [theta]
        if number_steps > 10000:
            cycler = tqdm(range(number_steps - 1))
        else:
            cycler = range(number_steps - 1)
        for _ in cycler:
            theta = self.next_chain_step(theta)
            theta_arr.append(theta)
        self.chain = np.array(theta_arr)
        t2 = time()
        print(f'Calculated {number_steps} steps in {t2-t1:.2f} seconds.')
        print(f'This means {number_steps / (t2-t1):.0f} steps per second.')

    def trim_chain(self, trim_number):
        self.chain = self.chain[trim_number:, ...]
        self.trim += trim_number

    @property
    def mean(self):
        return (np.average(self.chain, axis=0))

    def autocorrelation(self, tau):
        deviation = self.chain - self.mean

        numerator = np.sum(deviation[: - tau, ...] * deviation[tau:, ...])
        denominator = np.sum(deviation[: - tau, ...] ** 2)
        normalization = 1 / (self.number_steps - tau)

        return(normalization * numerator / denominator)

    def autocorrelation_array(self, max_tau=None):
        if not max_tau:
            max_tau = self.number_steps // 10
        taus = range(1, max_tau + 1)
        autocorrelations = [self.autocorrelation(tau) for tau in taus]

        return (taus, autocorrelations)

    def autocorrelation_time(self, *args):
        taus, autocorrelations = self.autocorrelation_array(*args)
        return 1 + 2 * np.sum(autocorrelations)

    @property
    def covariance(self):
        deviation = self.chain - self.mean[np.newaxis, ...]
        cov = np.sum(deviation[:, :, np.newaxis] *
                     deviation[:, np.newaxis, :], axis=0)
        return (cov / (self.number_steps - 1))

    def autocorrelation_plot(self, *args):
        taus, autocorrelations = self.autocorrelation_array(*args)

        plt.plot(taus, autocorrelations)
        plt.xlabel('$\\tau$')
        plt.ylabel('Autocorrelation')

    def posterior_plot(self):
        posterior_arr = self.posterior(self.chain)
        plt.plot(posterior_arr)
        plt.xlabel('Step number')
        plt.ylabel('Posterior evaluated at the chain step')

    def trace_plot(self):
        posterior_arr = self.posterior(self.chain)
        log_posterior = -np.log(posterior_arr)
        trace = np.cumsum(log_posterior) / np.arange(1,
                                                     1 + self.number_steps - self.trim)
        initial_str = ', '.join([f'{i:.1f}' for i in self.initial_position])
        plt.plot(trace, label=initial_str)
        plt.xlabel('Step number')
        plt.ylabel('Trace')

    @staticmethod
    def interval_from_samples(samples, percentage):
        
        samples = sorted(samples)
        n_points = int(len(samples) * percentage)
        
        dist = np.inf
        for i in range(len(samples) - n_points):
            new_dist = samples[i + n_points] - samples[i]
            if new_dist < dist:
                chosen_i = i
                dist = new_dist
        
        return(samples[chosen_i], samples[chosen_i + n_points])

class MetropolisHastings(Sampler):
    """
    Metropolis - Hastings algorithm: proposal and rejection
    """

    def __init__(self, proposal, *args, **kwargs):
        self.proposal = proposal
        try:
            self.calculate_acceptance_rate = kwargs['calculate_acceptance_rate']
        except(KeyError):
            self.calculate_acceptance_rate = False
        if self.calculate_acceptance_rate:
            self.rejections = []
        super().__init__(*args, **kwargs)

    def next_chain_step(self, theta):

        if not self.calculate_acceptance_rate:
            while True:
                new_theta = theta + self.proposal(theta)
                post_ratio = self.posterior(new_theta) / self.posterior(theta)
                acceptance_prob = min(post_ratio, 1)
                if(np.random.uniform(low=0., high=1.) < acceptance_prob):
                    return (new_theta)
        else:
            n = 0
            while True:
                new_theta = theta + self.proposal(theta)
                post_ratio = self.posterior(new_theta) / self.posterior(theta)
                acceptance_prob = min(post_ratio, 1)
                if(np.random.uniform(low=0., high=1.) < acceptance_prob):
                    self.rejections.append(n)
                    return (new_theta)
                n += 1

    @property
    def rejection_rate(self):
        mean_rejections = np.average(self.rejections)
        twice_mr = 2 * mean_rejections

        return(1 + 1/twice_mr - np.sqrt(1 + 2 * twice_mr) / twice_mr)


class Gibbs(Sampler):
    """
    Gibbs sampling: through the conditioned probability density.
    """

    def __init__(self, conditional, *args, **kwargs):
        self.conditional = conditional
        super().__init__(*args, **kwargs)

    def next_chain_step(self, theta):

        new_theta = np.zeros_like(theta)

        for i, t in enumerate(theta):
            new_theta[i] = conditional(i, theta)
        return(new_theta)
    

if __name__ == "__main__":

    mean_1 = np.array([4, 2])
    mean_2 = np.array([-1, 0])
    covariance = np.array([
        [1.44, -.702],
        [-.702, .81]
    ])

    from scipy.stats import multivariate_normal

    def gaussian_proposal(theta=None):
        return (np.random.normal(scale=1, size=2))

    def my_MVN(x):
        a = multivariate_normal(mean=mean_1, cov=covariance).pdf(x)
        b = multivariate_normal(mean=mean_2, cov=covariance).pdf(x)
        return(a+2*b)

    metropolis_hastings = MetropolisHastings(
        gaussian_proposal, my_MVN,  [0, 0], 10000, calculate_acceptance_rate=True)
