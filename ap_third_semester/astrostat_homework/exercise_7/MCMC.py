import numpy as np
import matplotlib.pyplot as plt
from time import time
from tqdm import tqdm
from KDEpy.FFTKDE import FFTKDE
import seaborn as sns


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
        self.number_steps = 0
        self.chain = None
        self.calculate_chain(number_steps, **kwargs)
    

    def calculate_chain(self, number_steps, **kwargs):
        
        self.number_steps += number_steps

        t1 = time()
        theta = self.chain[-1] if self.chain is not None else self.initial_position
        theta_arr = [theta]
        
        cycler = range(number_steps - 1)
        
        progress_bar = kwargs.pop('progress_bar', None)
        if progress_bar:
            cycler = tqdm(cycler)
        
        for _ in cycler:
            theta = self.next_chain_step(theta)
            theta_arr.append(theta)

        if self.chain is not None:
            self.chain = np.append(self.chain, theta_arr, axis=0)
        else:
            self.chain = np.array(theta_arr)

        t2 = time()
        verbose = kwargs.pop('verbose', False)
        if verbose:
            print(f'Calculated {number_steps} steps in {t2-t1:.2f} seconds.')
            print(f'This means {number_steps / (t2-t1):.0f} steps per second.')

    def trim_chain(self, trim_number):
        if trim_number > self.trim:
            self.chain = self.chain[trim_number-self.trim:, ...]
            self.trim += trim_number
    
    @property
    def effective_steps(self):
        return self.number_steps - self.trim 

    @property
    def mean(self):
        return (np.average(self.chain, axis=0))

    def autocorrelation(self, tau):
        deviation = self.chain - self.mean[np.newaxis,:]
        
        numerator = np.sum(deviation[: - tau, ...] * deviation[tau:, ...], axis=0)
        denominator = np.sum(deviation[: - tau, ...]** 2, axis=0)
        
        # This expression was wrong in the slides! 
        normalization = self.number_steps / (self.number_steps - tau)

        # We return an array, which is commonly denoted as rho in formulas
        # This is the autocorrelation corresponding to the two coordinate axes,
        # if one is computing a certain function of the parameters then a 
        # custom correlation should be used
        return(normalization * numerator / denominator)

    def autocorrelation_array(self, max_tau=None):
        if not max_tau:
            max_tau = self.number_steps // 10
        taus = range(1, max_tau + 1)
        autocorrelations = np.array([self.autocorrelation(tau) for tau in taus])
        return (taus, autocorrelations)

    def autocorrelation_time(self, *args, C=4):
        taus, autocorrelations = self.autocorrelation_array(*args)
        times = 1 + 2 * np.cumsum(autocorrelations, axis=0)

        # Implementing an idea from
        # https://dfm.io/posts/autocorr/
        # Basically, since the autocorrelations become noise-dominated fast,
        # we compare a linear function to the estimated correlation, 
        # which should taper off, to find out where to stop integrating.
        
        finished = False
        best_max_tau = taus[-1] - 1
        for tau, time in zip(taus, times):
            if finished:
                break
            for t in time:
                if tau > C * t:
                    best_max_tau = tau
                    finished = True
        
        return (np.average(times[best_max_tau, :]))
        
    @property
    def covariance(self):
        deviation = self.chain - self.mean[np.newaxis, ...]
        cov = np.sum(deviation[:, :, np.newaxis] *
                     deviation[:, np.newaxis, :], axis=0)
        return (cov / (self.number_steps - 1))

    def autocorrelation_plot(self, *args, **kwargs):
        taus, autocorrelations = self.autocorrelation_array(*args)
        ax = kwargs.pop('ax', None)
        if not ax:
            ax = plt.gca()        
        ax.plot(taus, autocorrelations, label='Autocorrelation', **kwargs)
        ax.set_xlabel('$\\tau$')

    def posterior_plot(self):
        posterior_arr = self.posterior(self.chain)
        plt.plot(posterior_arr)
        plt.xlabel('Step number')
        plt.ylabel('Posterior evaluated at the chain step')

    def steps_trace(self, every):
        posterior_arr = self.posterior(self.chain)
        log_posterior = -np.log(posterior_arr)
                
        plotted_points = self.effective_steps // every
        trace = np.zeros(plotted_points)
        for i in range(plotted_points):
            trace[i] = np.average(log_posterior[i * every:(i + 1) * every], axis=0)
        steps = np.arange(len(trace)) * every
        return(steps, trace)

    def trace_plot(self, **kwargs):
        every = kwargs.pop('every', 1)

        steps, trace = self.steps_trace(every)
        # log_posterior = -np.log(posterior_arr)
        # trace = np.cumsum(log_posterior) / np.arange(1,
                                                    #  1 + self.number_steps - self.trim)
        initial_str = ', '.join([f'{i:.1f}' for i in self.initial_position])
        plt.plot(steps, trace, label=initial_str, **kwargs)

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
    name = 'Metropolis-Hastings'

    def __init__(self, proposal, *args, **kwargs):
        self.proposal = proposal
        self.rejections = 0
        super().__init__(*args, **kwargs)
    
    def __str__(self):
        return('Metropolis-Hastings')

    def next_chain_step(self, theta):
        new_theta = theta + self.proposal(theta)
        post_ratio = self.posterior(new_theta) / self.posterior(theta)
        
        if(np.random.uniform(low=0., high=1.) < post_ratio):
            # if post_ratio > 1  this is automatically satisfied
            # in this case, we are accepting the proposal
            return (new_theta)
            
        # otherwise, we are rejecting it, therefore we must count it again
        self.rejections+=1
        return (theta)

    @property
    def rejection_rate(self):
        return self.rejections / self.number_steps
    
    @property
    def acceptance_rate(self):
        return (1 - self.rejection_rate)
    

class Gibbs(Sampler):
    """
    Gibbs sampling: through the conditioned probability density.
    """
    name = 'Gibbs'

    def __init__(self, conditional, *args, **kwargs):
        self.conditional = conditional
        super().__init__(*args, **kwargs)
    
    def __str__(self):
        return('Gibbs')

    def next_chain_step(self, theta):

        new_theta = np.copy(theta)

        for i, t in enumerate(theta):
            new_theta[i] = self.conditional(i, new_theta)
            
        return (new_theta)
        
class Cholesky(Sampler):
    """
    Cholesky sampling: not a MCMC, but included in this framework for convenience.
    """
    name = 'Cholesky'

    def __init__(self, cholesky_sample, *args, **kwargs):
        self.cholesky_sample = cholesky_sample
        super().__init__(*args, **kwargs)
    
    def __str__(self):
        return('Cholesky')

    def next_chain_step(self, theta):

        return(self.cholesky_sample(1).reshape(theta.shape))

class SampleSet2D():
    
    
    def __init__(self, samples, *args, **kwargs):
        self.samples = np.array(samples)
        if self.samples.shape[-1] != 2:
            raise (NotImplementedError('We only support bidimensional plots for now'))

    @property            
    def mean(self):
        return (np.average(self.samples, axis=0))
    
    @property
    def covariance(self):
        deviation = self.samples - self.mean[np.newaxis, :]
        sum_cov = np.sum(deviation[:, :, np.newaxis] *
                     deviation[:, np.newaxis, :], axis=0)
        return (sum_cov / (self.samples.shape[0] - 1))
        
    def samples_plot(self, CL=None, **kwargs):
        grid = sns.jointplot(x=self.samples[:, 0], y=self.samples[:, 1], **kwargs)

        if CL:
            interval_x, interval_y = self.intervals(CL)

            [grid.ax_marg_x.axvline(x) for x in interval_x]
            [grid.ax_marg_y.axhline(y) for y in interval_y]
    
    def samples_plot_kde(self, number_sigmas=3, **kwargs):
        x, y, z = self.kde()
        
        mean = self.mean
        sigmas = np.sqrt(np.diag(self.covariance))
        limits_lower = mean - number_sigmas * sigmas
        limits_upper = mean + number_sigmas * sigmas
        
        plt.contourf(x, y, z, **kwargs)
        plt.xlim(limits_lower[0], limits_upper[0])
        plt.ylim(limits_lower[1], limits_upper[1])

    def intervals(self, CL):
        interval_x = Sampler.interval_from_samples(self.samples[:, 0], CL)
        interval_y = Sampler.interval_from_samples(self.samples[:, 1], CL)
        return (interval_x, interval_y)
    
    def marginal(self, index):
        return (self.samples[:, index])
    
    def kde(self, num = 2**10):
        nn = (num, num)
        x, y = FFTKDE().fit(self.samples).evaluate((num, num))
        return(x[:,0].reshape(nn), x[:,1].reshape(nn), y.reshape(nn))
    
    def conditional_kde(self, index, other_parameter):
        x, y, z = self.kde()
        dx = x[1, 0] - x[0, 0]
        dy = y[0, 1] - y[0, 0]
        
        if index == 0:
            cond_index = np.isclose(other_parameter, y, atol=.75 * dy)
            return(x[cond_index], z[cond_index])
        elif index == 1:
            cond_index = np.isclose(other_parameter, x, atol=.75 * dx)
            return(y[cond_index], z[cond_index])

    def conditional_cut(self, index, other_parameter, thr=.1):
        s = self.samples
        if index == 0:
            cond_index = np.isclose(other_parameter, s[:,1], atol=thr)
        elif index == 1:
            cond_index = np.isclose(other_parameter, s[:,0], atol=thr)

        return(s[cond_index])
