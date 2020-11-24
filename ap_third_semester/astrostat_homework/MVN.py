import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable


class MultivariateNormal():
    """
    A Multivariate Normal Distribution.
    """

    number_points = 200
    number_sigmas = 3

    def __init__(self, mean, cov):

        self.mean = np.array(mean)
        self.cov = np.array(cov)
        self.dim = self.mean.shape[0]

        self.normalization = (2 * np.pi)**(-self.dim / 2) / \
            np.sqrt(np.linalg.det(self.cov))
        self.precision_matrix = np.linalg.inv(self.cov)

    def __str__(self):
        return (f'{self.dim}-dimensional MVN with mean {self.mean} and covariance {self.cov}')

    def __repr__(self):
        return (f'MultivariateNormal({self.mean}, {self.cov})')


    def pdf(self, x):
        shifted_arg = x - self.mean
        argument = np.sum(
            (shifted_arg) @ self.precision_matrix * (shifted_arg), axis=-1)

        return self.normalization * np.exp(- 1 / 2 * argument)

    def marginalize(self, index):
        # not very general! we are only able to marginalize down
        # to one parameter right now
        # may generalize later
        
        marginal_mean = self.mean[index, np.newaxis]
        marginal_cov = self.cov[index, index, np.newaxis, np.newaxis]
        return self.__class__(marginal_mean, marginal_cov)

    def condition(self, index, condition_other_params):
        precision_matrix = np.linalg.inv(self.cov)
        conditioned_cov = np.linalg.inv(
            precision_matrix[index, index, np.newaxis, np.newaxis])

        cross_precision = np.delete(
            precision_matrix, index, axis=-1)[index, np.newaxis, :]

        mean_other_params = np.delete(self.mean, index, -1)
        mean_correction = conditioned_cov @ cross_precision @ (
            condition_other_params - mean_other_params)

        conditioned_mean = self.mean[index, np.newaxis] - mean_correction

        return self.__class__(conditioned_mean, conditioned_cov)

    @staticmethod
    def _analytical_CI(mean, std, percentage):
        number_sigmas = norm.ppf((percentage + 1)/2)
        delta = number_sigmas * std
        return (mean - delta, mean + delta)
    
    def analytical_CI(self, percentage):
        
        if self.dim != 1:
            print('This only applies to univariate normals!')
        
        mean = self.mean.flatten()[0]
        std = np.sqrt(self.cov.flatten()[0])
        return self._analytical_CI(mean, std, percentage)
    
    @property
    def coordinate_arrays(self):
        number_points = type(self).number_points
        number_sigmas = type(self).number_sigmas
        
        sigmas = np.sqrt(np.diagonal(self.cov))
        starts = self.mean - number_sigmas * sigmas
        ends = self.mean + number_sigmas * sigmas
        return(np.linspace(starts, ends, num=number_points).T)
    
    @property
    def cholesky_L(self):
        return (np.linalg.cholesky(self.cov))
        
    def cholesky_sample(self, number_samples):
        independent_deviates = norm.rvs(size=(self.dim, number_samples))
        correlated_deviates = (self.cholesky_L @ independent_deviates).T
        return correlated_deviates + self.mean[np.newaxis, :]

    def plot_2d_analytical(self, chosen_x, chosen_y, percentage):

        if self.dim != 2:
            print('This only applies to bivariate normals!')
            return None

        marginal_x = self.marginalize(0)
        marginal_y = self.marginalize(1)

        conditioned_x = self.condition(0, chosen_y)
        conditioned_y = self.condition(1, chosen_x)

        x_array, y_array = self.coordinate_arrays
        x, y = np.meshgrid(x_array, y_array)
        positions = np.stack((x, y), axis=-1)
        z = self.pdf(positions)
        
        mx_interval = marginal_x.analytical_CI(percentage)
        my_interval = marginal_y.analytical_CI(percentage)
        cx_interval = conditioned_x.analytical_CI(percentage)
        cy_interval = conditioned_y.analytical_CI(percentage)
        
        fig = plt.figure(figsize=(10, 10))
        gs = gridspec.GridSpec(
            2, 2, width_ratios=[2, 1], height_ratios=[1, 2])

        ax1 = plt.subplot(gs[0])
        line_mx, = ax1.plot(x_array, marginal_x.pdf(
            x_array[:, np.newaxis]), label='Marginal')
        line_cx, = ax1.plot(x_array, conditioned_x.pdf(
            x_array[:, np.newaxis]), label='Conditioned')
        [ax1.axvline(x=mx_value, c=line_mx._color) for mx_value in mx_interval]
        [ax1.axvline(x=cx_value, c=line_cx._color) for cx_value in cx_interval]
        ax1.legend()

        ax4 = plt.subplot(gs[3])
        line_my, = ax4.plot(marginal_y.pdf(
            y_array[:, np.newaxis]), y_array, label='Marginal')
        line_cy, = ax4.plot(conditioned_y.pdf(
            y_array[:, np.newaxis]), y_array, label='Conditioned')
        [ax4.axhline(y=my_value, c=line_my._color) for my_value in my_interval]
        [ax4.axhline(y=cy_value, c=line_cy._color) for cy_value in cy_interval]
        ax4.legend()

        ax3 = plt.subplot(gs[2], sharex=ax1, sharey=ax4)
        colors = ax3.contourf(x, y, z)
        ax3.axvline(chosen_x, c='white')
        ax3.axhline(chosen_y, c='white')

        ax2 = plt.subplot(gs[1])
        ax2.set_visible(False)
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes('left', size='20%', pad=0.05)
        cbar = fig.colorbar(colors, cax=cax)
        cbar.ax.set_ylabel('Probability density')
