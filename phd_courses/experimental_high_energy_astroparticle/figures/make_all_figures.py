from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
          \usepackage{physics}
          \usepackage{siunitx}
          ''')
rc('figure', dpi=150)


def plot_and_save(plotting_func):
    plotting_func()
    plt.savefig(str(plotting_func.__name__).split(sep='.')[0] + '.pdf')

if __name__ == "__main__":
    
    from cosmic_rays_energies import cosmic_rays_energies
    
    plotter_list = [
        cosmic_rays_energies,
        angular_distribution_shift,
    ]
    
    for plotting_func in tqdm(plotter_list):
        plot_and_save(plotting_func)