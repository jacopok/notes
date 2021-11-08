import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
          \usepackage{physics}
          \usepackage{siunitx}
          ''')
rc('figure', dpi=150)

from cosmic_rays_energies import cosmic_rays_energies

if __name__ == "__main__":
    plotter_list = [
        cosmic_rays_energies,
    ]
    
    for plotting_func in plotter_list:
        plotting_func()
        plt.savefig(str(plotting_func.__name__).split(sep='.')[0] + '.pdf')