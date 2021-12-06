from tqdm import tqdm # type: ignore
import matplotlib.pyplot as plt # type: ignore 
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('text.latex', preamble=r'''\usepackage{amsmath}
\usepackage{physics}
\usepackage{siunitx}
''')
rc('figure', dpi=150)
from typing import Callable


def plot_and_save(plotting_func : Callable):
    plotting_func()
    plt.savefig(str(plotting_func.__name__).split(sep='.')[0] + '.pdf', bbox_inches='tight', pad_inches = 0)
    plt.close()

if __name__ == "__main__":
    plotter_list : list[Callable] = []
    for plotting_func in tqdm(plotter_list):
        plot_and_save(plotting_func)