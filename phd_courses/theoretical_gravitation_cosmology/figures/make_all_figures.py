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
    plt.savefig(str(plotting_func.__name__).split(sep='.')[0] + '.pdf', bbox_inches='tight', pad_inches = 0)
    plt.close()


if __name__ == "__main__":
    plotter_list = []
    for plotting_func in tqdm(plotter_list):
        plot_and_save(plotting_func)