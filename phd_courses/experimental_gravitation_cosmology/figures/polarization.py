import numpy as np
import matplotlib.pyplot as plt

def vector_field(x, y, h_p=1, h_c=0):
    
    # x_new = x * (1 + h_p * expon) + y * h_c * expon
    # y_new = y * (1 - h_p * expon) + x * h_c * expon
    
    matrix = np.array([
        [h_p, h_c],
        [h_c, - h_p]
    ])
    
    return np.einsum('ij, jkl->ikl', matrix, np.stack((x, y), axis=0))

def polarization():

    xs = np.linspace(-1, 1)
    ys = np.linspace(-1, 1)
    
    X, Y = np.meshgrid(xs, ys)
    
    plt.streamplot(X, Y, *vector_field(X, Y), density=1.4, color=np.linalg.norm(vector_field(X, Y), axis=0), cmap = plt.get_cmap('inferno'))
    

if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(polarization)