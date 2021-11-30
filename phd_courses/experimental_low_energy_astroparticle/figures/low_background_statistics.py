import numpy as np
import matplotlib.pyplot as plt

N=int(1e5)

def poisson():
    rng = np.random.default_rng(seed=100)
    
    t = np.linspace(0, 100)
    
    rate_sig = 1
    rate_bkg = .1
    
    sig_events = rng.poisson(lam=t * rate_sig, size=(N,) + t.shape)
    bkg_events = rng.poisson(lam=t * rate_bkg, size=(N,) + t.shape)
    
    tot_events = sig_events + bkg_events
    
    bkg_dev = np.std(bkg_events, axis=0)
    sig_avg = np.average(sig_events, axis=0)

    signal_over_avg_bkg = tot_events - t[np.newaxis, :] * rate_bkg 
    snr = np.average(signal_over_avg_bkg, axis=0) / np.std(signal_over_avg_bkg, axis=0)
    
    plt.plot(t, snr, label='SNR')
    plt.plot(t, t*rate_bkg, label='average $N_B$')
    
    plt.xlabel('Time')
    plt.title(f'Signal rate {rate_sig}, bkg rate {rate_bkg}')

    plt.legend()


if __name__ == "__main__":
    from make_all_figures import plot_and_save
    plot_and_save(poisson)