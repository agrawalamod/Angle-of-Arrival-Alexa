import numpy as np
import matplotlib.pyplot as plt

def aoa(y, Fs, azi_2d, ele_2d):
    V_SOUND = 343
    PAIR_DIST = 0.045
    PAIR_ANG = np.pi*(1/3)
    USB_DEG = 90 * np.pi / 180
    ch_pairs = [[0,1],[0,2],[0,3],[0,4],[0,5]]#[[0, 3],[1, 4],[2, 5]]#
    
    MIC_LOCS = np.zeros((6,3))
    for i in range(6):
        ang_from_x = -PAIR_ANG*(i-0.5) + USB_DEG
        MIC_LOCS[i,:] = -PAIR_DIST*np.array([np.cos(ang_from_x), np.sin(ang_from_x), 0])
    
    mic_disps = np.zeros((len(ch_pairs),3))
    for i in range(len(ch_pairs)):
        mic_disps[i,:] = MIC_LOCS[ch_pairs[i][0],:] - MIC_LOCS[ch_pairs[i][1],:]
    
    taus = []
    for pair in ch_pairs:
        taus.append(gccphat(y[:,pair[0]], y[:,pair[1]],Fs))
    taus = np.array(taus)
    dists = taus*V_SOUND
    print(dists*100)
    scores = np.zeros(np.shape(azi_2d))
    print('computing heat map...')
    for A in range(np.shape(azi_2d)[1]):
        for E in range(np.shape(azi_2d)[0]):
            azi = azi_2d[E,A] * np.pi / 180
            ele = ele_2d[E,A] * np.pi / 180
            #Compute incidence vector
            vec = np.array([[np.cos(azi)*np.cos(ele)],
                            [np.sin(azi)*np.cos(ele)],
                            [np.sin(ele)]])
            scores[E,A] = 1/np.linalg.norm(np.dot(mic_disps,vec).T - dists)
    print('Done.')
    return scores


def gccphat(y1, y2, Fs):
    EPSILON = 0.000001
    assert(len(y1) == len(y2))
    N = len(y1)*2
    Nd2 = len(y1)
    FFT_PROD = np.fft.rfft(y1,n=N) * np.conj(np.fft.rfft(y2,n=N))
    x = np.fft.irfft(FFT_PROD/(np.abs(FFT_PROD) + EPSILON),n=N)
    #x = np.fft.irfft(FFT_PROD/(EPSILON),n=N)
    x = np.concatenate((x[-Nd2:], x[:Nd2+1]))
    lag = np.array(range(N+1)) - Nd2
#    plt.plot(lag, x)
#    plt.show()
    
    #Interpolate to find max
    i_corr = np.argmax(x[Nd2-8:Nd2+8]) + (Nd2-8)
    tau_approx = lag[i_corr] / Fs
    idx = i_corr #+ Nd2
    deltak = 0.5*(x[idx-1] - x[idx+1]) / (x[idx+1] + x[idx-1] - 2*x[idx])
    if abs(deltak) >= 1:
        deltak = 0
    t_offset = deltak/Fs
    tau = tau_approx + t_offset
    return tau