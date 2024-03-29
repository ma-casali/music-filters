"""
Created on Sun Nov 26 20:57:33 2023
Author: Matthew Casali

"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal


mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['mathtext.fontset'] = 'cm' #LaTeX format used by overleaf

mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.grid'] = True

mpl.rcParams['lines.linewidth'] = 2

# current figsize corresponds to dpi, changing figsize requires change in dpi
mpl.rcParams['figure.figsize'] = [12,6]
mpl.rcParams['figure.autolayout'] = True # will result in fig.tight_layout()

mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['xtick.color'] = 'k'
mpl.rcParams['xtick.labelcolor'] = 'k'
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['ytick.color'] = 'k'
mpl.rcParams['ytick.labelcolor'] = 'k'

mpl.rcParams['legend.fontsize'] = 20

mpl.rcParams['savefig.dpi'] = 150 # set dpi to MATLAB default

#%% Save working wavfile

plt.close('all')
    
fs, data = wavfile.read(r"your-file.wav")
dt = 1/fs

n = 30*fs
start = 0*fs
wavfile.write(r"your-file_trimmed.wav", fs, data[start:start + n,0])

#%%

plt.close('all')
    
fs, data = wavfile.read(r"your-file_trimmed.wav")
dt = 1/fs
t_old = np.arange(len(data))*dt

df = .1 # Hz
T = 1/(2*df) # s
N = int(T/dt)

f, t_filter, Zxx = signal.stft(data, fs = fs, nperseg = N, noverlap = 3*N/4)#, return_onesided = False)
Z_new = np.empty_like(Zxx)
window_length = 3

win_wd = int((window_length - 1)/2)
freq = 1/(2*window_length*np.diff(t_filter)[0])

mult_arr = signal.windows.hamming(window_length)
ii_stop = 0

for ii in range(t_filter.shape[0]):
    
    if ii - win_wd < 0:
        mult_start = win_wd - ii
        mult_end = window_length + 1
        start = 0
        end = ii + win_wd + 1

    elif ii + win_wd + 1 > t_filter.shape[0]:
        mult_start = 0 
        mult_end = win_wd + (t_filter.shape[0] - ii)
        start = ii - win_wd
        end = t_filter.shape[0]
        
    else:
        mult_start = 0
        mult_end = window_length + 1
        start = ii - win_wd
        end = ii + win_wd + 1

    Zmultiplier = Zxx[:,start:end]*np.exp(1j*2*np.pi*np.matmul(f[:,np.newaxis], np.transpose(t_filter[ii] - t_filter[start:end,np.newaxis])))
    mult = Zmultiplier * mult_arr[mult_start:mult_end]
    div = np.sum(np.abs(mult), axis = 0)
    mult = mult / div

    Z_new[:,ii] = np.prod(mult, axis = 1)
    Z_new[:,ii]  = Z_new[:,ii] / np.sum(np.abs(Z_new[:,ii]))
    
t_Z, new_data = np.real(signal.istft(Z_new, fs = fs, nperseg = N, noverlap = 3*N/4))#, input_onesided = False))

    
new_data = new_data/np.max(np.abs(new_data))
new_data = np.int16(new_data*32767)
wavfile.write(r"your-file_output.wav", fs, new_data.astype(np.int16))

#%%

plt.close('all')
    
fs, data = wavfile.read(r"your-file_trimmed.wav")
dt = 1/fs
t_old = np.arange(len(data))*dt

df = 0.1 # Hz
T = 1/(2*df) # s
N = int(T/dt)

f, t_filter, Zxx = signal.stft(data, fs = fs, nperseg = N, noverlap = 3*N/4)#, return_onesided = False)
Z_new = np.empty_like(Zxx)
window_length = 3

win_wd = int((window_length - 1)/2)
freq = 1/(2*window_length*np.diff(t_filter)[0])

mult_arr = signal.windows.hamming(window_length)
# mult_arr = np.ones(window_length)
ii_stop = 0

for ii in range(t_filter.shape[0]):
    
    if ii - win_wd < 0:
        mult_start = win_wd - ii
        mult_end = window_length + 1
        start = 0
        end = ii + win_wd + 1

    elif ii + win_wd + 1 > t_filter.shape[0]:
        mult_start = 0 
        mult_end = win_wd + (t_filter.shape[0] - ii)
        start = ii - win_wd
        end = t_filter.shape[0]
        
    else:
        mult_start = 0
        mult_end = window_length + 1
        start = ii - win_wd
        end = ii + win_wd + 1

    Zmultiplier = Zxx[:,start:end]*np.exp(1j*2*np.pi*np.matmul(f[:,np.newaxis], np.transpose(t_filter[ii] - t_filter[start:end,np.newaxis])))
    mult = Zmultiplier * mult_arr[mult_start:mult_end]
    div = np.sum(np.abs(mult), axis = 0)
    mult = mult / div
    
    Z_new[:,ii] = np.prod(1 / mult[:,:win_wd ], axis = 1)*np.prod(1 / mult[:,win_wd:], axis = 1)
    Z_new[:,ii]  = Z_new[:,ii] / np.sum(np.abs(Z_new[:,ii]))
    
t_Z, new_data = np.real(signal.istft(Z_new, fs = fs, nperseg = N, noverlap = 3*N/4))#, input_onesided = False))

    
new_data = new_data/np.max(np.abs(new_data))
new_data = np.int16(new_data*32767)
wavfile.write(r"your-file_output.wav", fs, new_data.astype(np.int16))

