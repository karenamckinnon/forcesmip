# script written for ~/Save/ForceSMIP/evaluation.trends.1mem.plot.ncl
import os
import xarray as xr 
from scipy import signal

def lpf_modata(dat, period=10):
    fs=1/(30.41667*24*3600)        # 1 month in Hz (sampling frequency)
    nyquist = fs / 2          # 0.5 times the sampling frequency
    cutoff = fs/period            # cutoff frequency
    cutoff = cutoff/nyquist   # as fraction of nyquist  
    print('cutoff= ',(1/(cutoff*nyquist))/(30.41667*24*3600),' months') 
    filtsos = signal.butter(4, cutoff, 'lowpass', output='sos') #low pass filter
    filtb, filta = signal.butter(4, cutoff, 'lowpass')
    dat_out = xr.apply_ufunc(signal.sosfiltfilt, filtsos,dat.fillna(0),kwargs={'padtype':'even','axis':0}).where(dat.notnull())
    #dat_out = xr.apply_ufunc(signal.sosfiltfilt, filtsos,dat.fillna(0),kwargs={'padtype':'odd','axis':1}).where(dat.notnull())
    print('Function call completed')
    #return np.float32(dat_out)
    return dat_out


ds = xr.open_dataset('tmp.nc',decode_times=False)
arr = ds['aa']

arr_lpf = lpf_modata(arr, period=120)

arr_lpf.to_netcdf('tmp2.nc')
