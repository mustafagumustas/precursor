import matplotlib.pylab as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec

import numpy as np

# from mpl_toolkits.axes_grid1 import make_axes_locatable

from waveFunctions import wave_signif, wavelet


# WAVETEST Example Python script for WAVELET, using NINO3 SST dataset

from astropy.io import fits
import pandas as pd
from introduction import lc_reader

df = lc_reader("GRB090510/msec128.lc", plot=False, returns=True)

count = df["RATE"].values
error = df["ERROR"].values
points = df["VALUES"]
points = points[points.values > 0]
xx = [float(i) for i in points]

# wavelet(Y, 0.128, pad=0, dj=-1, s0=-1, J1=-1, mother="MORLET", param=6, freq=None):
# wave, period, scale, coi = wavelet(xx, 0.128)
# sst = (np.absolute(wave)) ** 2

sst = xx
sst = sst - np.mean(sst)
variance = np.std(sst, ddof=1) ** 2
print("variance = ", variance)

# ----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E---------------

# normalize by standard deviation (not necessary, but makes it easier
# to compare with plot on Interactive Wavelet page, at
# "http://paos.colorado.edu/research/wavelets/plot/"
if 0:
    variance = 1.0
    sst = sst / np.std(sst, ddof=1)
n = len(sst)
dt = 0.25
# time = np.arange(len(sst)) * dt  # construct time array
time = points.index
# xlim = [1870, 2000]  # plotting range
pad = 1  # pad the time series with zeroes (recommended)
dj = 0.25  # this will do 4 sub-octaves per octave
s0 = 2 * dt  # this says start at a scale of 6 months
j1 = 7 / dj  # this says do 7 powers-of-two with dj sub-octaves each
lag1 = 0.72  # lag-1 autocorrelation for red noise background
print("lag1 = ", lag1)
mother = "MORLET"

# Wavelet transform:
wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, mother)
power = (np.abs(wave)) ** 2  # compute wavelet power spectrum

# Significance levels:
signif = wave_signif(
    ([variance]), dt=dt, sigtest=0, scale=scale, lag1=lag1, mother=mother
)
# scale = np.linspace(0.256, 4, 13)
# signif = signif[:13]
# power = power[:13, :]
# period = period[:13]
# expand signif --> (J+1)x(N) array
sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
sig95 = power / sig95  # where ratio > 1, power is significant


# ------------------------------------------------------ Plotting

# --- Plot time series
fig = plt.figure(figsize=(9, 10))
gs = GridSpec(1, 1, hspace=0.4, wspace=0.75)
plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.95, wspace=0, hspace=0)
plt.subplot(gs[0, 0:3])
# plt.plot(time, sst, "k")

# --- Contour plot wavelet power spectrum
plt3 = plt.subplot(gs[0, 0])
levels = [0, 0.5, 1, 2, 4, 999]
# *** or use 'contour'
CS = plt.contourf(time, period, power, len(levels))
im = plt.contourf(
    CS, levels=levels, colors=["white", "bisque", "orange", "orangered", "darkred"]
)
# 95# significance contour, levels at -99 (fake) and 1 (95# signif)
print(sig95.shape)
plt.contour(time, period, sig95, [-99, 1], colors="k")
# cone-of-influence, anything "below" is dubious
# plt.fill_between(
#     time, coi * 0 + period[-1], coi, facecolor="none", edgecolor="#00000040", hatch="x"
# )
# plt.plot(time, coi, "k")
# format y-scale
plt3.set_yscale("log", base=2, subs=None)
plt.ylim([np.min(period), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(ticker.ScalarFormatter())
plt3.ticklabel_format(axis="y", style="plain")
plt3.invert_yaxis()

plt.show()
