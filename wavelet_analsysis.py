import matplotlib.pylab as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec

import numpy as np

# from mpl_toolkits.axes_grid1 import make_axes_locatable

from waveFunctions import wave_signif, wavelet

__author__ = "Evgeniya Predybaylo"

from astrosc import fits2df

# WAVETEST Example Python script for WAVELET
# READ THE DATA

df = fits2df("GRB090510/msec128.lc", real_value=True, time_index=True, TRIGGERTIME=True)
df = df.loc[-20:10]
sst = [float(i) for i in df["VALUES"]]
sst = sst - np.mean(sst)


# ----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E---------------

# normalize by standard deviation (not necessary, but makes it easier
# to compare with plot on Interactive Wavelet page, at
# "http://paos.colorado.edu/research/wavelets/plot/"

variance = 1.0
sst = sst / np.std(sst, ddof=1)
n = len(sst)
dt = 0.25
time = df.index

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
global_ws = np.sum(power, axis=1) / n  # time-average over all times

# Significance levels:
signif = wave_signif(
    ([variance]), dt=dt, sigtest=0, scale=scale, lag1=lag1, mother=mother
)
# expand signif --> (J+1)x(N) array
sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
sig95 = power / sig95  # where ratio > 1, power is significant

# Global wavelet spectrum & significance levels:
dof = n - scale  # the -scale corrects for padding at edges
global_signif = wave_signif(
    variance, dt=dt, scale=scale, sigtest=1, lag1=lag1, dof=dof, mother=mother
)

# Scale-average between El Nino periods of 2--8 years
avg = np.logical_and(scale >= 2, scale < 8)
Cdelta = 0.776  # this is for the MORLET wavelet
# expand scale --> (J+1)x(N) array
scale_avg = scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
scale_avg = power / scale_avg  # [Eqn(24)]
scale_avg = dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
scaleavg_signif = wave_signif(
    variance, dt=dt, scale=scale, sigtest=2, lag1=lag1, dof=([2, 7.9]), mother=mother
)

# ------------------------------------------------------ Plotting

# --- Plot time series
fig = plt.figure(figsize=(10, 8))
gs = GridSpec(3, 4, hspace=0.4, wspace=0.75)
plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.95, wspace=0, hspace=0)
plt.subplot(gs[0, 0:3])
plt.step(time, sst, "k")
print(sst)
plt.axhline(y=0, color="r")
plt.xlabel("Time ")
plt.ylabel("count")
plt.title("GRB090510 msec128 LC")


# --- Contour plot wavelet power spectrum
# plt3 = plt.subplot(3, 1, 2)
plt3 = plt.subplot(gs[1, 0:3])
levels = [0, 0.5, 1, 2, 4, 999]
# *** or use 'contour'
CS = plt.contourf(time, period, power, len(levels))
im = plt.contourf(
    CS, levels=levels, colors=["white", "bisque", "orange", "orangered", "darkred"]
)
plt.xlabel("scale")
plt.ylabel("Period")
plt.title("Wavelet Power Spectrum")
# 95# significance contour, levels at -99 (fake) and 1 (95# signif)
plt.contour(time, period, sig95, [-99, 1], colors="k")
# cone-of-influence, anything "below" is dubious
plt.fill_between(
    time, coi * 0 + period[-1], coi, facecolor="none", edgecolor="#00000040", hatch="x"
)
plt.plot(time, coi, "k")
# format y-scale
plt3.set_yscale("log", base=2, subs=None)
plt.ylim([np.min(period), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(ticker.ScalarFormatter())
plt3.ticklabel_format(axis="y", style="plain")
plt3.invert_yaxis()

# --- Plot global wavelet spectrum
plt4 = plt.subplot(gs[1, -1])
plt.plot(global_ws, period)
plt.plot(global_signif, period, "--")
plt.xlabel("Power")
plt.title("Wavelet Spectrum")
# format y-scale
plt4.set_yscale("log", base=2, subs=None)
plt.ylim([np.min(period), np.max(period)])
ax = plt.gca().yaxis
ax.set_major_formatter(ticker.ScalarFormatter())
plt4.ticklabel_format(axis="y", style="plain")
plt4.invert_yaxis()

# --- Plot 2--8 yr scale-average time series
plt.subplot(gs[2, 0:3])
plt.plot(time, scale_avg, "k")

plt.xlabel("Time")
plt.ylabel("Avg variance")
plt.title("Scale-average Time Series")

plt.show()
