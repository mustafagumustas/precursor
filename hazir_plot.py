from numpy import *
from scipy import *
import pandas as pd
from pylab import *
import pywt
import scaleogram as scg
from scipy.stats import chi2
import matplotlib.gridspec as gridspec
from matplotlib import ticker, cm
from matplotlib.colors import LogNorm
import time
from tqdm import tqdm
import urllib


def plot_wavelet(
    times,
    signal,
    scales,
    variance,
    alpha,
    waveletname,
    p_range=[0, 10],
    ylabel="Period (s)",
    xlabel="Time (s)",
    power_scale="log",
):

    # dt = times[1] - times[0]
    [coefficients, frequencies] = pywt.cwt(signal, scales, waveletname, dt)
    power = (abs(coefficients)) ** 2

    time1 = []
    signal1 = []
    power1 = []
    period = 1.0 / frequencies
    for i in range(0, len(times)):
        if (times[i] >= p_range[0]) and (times[i] <= p_range[1]):
            time1.append(times[i])
            signal1.append(signal[i])
            power1.append(power[:, i])
    power = np.array(power1)
    signal = np.array(signal1)
    times = np.array(time1)
    power = power.transpose()

    significance_level = 0.997
    n0 = len(signal)
    freq = frequencies
    dofmin = 2
    pk = lambda k, a, N: (1 - a**2) / (1 + a**2 - 2 * a * cos(2 * pi * k / N))
    fft_theor = pk(freq, alpha, n0)
    fft_theor = variance * fft_theor  # Including time-series variance
    signif = fft_theor
    dof = dofmin
    chisquare = chi2.ppf(significance_level, dof) / dof
    signif = fft_theor * chisquare

    sig95 = (signif * np.ones((len(signal), 1))).transpose()
    sig95 = power / sig95

    global_ws = np.sum(power.conj().transpose(), axis=0) / n0

    global_signif = fft_theor * chisquare

    for i in range(0, len(power)):
        power[i, :] = power[i, :] / (np.mean(power[i, :]))

    fig = plt.figure(figsize=(10, 8))

    gs = gridspec.GridSpec(12, 12)
    ax = fig.add_subplot(gs[4:12, :8])

    ax_1 = fig.add_subplot(gs[0:4, :8])
    plt.title(title)
    ax_1.xaxis.set_visible(False)
    #     levels = [1, 1e1, 1e2, 1e3, 1e4, 1e5]
    #     im = ax.contourf(times, frequencies, log10(power), log10(levels),locator=ticker.LogLocator(), extend='both',cmap=cmap)
    if power_scale == "log":
        im = ax.pcolor(
            times,
            period,
            power,
            norm=LogNorm(vmin=power.min(), vmax=power.max()),
        )
    else:
        im = ax.pcolor(times, period, power)
    #     ax.contour(times,period, sig95, [-99, 1], colors='black',
    #            linewidths=2.)
    ax.set_ylabel(ylabel, fontsize=18)
    ax_1.set_ylabel("Counts", fontsize=18)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_yscale("log")
    ax.set_ylim(max(2 * dt, min(period)), max(period))
    #     ax.invert_yaxis()

    cbar_ax = fig.add_axes([0, 0.15, 0.01, 0.4])
    fig.colorbar(im, cax=cbar_ax, orientation="vertical")

    ax_1.set_xlim(times[0], times[-1])
    pic = ax_1.plot(times, signal)

    ax_2 = fig.add_subplot(gs[0:3, 9:12])
    wav = pywt.ContinuousWavelet(waveletname)
    fun_wav, times = wav.wavefun()
    ax_2.plot(times, fun_wav.real)
    ax_2.set_xlabel("Time (s)")
    ax_2.set_ylabel("Amplitude")
    ax_2.set_ylim(-1, 1)
    ax_2.set_title("Wavelet " + waveletname)

    xf = fftfreq(n0, dt)
    fft_power = np.abs(np.fft.fft(signal)) ** 2

    ax_3 = fig.add_subplot(gs[4:12, 8:12])
    plt.setp(ax_3.get_yticklabels(), visible=False)
    plt.setp(ax_1.get_xticklabels(), visible=False)
    #     ax_3.plot(fft_power, xf, 'gray', label='Fourier spectrum')
    ax_3.plot(global_ws, period, "b", label="Wavelet spectrum")
    ax_3.plot(
        global_signif,
        period,
        "r--",
        label=str(significance_level * 100) + "% confidence spectrum",
    )
    ax_3.legend(loc=0)
    ax_3.set_yscale("log")
    ax_3.set_xlabel("Power", fontsize=18)
    ax_3.set_ylim(max(2 * dt, min(period)), max(period))
    ax_3.set_xscale("log")


#     ax_3.invert_yaxis()
#     ax_3.set_xlim(0,max(max(global_ws),max(fft_power)))


#  with open("BBH_events.txt") as mytxt:
# for line in mytxt:

from astrosc import fits2df
from waveFunctions import wave_signif, wavelet

# df = fits2df("GRB090510/msec128.lc", real_value=True, time_index=True, TRIGGERTIME=True)/
# df = fits2df("GRB081024A/msec128.lc", real_value=True, time_index=True, TRIGGERTIME=True)
df = fits2df(
    "GRB081024A/msec128_yeni.lc", real_value=True, time_index=True, TRIGGERTIME=True
)

df = df.loc[-240:-0.128]
sst = [float(i) for i in df["VALUES"]]
sst = sst - np.mean(sst)
final = sst
# variance = 1.0
sst = sst / np.std(sst, ddof=1)
n = len(sst)
dt = 0.25


times = df.index
signal = sst
binning = 1
binned_time = []
binned_signal = []
for i in range(0, len(signal) // (binning * 20) - 1):
    binned_time.append(
        np.mean(times[(i * binning * 20) : ((i + 1) * binning * 20) + 1])
    )
    binned_signal.append(
        np.sum(signal[(i * binning * 20) : ((i + 1) * binning * 20) + 1])
    )

variance = np.var(binned_signal)
alpha = np.corrcoef(binned_signal[0:-1], binned_signal[1:])[0, 1]
# alpha = 0.72
# wavelet_name = "cmor1-1"
# scales = scg.periods2scales(logspace(np.log10(1), log10(100), 200))
dj = 0.25
wavelet_name = "morl"

dt = 0.25
s0 = 2 * dt  # this says start at a scale of 6 months
j1 = np.fix((np.log(n * dt / s0) / np.log(2)) / dj)
wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, "MORLET")
plot_wavelet(
    binned_time,
    binned_signal,
    scale,
    variance,
    alpha,
    wavelet_name,
    p_range=[-500, 500],
    power_scale="linear",
)
# plt.savefig('wavelet_BBH_'+trigger_time[0:11]+'_'+wavelet_name)
plt.show()
# time.sleep(60)
