from astropy.io import fits
import astropy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from waveFunctions import wavelet, wave_signif

import matplotlib
from astrosc import fits2df, lc_plotter


def lc_reader(path, plot=True, returns=False):
    df = fits2df(path, real_value=True, time_index=True, TRIGGERTIME=True)
    df = df.loc[-20:20]
    # creating new df, thats for picking the below point of each data
    values = pd.DataFrame(df["RATE"] - df["ERROR"], columns=["Values"])
    values = values[values["Values"] > 0]

    esik_sayaci = 0
    count_list = []

    if plot == True:
        for i in range(len(values.index)):
            try:
                fark = values.index.values[i + 1] - values.index.values[i]
                if fark < 0.130:
                    esik_sayaci += 1
                    if values.index.values[i] not in count_list:

                        count_list.append(values.index.values[i])
                    # print(values.index.values[i + 1], values.index.values[i])
                    if esik_sayaci > 1:

                        count_list = [
                            i + 0.128 if i == max(count_list) else i for i in count_list
                        ]
                        count_list = [
                            i - 0.128 if i == min(count_list) else i for i in count_list
                        ]
                        for i in count_list:
                            if i < 0:
                                # print(i)
                                pass
                            else:
                                count_list = []
                        x = [plt.axvline(x=i, color="r") for i in count_list]
                        plt.axvspan(
                            min(count_list),
                            max(count_list),
                            color="r",
                            alpha=0.2,
                        )
                else:
                    esik_sayaci = 0
                    count_list = []
            except:
                pass

    lc_plotter(
        df,
        df["RATE"],
        error=True,
    )


# df = lc_reader("GRB090510/msec128.lc", plot=True, returns=True)


def wave_analysis(GRB):
    df = fits2df(GRB, real_value=True, time_index=True, TRIGGERTIME=True)
    df = df.loc[-20:20]
    points = df["VALUES"]
    points = points[points.values > 0]
    xx = [float(i) for i in points]
    sst = xx
    sst = sst - np.mean(sst)

    # ----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E---------------

    # variance = np.std(sst, ddof=1) ** 2
    # print("variance = ", variance)
    variance = 1.0
    sst = sst / np.std(sst, ddof=1)

    n = len(sst)
    dt = 0.25
    # time = np.arange(len(sst)) * dt  # construct time array
    time = points.index
    pad = 1  # 0 is default
    dj = 0.25
    s0 = 2 * dt
    j1 = np.fix((np.log(n * dt / s0) / np.log(2)) / dj)
    # j1 = 7 / dj
    lag1 = 0.72  # lag-1 autocorrelation for red noise background
    mother = "MORLET"

    # Wavelet transform:
    # wavelet(Y, 0.128, pad=0, dj=-1, s0=-1, J1=-1, mother="MORLET", param=6, freq=None):
    wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, mother)

    power = (np.abs(wave)) ** 2  # compute wavelet power spectrum
    # scale = np.linspace(0.256, 4, 13)

    # Significance levels:
    signif = wave_signif(
        ([variance]), dt=dt, sigtest=0, scale=scale, lag1=lag1, mother=mother
    )  # try with and without lag1 value

    # signif = signif[:13]
    # power = power[:13, :]
    # period = period[:13]

    sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
    sig95 = power / sig95  # where ratio > 1, power is significant

    # ------------------------------------------------------ Plotting

    # --- Contour plot wavelet power spectrum
    plt.contourf(time, period, power)

    plt.contour(time, period, sig95, [-99, 1], colors="k")
    # cone-of-influence, anything "below" is dubious
    plt.fill_between(
        time,
        coi * 0 + period[-1],
        coi,
        facecolor="none",
        edgecolor="#00000040",
        hatch="x",
    )
    plt.plot(time, coi, "k")
    # format y-scale
    plt.yscale("log", base=2, subs=None)
    plt.ylim([np.min(period), np.max(period)])
    ax = plt.gca()
    ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.ticklabel_format(axis="y", style="plain")
    ax.invert_yaxis()

    plt.show()


wave_analysis("GRB090510/msec128.lc")
