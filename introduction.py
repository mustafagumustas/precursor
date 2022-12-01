from astropy.io import fits
import astropy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from waveFunctions import wavelet, wave_signif

import matplotlib


def fits2df(fits_path, error=True, time_index=True, NGOODPIX=True, TRIGGERTIME=True):
    # this function created to read fits functions and make DataFrame out of it
    # it has optimized for Astronomical satellite datas, it looks for Trigger
    # time and NGOODPIX from header then subtracks the related data, user can block that functions
    # function also subtracks the error out of the count by default

    # read fits file with given path
    hdul = fits.open(fits_path)
    # creating dataframe to work on the data better
    df = pd.DataFrame(hdul[1].data)

    if TRIGGERTIME:
        trigger_t = hdul[1].header["TRIGTIME"]  # trigger time
        # fixing time values, making trigger at 0
        df["TIME"] = df["TIME"] - trigger_t  # time values

    if NGOODPIX:
        NGOODPIX = hdul[1].header["NGOODPIX"]  # NGOODPIX value
        # we need to multuply the count rates with NGOODPIX value
        # in order to get the real values, its a must
        df["RATE"] = df["RATE"] * NGOODPIX

    if error:
        # subtract error from count values and add new col to df as values
        df["ERROR"] = df["ERROR"] * NGOODPIX
        df["VALUES"] = df["RATE"] - df["ERROR"]

    if time_index:
        # set time as index of DF
        df = df.set_index("TIME")

    return df


def lc_reader(path, plot=True, returns=False):
    df = fits2df(path)

    count = df["RATE"]
    error = df["ERROR"]
    count = count[(-20 < count.index) & (count.index < 0)]
    error = error[(-20 < error.index) & (error.index < 0)]

    # creating new df, thats for picking the below point of each data
    values = pd.DataFrame(count - error, columns=["Values"])
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
        plt.step(count.index, count.values)
        # print(count_list)

        # print(min(i))
        plt.axhline(y=0, color="r")
        plt.errorbar(
            (count.index - 0.065),
            count.values,
            yerr=error.values,
            ls="none",
            ecolor="black",
            elinewidth=0.5,
        )
        plt.show()
    if returns == True:
        return (
            pd.DataFrame()
            .assign(
                RATE=count.values,
                ERROR=error.values,
                VALUES=(count.values - error.values),
            )
            .set_index(count.index)
        )


# df = lc_reader("GRB090510/msec128.lc", plot=False, returns=True)


def wave_analysis(GRB):
    df = lc_reader(GRB, plot=False, returns=True)
    points = df["VALUES"]
    points = points[points.values > 0]
    xx = [float(i) for i in points]
    sst = xx
    count = df["RATE"].values
    error = df["ERROR"].values
    points = df["VALUES"]
    points = points[points.values > 0]
    xx = [float(i) for i in points]

    # wavelet(Y, 0.128, pad=0, dj=-1, s0=-1, J1=-1, mother="MORLET", param=6, freq=None):
    sst = sst - np.mean(sst)
    variance = np.std(sst, ddof=1) ** 2
    print("variance = ", variance)

    # ----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E---------------

    variance = 1.0
    sst = sst / np.std(sst, ddof=1)
    n = len(sst)
    dt = 0.25
    # time = np.arange(len(sst)) * dt  # construct time array
    time = points.index
    # xlim = [1870, 2000]  # plotting range
    pad = 1  # 0 is default
    dj = 0.25
    s0 = 2 * dt
    j1 = np.fix((np.log(n * dt / s0) / np.log(2)) / dj)
    # j1 = 7 / dj
    lag1 = 0.72  # lag-1 autocorrelation for red noise background
    mother = "MORLET"

    # Wavelet transform:
    wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, mother)
    power = (np.abs(wave)) ** 2  # compute wavelet power spectrum

    # Significance levels:
    signif = wave_signif(
        ([variance]), dt=dt, sigtest=0, scale=scale, lag1=lag1, mother=mother
    )  # try with and without lag1 value
    # scale = np.linspace(0.256, 4, 13)
    # signif = signif[:13]
    # power = power[:13, :]
    # period = period[:13]
    sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
    sig95 = power / sig95  # where ratio > 1, power is significant

    # ------------------------------------------------------ Plotting

    # --- Contour plot wavelet power spectrum
    CS = plt.contourf(time, period, power)
    # im = plt.contourf(CS, colors=["white", "bisque", "orange", "orangered", "darkred"])
    # 95# significance contour, levels at -99 (fake) and 1 (95# signif)
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
