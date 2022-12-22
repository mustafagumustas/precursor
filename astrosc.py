from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd
from waveFunctions import wavelet, wave_signif
import numpy as np
import matplotlib

__author__ = "Mustafa Gumustas"


def fits2df(
    fits_path,
    table_n=1,
    real_value=False,
    time_index=False,
    TRIGGERTIME=False,
):
    # Reads the fits file and creates DataFrame (df) out of it
    # It could read any fits file and return as df by default.

    # table_n is the index of value table inside fits file, usually its 1 in most of the light curve files, but user can change that.

    # If real_value is True, there will be "VALUES" column added to df that is created by value - error values.

    # If time_index True, the result df has the time from fits file as index to result df

    # NGOODPIX and TRIGGERTIME values taken from header if given as true parameter.

    # read fits file with given path
    hdul = fits.open(fits_path)
    # creating dataframe to work on the data better
    df = pd.DataFrame(hdul[table_n].data)

    if TRIGGERTIME:
        trigger_t = hdul[1].header["TRIGTIME"]  # trigger time
        # fixing time values, making trigger at 0
        df["TIME"] = df["TIME"] - trigger_t  # time values

    # NGOODPIX = hdul[1].header["NGOODPIX"]  # NGOODPIX value
    # we need to multuply the count rates with NGOODPIX value
    # in order to get the real values, its a must
    # df["RATE"] = df["RATE"] * NGOODPIX

    if real_value:
        # subtract error from count values and add new col to df as values
        # df["ERROR"] = df["ERROR"] * NGOODPIX
        df["VALUES"] = df["RATE"] - df["ERROR"]

    if time_index:
        # set time as index of DF
        df = df.set_index("TIME")

    return df


def lc_plotter(df, x, error=False):
    plt.step(df.index, x.values)
    plt.axhline(y=0, color="r")
    if error == True:
        col_names = df.columns
        error_col = input(
            f"Please select the error column name from your df \n{[i for i in col_names]}:"
        )
        time_span = float(input("Plese give the time span: "))
        print(df[error_col])
        plt.errorbar(
            (df.index.values - (time_span / 2)),
            x,
            yerr=df[error_col].values,
            ls="none",
            ecolor="black",
            elinewidth=0.5,
        )
    plt.show()


def waveanalysis(df, time, mother="MORLET"):
    # gets values data from user and makes wave analysis using waveFunctions library

    # values_col = input("Please enter value column name: ")
    # points = df[values_col]
    points = df
    # points = points[points.values > 0]
    sst = [float(i) for i in points]
    sst = sst - np.mean(sst)

    # variance = np.std(sst, ddof=1) ** 2
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
    lag1 = 0.72  # lag-1 autocorrelation for red noise background
    mother = "MORLET"

    # Wavelet transform:
    wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, mother)

    power = (np.abs(wave)) ** 2  # compute wavelet power spectrum

    # on article they use from 2 dt to 4 with 13 time scales
    # scale = np.linspace(0.256, 4, 13)

    # Significance levels:
    signif = wave_signif(
        ([variance]), dt=dt, sigtest=0, scale=scale, lag1=lag1, mother=mother
    )  # try with and without lag1 value

    sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
    sig95 = power / sig95  # where ratio > 1, power is significant

    return time, period, power, sig95, coi


def wave_contour_plot(time, period, power, sig95, coi):
    # user must plot after function
    levels = [0, 0.5, 1, 2, 4, 999]
    CS = plt.contourf(time, period, power, len(levels))
    im = plt.contourf(
        CS, levels=levels, colors=["white", "bisque", "orange", "orangered", "darkred"]
    )
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


def new_wave(df, signal, scales, waveletname, dt):
    time = df.index


# hdul = fits.open("GRB090510/msec128_ucuncu.lc")
# count = pd.DataFrame([i[1] for i in hdul[1].data])
# time = pd.DataFrame([i[0] for i in hdul[1].data])
# TRIGTIME = hdul[1].header["TRIGTIME"]
# time = time - TRIGTIME

# # creating count and error values in DF
# count = pd.DataFrame([i[1] for i in hdul[1].data], index=time[0])

# # plt.step(count.index, count[0])
# total = count[0]+count[1]+count[2]+count[3]
# # print(total)
# total = total.loc[-20:0]
# plt.step(total.index, total)
# plt.show()
