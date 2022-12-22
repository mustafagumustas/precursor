import matplotlib.gridspec as gridspec
from astropy.io import fits
import astropy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from waveFunctions import wavelet, wave_signif

import matplotlib
from astrosc import fits2df, lc_plotter, waveanalysis, wave_contour_plot


def lc_reader(path, plot=True, returns=False):
    df = fits2df(path, real_value=True, time_index=True, TRIGGERTIME=True)
    # df = df.loc[-240:0]

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
                    if 4 > esik_sayaci > 1:

                        count_list = [
                            i + 0.128 if i == max(count_list) else i for i in count_list
                        ]
                        count_list = [
                            i - 0.128 if i == min(count_list) else i for i in count_list
                        ]
                        for i in count_list:
                            if i < 0:
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


df = lc_reader("GRB090510/msec128.lc", plot=True, returns=True)


def wave_analysis(GRB):
    df = fits2df(GRB, real_value=True, time_index=True, TRIGGERTIME=True)
    df = df.loc[-200:-0.128]
    points = df["RATE"]
    # points = points[points.values > 0]
    sst = [float(i) for i in points]
    sst = sst - np.mean(sst)

    # ----------C-O-M-P-U-T-A-T-I-O-N------S-T-A-R-T-S------H-E-R-E---------------

    # variance = np.std(sst, ddof=1) ** 2
    # print("variance = ", variance)
    variance = 1.0
    sst = sst / np.std(sst, ddof=1)

    n = len(sst)
    dt = 0.128
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

    scale = scale[:13]
    signif = signif[:13]
    power = power[:13, :]
    period = period[:13]

    sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
    sig95 = power / sig95  # where ratio > 1, power is significant

    # ------------------------------------------------------ Plotting
    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(8, 16)
    ax_1 = fig.add_subplot(gs[0:8, :12])
    # --- Contour plot wavelet power spectrum
    ax_1.contourf(time, period, power)
    ax_1.contour(time, period, sig95, [-99, 1], colors="k")
    # cone-of-influence, anything "below" is dubious
    ax_1.fill_between(
        time,
        coi * 0 + period[-1],
        coi,
        facecolor="none",
        edgecolor="#00000040",
        hatch="x",
    )
    ax_1.plot(time, coi, "k")
    # format y-scale
    ax_1.set_yscale("log", base=2, subs=None)
    ax_1.set_ylim([np.min(period), np.max(period)])
    ax_1.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax_1.ticklabel_format(axis="y", style="plain")
    ax_1.invert_yaxis()

    ax_2 = fig.add_subplot(gs[0:8, 12:16])
    global_ws = np.sum(power, axis=1) / n  # time-average over all times
    global_ws = global_ws[:13]
    dof = n - scale  # the -scale corrects for padding at edges
    global_signif = wave_signif(
        variance, dt=dt, scale=scale, sigtest=1, lag1=lag1, dof=dof, mother=mother
    )
    global_signif = global_signif[:13]
    ax_2.plot(global_ws, period, "b", label="Wavelet spectrum")
    ax_2.plot(
        global_signif,
        period,
        "r--",
        label=str(0.997 * 100) + "% confidence spectrum",
    )
    ax_2.legend(loc=0)
    ax_2.set_yscale("log")
    ax_2.set_xlabel("Power", fontsize=18)
    ax_2.set_ylim(max(2 * dt, min(period)), max(period))
    ax_2.set_xscale("log")
    plt.show()


# df = fits2df(
#     "GRB081024A/msec128.lc", real_value=True, time_index=True, TRIGGERTIME=True
# )
# df = df.loc[-240:5]
# points = df["VALUES"]
# time, period, power, sig95, coi = waveanalysis(points, points.index)
# # wave_contour_plot(time, period, power, sig95, coi)

# wave_analysis("GRB090510/msec128.lc")
# wave_analysis("GRB081024A/msec128.lc")
# wave_analysis("GRB081024A/msec128_yeni.lc")
plt.show()

# df = fits2df("GRB090510/msec128.lc", real_value=True, time_index=True, TRIGGERTIME=True)
# df = df.loc[-240:0]
# # df["ERROR"] = df["ERROR"] / 2
# df["ERROR"] = np.sqrt(df["ERROR"])
# df["new"] = df["RATE"] - ((df["ERROR"]))

# x = df["new"][df["new"] > 0]
# x = x.dropna()
# # print(x)
# aa = [
#     x.index[i] for i in range(len(x.index) - 1) if (x.index[i + 1] - x.index[i]) < 0.129
# ]
# b = df["VALUES"].loc[aa]
# # print(b)
# c = b.index.to_series().diff()
# c = pd.DataFrame(c[c < 0.129])
# d = c.index.to_series().diff()
# d = d[d<0.129]
# print(c[-20:])
# print(x.loc[i])
# print(x.iloc[i + 1] - x.iloc[i])

# plt.scatter(df.index, df.new.values)

# lc_plotter(df, df["RATE"], error=True)
# plt.show()

# lc_reader("GRB090510/msec128.lc", plot=True)
