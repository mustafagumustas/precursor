from astropy.io import fits
import astropy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from waveFunctions import wavelet, wave_signif


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
    count = count[(-20 < count.index) & (count.index < 20)]
    error = error[(-20 < error.index) & (error.index < 20)]

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
    wave, period, scale, coi = wavelet(xx, 0.128)
    wave = (np.absolute(wave)) ** 2
    n = len(sst)
    # wave_signif(Y,dt,scale,sigtest=0,lag1=0.0,siglvl=0.95,dof=None,mother="MORLET",param=None,gws=None)
    signif = wave_signif(wave, 0.128, scale=scale)
    print(signif)
    # sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])
    # sig95 = wave / sig95

    # limitations
    # sig95 = sig95[:13]
    # wave = wave[:13, :]
    # period = period[:13]
    # coi = coi[:13]
    # scale = np.linspace(0.256, 4, 13)
    # time = points.index

    contour = plt.contourf(points.index, scale, wave)
    # plt.fill_between(
    #     time,
    #     coi * 0 + period[-1],
    #     coi,
    #     facecolor="none",
    #     edgecolor="#00000040",
    #     hatch="x",
    # )
    plt.gca().invert_yaxis()
    plt.plot(points.index, coi, "k")
    plt.ylabel("scale")
    plt.xlabel("time")
    plt.title("GRB090510 Wavelet Analysis \n(Morlet Mothre Function)")
    # plt.plot(signif)
    # plt.ylim(top=4)
    # plt.show()


# plot between 2*0.128 and 4

# sig %95
wave_analysis("GRB090510/msec128.lc")
