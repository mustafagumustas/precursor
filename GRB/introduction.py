from astropy.io import fits
import astropy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from waveFunctions import wavelet, wave_signif


def lc_reader(path, plot=True, returns=False):
    # read fits file with given path
    hdul = fits.open(path)

    # Getting info from header
    trigger_t = hdul[1].header["TRIGTIME"]  # trigger time
    NGOODPIX = hdul[1].header["NGOODPIX"]  # NGOODPIX value

    # creating dataframe to work on the data better
    df = pd.DataFrame(hdul[1].data)

    # fixing time values, making trigger at 0
    df["TIME"] = df["TIME"] - trigger_t  # time values

    count = df["RATE"]  # count values
    error = df["ERROR"]  # error values

    # set time as index of DF
    count.index = df["TIME"]
    error.index = df["TIME"]

    # we need to multuply the count rates with NGOODPIX value
    # in order to get the real values, its a must
    count = count * NGOODPIX
    error = error * NGOODPIX

    # filtering between -20 and 20 sec
    count = count[(-20 < count.index) & (count.index < 20)]
    error = error[(-20 < error.index) & (error.index < 20)]

    # creating new df, thats for picking the below point of each data
    points = pd.DataFrame(count - error, columns=["foo"])
    points = points[points["foo"] > 0]

    esik_sayaci = 0
    count_list = []
    for i in range(len(points.index)):
        try:
            fark = points.index.values[i + 1] - points.index.values[i]
            if fark < 0.130:
                esik_sayaci += 1
                if points.index.values[i] not in count_list:

                    count_list.append(points.index.values[i])
                # print(points.index.values[i + 1], points.index.values[i])
                if 4 > esik_sayaci > 1:

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

    # plt.scatter((higher_p.index - 0.065), higher_p.values, color="r")
    if plot == True:
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
        return pd.DataFrame(count - error)


# df = lc_reader("GRB090510/msec128.lc", plot=False, returns=True)


def wave_analysis(GRB):
    hdul = fits.open(GRB)
    trigger_t = hdul[1].header["TRIGTIME"]
    NGOODPIX = hdul[1].header["NGOODPIX"]

    df = pd.DataFrame(hdul[1].data)

    df["TIME"] = df["TIME"] - trigger_t  # time values

    count = df["RATE"]  # count values
    error = df["ERROR"]  # error values

    count.index = df["TIME"]
    error.index = df["TIME"]

    count = count * NGOODPIX
    error = error * NGOODPIX

    count = count[(-20 < count.index) & (count.index < 20)]
    error = error[(-20 < error.index) & (error.index < 20)]

    points = pd.DataFrame(count - error, columns=["foo"])
    points = points[points["foo"] > 0]
    xx = [float(i) for i in points.values]

    # wavelet(Y, 0.128, pad=0, dj=-1, s0=-1, J1=-1, mother="MORLET", param=6, freq=None):
    wave, period, scale, coi = wavelet(xx, 0.128)

    # wave_signif(Y,dt,scale,sigtest=0,lag1=0.0,siglvl=0.95,dof=None,mother="MORLET",param=None,gws=None)
    signif = wave_signif(wave, 0.128, scale)

    # print(signif.shape)
    # scale = np.linspace(0.256, 4, 13)

    wave = (np.absolute(wave)) ** 2
    signif = wave_signif(points.index, 0.128, scale)
    # print(signif)
    # print(wave)

    plt.contour(points.index, scale, wave)
    # plt.plot(signif)
    # plt.ylim(top=4)
    plt.show()


# plot between 2*0.128 and 4

# sig %95
wave_analysis("GRB090510/msec128.lc")
