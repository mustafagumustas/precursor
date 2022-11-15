from astropy.io import fits
import astropy
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class lc_reader:
    def __init__(self, path):
        # hdul = fits.open("GRB090510/msec128.lc")
        self.hdul = fits.open(path)

        # getting trigger time from header
        self.trigger_t = self.hdul[1].header["TRIGTIME"]

        # we need to multuply the count rates with NGOODPIX value
        # its a must EXPLAIN IT LATER

        self.NGOODPIX = self.hdul[1].header["NGOODPIX"]

        # creating dataframe to work on the data better
        self.df = pd.DataFrame(self.hdul[1].data)

        # time values
        self.df["TIME"] = self.df["TIME"] - self.trigger_t

        # count values
        self.count = self.df["RATE"]
        self.count.index = self.df["TIME"]

        # error values
        self.error = self.df["ERROR"]
        self.error.index = self.df["TIME"]

        self.count = self.count * self.NGOODPIX
        self.error = self.error * self.NGOODPIX

        # filtering to see the precursors for GRB090510
        self.count = self.count[(-20 < self.count.index) & (self.count.index < 20)]
        self.error = self.error[(-20 < self.error.index) & (self.error.index < 20)]

        self.points = pd.DataFrame(self.count - self.error, columns=["foo"])
        self.points = self.points[self.points["foo"] > 0]
        # print(self.points)
        esik_sayaci = 0
        count_list = []
        for i in range(len(self.points.index)):
            try:
                fark = self.points.index.values[i + 1] - self.points.index.values[i]
                if fark < 0.130:
                    esik_sayaci += 1
                    if self.points.index.values[i] not in count_list:

                        count_list.append(self.points.index.values[i])
                    # print(self.points.index.values[i + 1], self.points.index.values[i])
                    if 4 > esik_sayaci > 1:
                        print(count_list, end="\n\n")
                else:
                    esik_sayaci = 0
                    count_list = []
            except:
                pass

        # plt.scatter((self.higher_p.index - 0.065), self.higher_p.values, color="r")
        # plt.step(self.count.index, self.count.values)
        # plt.axhline(y=0, color="r")
        # plt.errorbar(
        #     (self.count.index - 0.065),
        #     self.count.values,
        #     yerr=self.error.values,
        #     ls="none",
        #     ecolor="black",
        #     elinewidth=0.5,
        # )
        # plt.show()


# lc_reader("GRB090510/msec128.lc")

from waveFunctions import wavelet


hdul = fits.open("GRB090510/msec128.lc")

# getting trigger time from header
trigger_t = hdul[1].header["TRIGTIME"]

# we need to multuply the count rates with NGOODPIX value
# its a must EXPLAIN IT LATER

NGOODPIX = hdul[1].header["NGOODPIX"]

# creating dataframe to work on the data better
df = pd.DataFrame(hdul[1].data)

# time values
df["TIME"] = df["TIME"] - trigger_t

# count values
count = df["RATE"]
count.index = df["TIME"]

# error values
error = df["ERROR"]
error.index = df["TIME"]

count = count * NGOODPIX
error = error * NGOODPIX

# filtering to see the precursors for GRB090510
count = count[(-20 < count.index) & (count.index < 20)]
error = error[(-20 < error.index) & (error.index < 20)]

points = pd.DataFrame(count - error, columns=["foo"])
points = points[points["foo"] > 0]
xx = [float(i) for i in points.values]
# wavelet(Y, 0.128, pad=0, dj=-1, s0=-1, J1=-1, mother="MORLET", param=6, freq=None):
yy = wavelet(xx, 0.128)

# print(yy[3])
# for i in yy:
#     print(i.shape, end="\n\n\n")
plt.plot(yy[0])
plt.show()
