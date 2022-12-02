from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd

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

    NGOODPIX = hdul[1].header["NGOODPIX"]  # NGOODPIX value
    # we need to multuply the count rates with NGOODPIX value
    # in order to get the real values, its a must
    df["RATE"] = df["RATE"] * NGOODPIX

    if real_value:
        # subtract error from count values and add new col to df as values
        df["ERROR"] = df["ERROR"] * NGOODPIX
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


# df = fits2df("GRB090510/msec128.lc", real_value=True, time_index=True, TRIGGERTIME=True)


# lc_plotter(
#     df,
#     df["RATE"],
#     error=True,
# )
