import matplotlib.pyplot as plt
import matplotlib

import numpy as np
import argparse
import warnings

import astropy.units as u
from astropy.modeling import Fittable1DModel, Parameter

# from astropy.modeling.models import PowerLaw1D, Drude1D
from astropy.modeling.fitting import LevMarLSQFitter

from measure_extinction.extdata import ExtData
from dust_extinction.helpers import _get_x_in_wavenumbers, _test_valid_x_range
from dust_extinction.conversions import AxAvToExv

from models_mcmc_extension import EmceeFitter


def drude(x, scale, x_o, gamma_o):
    """
    Drude to play with
    """
    y = (
        scale
        * ((gamma_o / x_o) ** 2)
        / ((x / x_o - x_o / x) ** 2 + (gamma_o / x_o) ** 2)
    )
    return y


def drude_asym(x, scale, x_o, gamma_o, asym):
    """
    Drude to play with
    """
    gamma = 2.0 * gamma_o / (1.0 + np.exp(asym * (x - x_o)))
    y = scale * ((gamma / x_o) ** 2) / ((x / x_o - x_o / x) ** 2 + (gamma / x_o) ** 2)
    return y


class G20(Fittable1DModel):
    """
    Powerlaw plus Drude profiles for the silicate features for the
    1 to 40 micron A(lambda)/A(V) extinction curve.

    Powerlaw portion based on Martin & Whittet (1990).

    Parameters
    ----------
    scale: float
        amplitude of the curve

    alpha : float
        power of powerlaw

    TBD: parameters for silicate features
    """

    # inputs = ("x",)
    # outputs = ("axav",)

    scale = Parameter(description="amplitude", default=0.5)
    alpha = Parameter(description="alpha (power of powerlaw)", default=1.8)
    sil1_amp = Parameter(default=0.05, min=0.01)
    sil1_center = Parameter(default=10.0, bounds=(8.0, 12.0))
    sil1_fwhm = Parameter(default=2.0, bounds=(0.01, 5.0))
    sil2_amp = Parameter(default=0.1, min=0.001)
    sil2_center = Parameter(default=20.0, bounds=(16.0, 24.0))
    sil2_fwhm = Parameter(default=3.0, bounds=(0.01, 10.0))

    x_range = [1.0 / 40.0, 1.0]

    def evaluate(
        self,
        in_x,
        scale,
        alpha,
        sil1_amp,
        sil1_center,
        sil1_fwhm,
        sil2_amp,
        sil2_center,
        sil2_fwhm,
    ):
        """
        G20 function

        Parameters
        ----------
        in_x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in wavenumbers [1/micron]

        Returns
        -------
        axav: np array (float)
            A(x)/A(V) extinction curve [mag]

        Raises
        ------
        ValueError
           Input x values outside of defined range
        """
        x = _get_x_in_wavenumbers(in_x)

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, self.x_range, "G20")

        # powerlaw
        axav = scale * ((1.0 / x) ** (-1.0 * alpha))

        # silicate feature drudes
        wave = 1 / x
        axav += drude(wave, sil1_amp, sil1_center, sil1_fwhm)
        axav += drude(wave, sil2_amp, sil2_center, sil2_fwhm)

        return axav


class G20_drude_asym(Fittable1DModel):
    """
    Powerlaw plus Drude profiles for the silicate features for the
    1 to 40 micron A(lambda)/A(V) extinction curve.

    Powerlaw portion based on Martin & Whittet (1990).

    Parameters
    ----------
    scale: float
        amplitude of the curve

    alpha : float
        power of powerlaw

    TBD: parameters for silicate features
    """

    # inputs = ("x",)
    # outputs = ("axav",)

    scale = Parameter(description="amplitude", default=0.5)
    alpha = Parameter(description="alpha (power of powerlaw)", default=1.8)
    sil1_amp = Parameter(default=0.05, min=0.01)
    sil1_center = Parameter(default=10.0, bounds=(8.0, 12.0))
    sil1_fwhm = Parameter(default=2.0, bounds=(0.01, 10.0))
    sil1_asym = Parameter(default=0.0)
    sil2_amp = Parameter(default=0.1, min=0.001)
    sil2_center = Parameter(default=20.0, bounds=(16.0, 24.0))
    sil2_fwhm = Parameter(default=3.0, bounds=(0.0, 10.0))
    sil2_asym = Parameter(default=0.0)

    x_range = [1.0 / 40.0, 1.0]

    def evaluate(
        self,
        in_x,
        scale,
        alpha,
        sil1_amp,
        sil1_center,
        sil1_fwhm,
        sil1_asym,
        sil2_amp,
        sil2_center,
        sil2_fwhm,
        sil2_asym,
    ):
        """
        G20 function

        Parameters
        ----------
        in_x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in wavenumbers [1/micron]

        Returns
        -------
        axav: np array (float)
            A(x)/A(V) extinction curve [mag]

        Raises
        ------
        ValueError
           Input x values outside of defined range
        """
        x = _get_x_in_wavenumbers(in_x)

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, self.x_range, "G20")

        # powerlaw
        axav = scale * ((1.0 / x) ** (-1.0 * alpha))

        # silicate feature drudes
        wave = 1 / x
        axav += drude_asym(wave, sil1_amp, sil1_center, sil1_fwhm, sil1_asym)
        axav += drude_asym(wave, sil2_amp, sil2_center, sil2_fwhm, sil2_asym)

        return axav


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="file with the extinction curve to fit")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    parser.add_argument("--path", help="path for the extinction curves")
    args = parser.parse_args()

    if args.path:
        locpath = args.path + "/"
    else:
        locpath = ""

    file = args.file
    ofile = file.replace(".fits", "_POWLAW2DRUDE.fits")

    # read in the observed E(l-V) or A(l)/A(V) extinction curve
    obsext = ExtData(filename=locpath + file)

    # get an observed extinction curve to fit
    (wave, y, y_unc) = obsext.get_fitdata(["BAND", "IRS"])
    # remove units as fitting routines often cannot take numbers with units
    x = wave.to(1.0 / u.micron, equivalencies=u.spectral()).value

    if obsext.type == "elx":
        # determine the initial guess at the A(V) values
        #  just use the average at wavelengths > 5
        #  limit as lambda -> inf, E(lamda-V) -> -A(V)
        (indxs,) = np.where(1.0 / x > 5.0)
        av_guess = -1.0 * np.average(y[indxs])
        if not np.isfinite(av_guess):
            av_guess = 1.0

        g20_init = G20() | AxAvToExv(Av=av_guess)
        g20_asym_init = G20_drude_asym() | AxAvToExv(Av=av_guess)
        # g20_init = G20_x() | AxAvToExv(Av=av_guess)
    elif obsext.type == "alav":
        g20_init = G20()
        g20_asym_init = G20_drude_asym()

    # fit the extinction only using data between 1 and 40 micron
    gvals = (1.0 < 1.0 / x) & (1.0 / x < 40.0)
    fit = LevMarLSQFitter()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        g20_fit = fit(
            g20_init,
            x[gvals],
            y[gvals],
            # weights=1.0 / y_unc[gvals],
            maxiter=10000,
            # epsilon=0.001,
        )

        g20_asym_fit = fit(
            g20_asym_init,
            x[gvals],
            y[gvals],
            # weights=1.0 / y_unc[gvals],
            maxiter=10000,
            # epsilon=0.1,
        )
        print(g20_fit.parameters)
        print(g20_asym_fit.param_names)
        print(g20_asym_fit.parameters)

    # setup the plot
    fontsize = 18
    # fontsize = 10
    font = {"size": fontsize}
    matplotlib.rc("font", **font)
    matplotlib.rc("lines", linewidth=1.5)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    fig, ax = plt.subplots(
        nrows=2, figsize=(12, 8), sharex=True, gridspec_kw={"height_ratios": [5, 1]}
    )

    obsext.plot(ax[0], color="k")
    g20_fit_y = g20_fit(wave[gvals])
    g20_asym_fit_y = g20_asym_fit(wave[gvals])
    # ax.plot(wave[gvals], g20_init(wave[gvals]), "c-", label="Initial Guess")
    ax[0].plot(wave[gvals], g20_fit_y, "g-", label="Best Fit", alpha=0.7)
    # ax[0].plot(wave[gvals], g20_asym_init(wave[gvals]), "b-", label="Asym Init")
    ax[0].plot(wave[gvals], g20_asym_fit_y, "b-", label="Asym Best Fit")
    if obsext.type == "elx":
        ax[0].plot(
            wave[gvals],
            g20_fit.Av_1.value * np.full((len(wave[gvals])), -1.0),
            "g:",
            label="-A(V)",
        )
        ax[0].plot(
            wave[gvals],
            g20_asym_fit.Av_1.value * np.full((len(wave[gvals])), -1.0),
            "b:",
            label="-A(V)",
        )
        ax[0].set_ylabel(r"$E(\lambda - V)$", fontsize=1.3 * fontsize)
    else:
        ax[0].set_ylabel(r"$A(\lambda)/A(V)$", fontsize=1.3 * fontsize)

    ax[1].set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)

    mmy = np.array([min(g20_fit_y), max(g20_fit_y)])
    if obsext.type == "elx":
        mmy[0] = min([mmy[0], -1.0 * g20_fit.Av_1.value, -1.0 * g20_asym_fit.Av_1.value])
    mmd = 0.1 * (mmy[1] - mmy[0])
    ax[0].set_ylim(mmy + np.array([-1.0, 1.0]) * mmd)
    ax[0].set_xlim(1.0, 40.0)
    ax[0].set_xscale("log")
    ax[0].set_title(file)

    ax[0].legend(loc="best")

    # residuals
    ax[1].plot(wave[gvals], np.zeros((len(wave[gvals]))), "k--")
    ax[1].plot(wave[gvals], y[gvals] - g20_fit_y, "g-", alpha=0.7)
    ax[1].plot(wave[gvals], y[gvals] - g20_asym_fit_y, "b-")
    ax[1].set_ylim(np.array([-1.0, 1.0]) * mmd)

    plt.tight_layout()

    # plot or save to a file
    outname = ofile.replace(".fits", "")
    if args.png:
        fig.savefig(outname + ".png")
    elif args.pdf:
        fig.savefig(outname + ".pdf")
    else:
        plt.show()
