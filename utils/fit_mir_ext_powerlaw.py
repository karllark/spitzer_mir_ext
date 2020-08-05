import matplotlib.pyplot as plt
import matplotlib

import numpy as np
import argparse

import astropy.units as u
from astropy.modeling import Fittable1DModel, Parameter
from astropy.modeling.fitting import LevMarLSQFitter

from measure_extinction.extdata import ExtData
from dust_extinction.shapes import P92
from dust_extinction.helpers import _get_x_in_wavenumbers, _test_valid_x_range
from dust_extinction.conversions import AxAvToExv

from models_mcmc_extension import EmceeFitter

# from pahfit.component_models import Drude1D


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

    inputs = ("x",)
    outputs = ("axav",)

    scale = Parameter(description="amplitude", default=0.5)
    alpha = Parameter(description="alpha (power of powerlaw)", default=1.8)
    sil1_amp = Parameter(default=0.05, min=0.01)
    sil1_center = Parameter(default=10.0, bounds=(8.0, 12.0))
    sil1_fwhm = Parameter(default=0.0, bounds=(0.0, 2.0))
    sil2_amp = Parameter(default=0.1, min=0.001)
    sil2_center = Parameter(default=20.0, bounds=(16.0, 24.0))
    sil2_fwhm = Parameter(default=3.0, bounds=(0.0, 5.0))
    fir_amp = Parameter(default=0.0, min=0.00, fixed=True)
    fir_center = Parameter(default=25.0, bounds=(20.0, 30.0))
    fir_fwhm = Parameter(default=10.0, min=5.0, max=20.0)

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
        fir_amp,
        fir_center,
        fir_fwhm,
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
        axav += (
            sil1_amp
            * ((sil1_fwhm / sil1_center) ** 2)
            / (
                (wave / sil1_center - sil1_center / wave) ** 2
                + (sil1_fwhm / sil1_center) ** 2
            )
        )

        axav += (
            sil2_amp
            * ((sil2_fwhm / sil2_center) ** 2)
            / (
                (wave / sil2_center - sil2_center / wave) ** 2
                + (sil2_fwhm / sil2_center) ** 2
            )
        )

        axav += (
            fir_amp
            * ((fir_fwhm / fir_center) ** 2)
            / (
                (wave / fir_center - fir_center / wave) ** 2
                + (fir_fwhm / fir_center) ** 2
            )
        )

        return axav


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="file with the extinction curve to fit")
    parser.add_argument(
        "-p", "--png", help="save figure as a png file", action="store_true"
    )
    parser.add_argument("--path", help="path for the extinction curves")
    args = parser.parse_args()

    if args.path:
        locpath = args.path + "/"
    else:
        locpath = ""

    # read in the observed E(l-V) extinction curve
    obsext = ExtData(filename=locpath + args.file)

    # get an observed extinction curve to fit
    (wave, y, y_unc) = obsext.get_fitdata(["BAND", "IRS"])
    # remove units as fitting routines often cannot take numbers with units
    x = wave.to(1.0 / u.micron, equivalencies=u.spectral()).value

    # determine the initial guess at the A(V) values
    #  just use the average at wavelengths > 5
    #  limit as lambda -> inf, E(lamda-V) -> -A(V)
    indxs, = np.where(1.0 / x > 5.0)
    av_guess = -1.0 * np.average(y[indxs])
    if not np.isfinite(av_guess):
        av_guess = 1.0

    # initialize the model
    g20_init = G20() | AxAvToExv(Av=av_guess)

    # fit the extinction only using data between 1 and 40 micron
    gvals = (1.0 < 1.0 / x) & (1.0 / x < 40.0)
    fit = LevMarLSQFitter()
    g20_fit = fit(g20_init, x[gvals], y[gvals])  # , weights=1.0 / y_unc[gvals])

    print(g20_fit)

    # setup the plot
    fontsize = 18
    fontsize = 10
    font = {"size": fontsize}
    matplotlib.rc("font", **font)
    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    fig, ax = plt.subplots(figsize=(12, 8))

    obsext.plot(ax)
    g20_fit_y = g20_fit(wave[gvals])
    # ax.plot(wave[gvals], g20_init(wave[gvals]), "c-", label="Initial Guess")
    ax.plot(wave[gvals], g20_fit_y, "r-", label="Best Fit")
    ax.plot(
        wave[gvals],
        g20_fit.Av_1.value * np.full((len(wave[gvals])), -1.0),
        "-",
        label="-A(V)",
    )

    ax.set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)
    ax.set_ylabel(r"$E(\lambda - V)$", fontsize=1.3 * fontsize)

    ax.set_ylim(min(g20_fit_y) - 0.5, max(g20_fit_y) + 0.5)
    ax.set_xlim(1.0, 40.0)
    ax.set_xscale("log")

    ax.legend(loc="best")

    plt.tight_layout()

    plt.show()
