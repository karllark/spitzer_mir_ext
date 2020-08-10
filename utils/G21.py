import numpy as np

from astropy.modeling import Fittable1DModel, Parameter

from dust_extinction.helpers import _get_x_in_wavenumbers, _test_valid_x_range


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

    scale = Parameter(description="amplitude", default=0.5, bounds=(0.0, 1.0))
    alpha = Parameter(
        description="alpha (power of powerlaw)", default=1.8, bounds=(0.5, 5.0)
    )
    sil1_amp = Parameter(default=0.07, bounds=(0.001, 0.3))
    sil1_center = Parameter(default=10.0, bounds=(8.0, 12.0))
    sil1_fwhm = Parameter(default=2.5, bounds=(1.0, 10.0))
    sil1_asym = Parameter(default=-0.3, bounds=(-2.0, 2.0))
    sil2_amp = Parameter(default=0.03, bounds=(0.001, 0.3))
    sil2_center = Parameter(default=20., bounds=(16.0, 24.0))
    sil2_fwhm = Parameter(default=10.0, bounds=(3.0, 15.0))
    sil2_asym = Parameter(default=0.0, bounds=(-2.0, 2.0))

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
