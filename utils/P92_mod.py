import numpy as np

from astropy.modeling import Fittable1DModel, Parameter
from dust_extinction.helpers import _get_x_in_wavenumbers, _test_valid_x_range


class P92_mod(Fittable1DModel):
    r"""
    P92 extinction model calculation

    Parameters
    ----------
    BKG_amp : float
      background term amplitude
    BKG_lambda : float
      background term central wavelength
    BKG_b : float
      background term b coefficient
    BKG_n : float
      background term n coefficient [FIXED at n = 2]

    FUV_amp : float
      far-ultraviolet term amplitude
    FUV_lambda : float
      far-ultraviolet term central wavelength
    FUV_b : float
      far-ultraviolet term b coefficent
    FUV_n : float
      far-ultraviolet term n coefficient

    NUV_amp : float
      near-ultraviolet (2175 A) term amplitude
    NUV_lambda : float
      near-ultraviolet (2175 A) term central wavelength
    NUV_b : float
      near-ultraviolet (2175 A) term b coefficent
    NUV_n : float
      near-ultraviolet (2175 A) term n coefficient [FIXED at n = 2]

    SIL1_amp : float
      1st silicate feature (~10 micron) term amplitude
    SIL1_lambda : float
      1st silicate feature (~10 micron) term central wavelength
    SIL1_b : float
      1st silicate feature (~10 micron) term b coefficent
    SIL1_n : float
      1st silicate feature (~10 micron) term n coefficient [FIXED at n = 2]

    SIL2_amp : float
      2nd silicate feature (~18 micron) term amplitude
    SIL2_lambda : float
      2nd silicate feature (~18 micron) term central wavelength
    SIL2_b : float
      2nd silicate feature (~18 micron) term b coefficient
    SIL2_n : float
      2nd silicate feature (~18 micron) term n coefficient [FIXED at n = 2]

    FIR_amp : float
      far-infrared term amplitude
    FIR_lambda : float
      far-infrared term central wavelength
    FIR_b : float
      far-infrared term b coefficent
    FIR_n : float
      far-infrared term n coefficient [FIXED at n = 2]

    Notes
    -----
    P92 extinction model

    From Pei (1992)

    Applicable from the extreme UV to far-IR

    Example showing a P92 curve with components identified.

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        import astropy.units as u

        from dust_extinction.shapes import P92

        fig, ax = plt.subplots()

        # generate the curves and plot them
        lam = np.logspace(-3.0, 3.0, num=1000)
        x = (1.0/lam)/u.micron

        ext_model = P92()
        ax.plot(1/x,ext_model(x),label='total')

        ext_model = P92(FUV_amp=0., NUV_amp=0.0,
                        SIL1_amp=0.0, SIL2_amp=0.0, FIR_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG only')

        ext_model = P92(NUV_amp=0.0,
                        SIL1_amp=0.0, SIL2_amp=0.0, FIR_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG+FUV only')

        ext_model = P92(FUV_amp=0.,
                        SIL1_amp=0.0, SIL2_amp=0.0, FIR_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG+NUV only')

        ext_model = P92(FUV_amp=0., NUV_amp=0.0,
                        SIL2_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG+FIR+SIL1 only')

        ext_model = P92(FUV_amp=0., NUV_amp=0.0,
                        SIL1_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG+FIR+SIL2 only')

        ext_model = P92(FUV_amp=0., NUV_amp=0.0,
                        SIL1_amp=0.0, SIL2_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG+FIR only')

        # Milky Way observed extinction as tabulated by Pei (1992)
        MW_x = [0.21, 0.29, 0.45, 0.61, 0.80, 1.11, 1.43, 1.82,
                2.27, 2.50, 2.91, 3.65, 4.00, 4.17, 4.35, 4.57, 4.76,
                5.00, 5.26, 5.56, 5.88, 6.25, 6.71, 7.18, 7.60,
                8.00, 8.50, 9.00, 9.50, 10.00]
        MW_x = np.array(MW_x)
        MW_exvebv = [-3.02, -2.91, -2.76, -2.58, -2.23, -1.60, -0.78, 0.00,
                     1.00, 1.30, 1.80, 3.10, 4.19, 4.90, 5.77, 6.57, 6.23,
                     5.52, 4.90, 4.65, 4.60, 4.73, 4.99, 5.36, 5.91,
                     6.55, 7.45, 8.45, 9.80, 11.30]
        MW_exvebv = np.array(MW_exvebv)
        Rv = 3.08
        MW_axav = MW_exvebv/Rv + 1.0
        ax.plot(1./MW_x, MW_axav, 'o', label='MW Observed')

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_ylim(1e-3,10.)

        ax.set_xlabel(r'$\lambda$ [$\mu$m]')
        ax.set_ylabel(r'$A(x)/A(V)$')

        ax.legend(loc='best')
        plt.show()
    """

    # ['FUV_width', 'NUV_width', 'SIL1_width', 'SIL2_width', 'FIR_width']
    # [0.14696938456699066, 0.051719412985268637, 3.7049250106326643, 8.3409029232536902, 35.355339059327378]

    # inputs = ("x",)
    # outputs = ("axav",)

    # constant for conversion from Ax/Ab to (more standard) Ax/Av
    AbAv = 1.0 / 3.08 + 1.0

    BKG_amp = Parameter(
        description="BKG term: amplitude", default=165.0 * AbAv, min=0.0
    )
    BKG_lambda = Parameter(description="BKG term: center wavelength", default=0.047)
    BKG_width = Parameter(description="BKG term: width", default=0.452, min=0.0)

    FUV_amp = Parameter(description="FUV term: amplitude", default=14.0 * AbAv, min=0.0)
    FUV_lambda = Parameter(
        description="FUV term: center wavelength", default=0.07, bounds=(0.06, 0.08)
    )
    FUV_b = Parameter(description="FUV term: b coefficient", default=4.0)
    FUV_n = Parameter(description="FUV term: n coefficient", default=6.5)

    NUV_amp = Parameter(
        description="NUV term: amplitude", default=0.045 * AbAv, min=0.0
    )
    NUV_lambda = Parameter(
        description="NUV term: center wavelength",
        default=0.2175,
        bounds=(0.2100, 0.2250),
    )
    NUV_width = Parameter(description="NUV term: width", default=0.05, min=0.0)

    SIL1_amp = Parameter(
        description="SIL1 term: amplitude", default=0.002 * AbAv, min=0.0
    )
    SIL1_lambda = Parameter(
        description="SIL1 term: center wavelength", default=9.7, bounds=(7.0, 13.0)
    )
    SIL1_width = Parameter(
        description="SIL1 term: width", default=2.0, bounds=(0.0, 10.0)
    )

    SIL2_amp = Parameter(
        description="SIL2 term: amplitude", default=0.002 * AbAv, min=0.00
    )
    SIL2_lambda = Parameter(
        description="SIL2 term: center wavelength", default=18.0, bounds=(15.0, 21.0)
    )
    SIL2_width = Parameter(description="SIL2 term: width", default=8.3, min=0.0)

    FIR_amp = Parameter(
        description="FIR term: amplitude", default=0.012 * AbAv, min=0.0
    )
    FIR_lambda = Parameter(
        description="FIR term: center wavelength", default=25.0, bounds=(20.0, 30.0)
    )
    FIR_width = Parameter(description="FIR term: width", default=35.0)

    x_range = [1.0 / 1e3, 1.0 / 1e-3]

    @staticmethod
    def _p92_single_term(in_lambda, amplitude, cen_wave, b, n):
        r"""
        Function for calculating a single P92 term

        .. math::

           \frac{a}{(\lambda/cen_wave)^n + (cen_wave/\lambda)^n + b}

        when n = 2, this term is equivalent to a Drude profile

        Parameters
        ----------
        in_lambda: vector of floats
           wavelengths in same units as cen_wave

        amplitude: float
           amplitude

        cen_wave: flaot
           central wavelength

        b : float
           b coefficient

        n : float
           n coefficient
        """
        l_norm = in_lambda / cen_wave

        return amplitude / (np.power(l_norm, n) + np.power(l_norm, -1 * n) + b)

    def evaluate(
        self,
        in_x,
        BKG_amp,
        BKG_lambda,
        BKG_width,
        FUV_amp,
        FUV_lambda,
        FUV_b,
        FUV_n,
        NUV_amp,
        NUV_lambda,
        NUV_width,
        SIL1_amp,
        SIL1_lambda,
        SIL1_width,
        SIL2_amp,
        SIL2_lambda,
        SIL2_width,
        FIR_amp,
        FIR_lambda,
        FIR_width,
    ):
        """
        P92 function

        Parameters
        ----------
        in_x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in wavenumbers [1/micron]

           internally wavenumbers are used

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
        _test_valid_x_range(x, self.x_range, "P92_mod")

        # compute b from lambda and width
        BKG_b = np.power((BKG_width / BKG_lambda), 2.0) - 2.0
        NUV_b = np.power((NUV_width / NUV_lambda), 2.0) - 2.0
        SIL1_b = np.power((SIL1_width / SIL1_lambda), 2.0) - 2.0
        SIL2_b = np.power((SIL2_width / SIL2_lambda), 2.0) - 2.0
        FIR_b = np.power((FIR_width / FIR_lambda), 2.0) - 2.0

        # calculate the terms
        lam = 1.0 / x
        axav = (
            self._p92_single_term(lam, BKG_amp, BKG_lambda, BKG_b, 2.0)
            + self._p92_single_term(lam, FUV_amp, FUV_lambda, FUV_b, FUV_n)
            + self._p92_single_term(lam, NUV_amp, NUV_lambda, NUV_b, 2.0)
            + self._p92_single_term(lam, SIL1_amp, SIL1_lambda, SIL1_b, 2.0)
            + self._p92_single_term(lam, SIL2_amp, SIL2_lambda, SIL2_b, 2.0)
            + self._p92_single_term(lam, FIR_amp, FIR_lambda, FIR_b, 2.0)
        )

        # return A(x)/A(V)
        return axav

    # use numerical derivaties (need to add analytic)
    fit_deriv = None
