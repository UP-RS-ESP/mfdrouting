import numpy as np

__version__ = 0.1
__author__ = 'Aljoscha Rheinwalt <aljoscha.rheinwalt@uni-potsdam.de>'

cdef extern from "cmfd.c":
    void mfdtda(double *a, const double *z,
        const int xlen, const int ylen,
        const double w, const double ppexp);

def sca(dem, const double cellwidth, const double pexp = 1.1):
    """
    sca(dem, cellwidth, exponent = 1.1)

    Returns the specific catchment area (SCA) for a given
    digital elevation model (DEM) using the multiple flow
    direction (MFD) algorithm as proposed by Freeman (1991).

    https://doi.org/10.1016/0098-3004(91)90048-I

    This code is based on the implementation by Pelletier (2008).
    
    http://dx.doi.org/10.1594/IEDA/100145

    Parameters
    ----------
    dem : ndarray
        A 2D Numpy array containing the elevations of the topography (DEM).
    cellwidth : float
        The cell width of the DEM in meters.
    exponent : float, optional
        The exponent determining the flow divergence in relation to the slope.
        See http://dx.doi.org/10.1016/0098-3004(91)90048-i.

    Returns
    -------
    a : ndarray
        The specific catchment area Numpy array on the same grid.
    """
    cdef double[:, :] dv, av
    cdef unsigned int n, m

    # NaN value mask
    ndem = np.asarray(dem)
    mask = np.isnan(ndem)
    emin = ndem[~mask].min()
    ndem -= emin
    ndem += 1
    ndem[mask] = 0
    dv = ndem

    # Numpy output array
    n, m = ndem.shape[0], ndem.shape[1]
    out = np.zeros((n, m))
    av = out

    # C call
    mfdtda(&av[0,0], &dv[0,0],
            m, n, cellwidth, pexp)

    # recover dem nans
    ndem += emin
    ndem -= 1
    ndem[mask] = np.nan

    # TDA -> SCA
    out /= cellwidth

    # fill in NaNs
    out[mask] = np.nan

    return out

