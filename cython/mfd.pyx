import numpy as np

__version__ = 0.1
__author__ = 'Aljoscha Rheinwalt <aljoscha.rheinwalt@uni-potsdam.de>'

cdef extern from "cmfd.c":
    void mfdtda(double *a, const double *z,
        const int xlen, const int ylen,
        const double w, const double ppexp);

def sca(const double[:,:] dem, const double cellwidth, const double pexp = 1.1):
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
    cdef double[:, :] mv
    cdef unsigned int n, m

    # NaN value mask
    mask = dem == np.nan
    emin = dem[~mask].min()
    dem -= emin
    dem[mask] = 0

    # Numpy output array
    n, m = dem.shape[0], dem.shape[1]
    out = np.zeros((n, m))
    mv = out

    # C call
    mfdtda(&mv[0,0], &dem[0,0],
            m, n, cellwidth, pexp)

    # TDA -> SCA
    out /= cellwidth

    # fill in NaNs
    out[mask] = np.nan

    return out

