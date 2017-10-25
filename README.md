# mfdrouting - Multiple Flow Direction (MFD) flow routing

This Python 3 module implements the Multiple Flow Direction (MFD) flow routing
as proposed by [Freeman (1991)][id1]. As suggested by Freeman the free
parameter p, the exponent that scales the fractions of flows, is tuned to the
divergent flows of cones with p = 1.1.  It is written in C and based on the
code by [Pelletier (2008)][id2].

## Install

	git clone https://github.com/UP-RS-ESP/mfdrouting.git
	cd mfdrouting
	sudo python3 setup.py install

## Usage

    import numpy as np
    import gdal
    import mfdrouting
    from matplotlib import pyplot as pl

[id1]: http://dx.doi.org/10.1016/0098-3004(91)90048-i "Calculating catchment area with divergent flow based on a regular grid. T. Graham Freeman, Computers & Geosciences (1991)."

[id2]: http://dx.doi.org/10.1594/IEDA/100145 "MFDrouting, version 0.1. J. D. Pelletier (2008)."
