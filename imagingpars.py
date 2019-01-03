#!/usr/bin/python
#-*- coding: latin-1 -*-
"""Contains imaging parameters like effective pixelsize, absorption
coefficient, ..."""
class ImagingPars(object):
    """Base class for parameters of imaging system.

    @cvar description: descriptive name of settings, used for GUI selection.
    @cvar pixelsize: Size of area (in µm) which is imaged to one pixel of the cam.
    @cvar sigma0: cross section for light absorption
    """
    description = None
    pixelsize = 1
    sigma0 = 1.5/3.14*(589e-9)**2
    mass = 23.0 * 1.66054e-27
    expansion_time = 0
    ODmax = 0 #maximum optical density

    def __str__(self):
        s = "%s, t_exp = %.1fms, OD_max = %.1f"%(
            self.description,
            self.expansion_time,
            self.ODmax)
        return s

class ImagingParsCMOS(ImagingPars):
    description = "cmos"
    pixelsize = 6.50/2.0 #0.972  #pixelsize in µm / imaging magnification
    palette = "gist_stern"
        
class ImagingParsVertical(ImagingPars):
    description = "vertical"
    pixelsize = 4.40*2.5 #0.972  #pixelsize in µm / imaging magnification
    palette = "gist_stern"

class ImagingParsHorizontal(ImagingPars):
    description = "horizontal"
    pixelsize = 4.40*1.0 #1.002  #pixelsize in µm / imaging magnification
    palette = "gist_stern"

class ImagingParsHorizontalDemag3(ImagingPars):
    description = "horizontal demagnified 3"
    pixelsize = 4.40*3 #1.002  #pixelsize in µm / imaging magnification
    palette = "gist_stern"

class ImagingParsHorizontalOld(ImagingPars):
    description = "horizontal"
    pixelsize = 4.40*2.6 #1.002  #pixelsize in µm / imaging magnification
    palette = "gist_stern"

class ImagingParsDark(ImagingPars):
    description = "dark"
    pixelsize = 4.40/2.0 #1.002  #pixelsize in µm / imaging magnification
    palette = "gist_stern"

class ImagingParsPixelSize(ImagingPars):
    description = "PixelSize"
    pixelsize = 4.40 #1.002  #pixelsize in µm / imaging magnification
    palette = "gist_stern"

