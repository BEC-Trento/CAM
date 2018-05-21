#!/usr/bin/python
#-*- coding: latin-1 -*-
"""Main file for camera program."""

from __future__ import with_statement

import matplotlib
#matplotib.use('WxAgg')

import os, time, shutil, pickle, re, copy
import pylab
import numpy
import scipy.ndimage.filters
import scipy.signal
#import math

from numpy import ma

#gui
import wx, wx.aui, wx.grid
import wx.lib.delayedresult as delayedresult

from matplotlib.widgets import Button#, Cursor 
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import StatusBarWx

from matplotlib.backends.backend_wx import NavigationToolbar2Wx
#from toolbar import NavigationToolbar2Wx

from matplotlib.figure import Figure

#imports
from readsis import loadimg, loadimg3, read, write_raw_image
#loadimg = loadimg3

import settings
import fitting
import FitResultTableGrid
import DataPanel
import filewatch
import ding
from roi import ROI
import imagingpars
import ImageTree
import ImagePanel
from custom_events import *
#from profiling import Tic
import mean_sis


#force reload of modules (for development)
reload(FitResultTableGrid)
reload(settings)
reload(fitting)
reload(DataPanel)
reload(filewatch)
reload(ding)
reload(ImagePanel)
reload(ImageTree)


class CamCursor(object):
    """
    A horizontal and vertical line span the axes that and move with
    the pointer.  You can turn off the hline or vline spectively with
    the attributes

      horizOn =True|False: controls visibility of the horizontal line
      vertOn =True|False: controls visibility of the horizontal line

    And the visibility of the cursor itself with visible attribute
    """
    def __init__(self, ax, useblit=False, **lineprops):
        """
        Add a cursor to ax.  If useblit=True, use the backend
        dependent blitting features for faster updates (GTKAgg only
        now).  lineprops is a dictionary of line properties.  See
        examples/widgets/cursor.py.
        """
        self.ax = ax
        self.canvas = ax.figure.canvas

        self.canvas.mpl_connect('motion_notify_event', self.onmove)
        self.canvas.mpl_connect('draw_event', self.ondraw)

        self._visible = True
        self.horizOn = True
        self.vertOn = True
        self.useblit = useblit

        self.lineh = ax.axhline(ax.get_ybound()[0], visible=False, **lineprops)
        self.linev = ax.axvline(ax.get_xbound()[0], visible=False, **lineprops)

        self.background = None
        self.needclear = False
        
    def save_background(self):
        #print "save background"
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)

    def set_visible(self, visible):
        if visible == self._visible:
            return
        
        if visible == True:
            self.save_background()
        else:
            self.linev.set_visible(visible and self.vertOn)
            self.lineh.set_visible(visible and self.horizOn)

            pass
            #self._update()

        self._visible = visible
        

    def get_visible(self):
        return self._visible

    visible = property(get_visible, set_visible)

    def ondraw(self, event):
        #print "draw event"
        if self.visible:
            self.clear(event)
        else:
            pass
        #self.clear(event)

    def clear(self, event):
        'clear the cursor'
        #print "Cursor/clear"
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.linev.set_visible(False)
        self.lineh.set_visible(False)

    def onmove(self, event):
        'on mouse motion draw the cursor if visible'
        if not self.visible: return #moved to beginning, this is the most important improvement

        if event.inaxes != self.ax:
            self.linev.set_visible(False)
            self.lineh.set_visible(False)

            if self.needclear:
                self._update() #instead of draw
                self.needclear = False
            return
        self.needclear = True
        self.linev.set_xdata((event.xdata, event.xdata))
        self.lineh.set_ydata((event.ydata, event.ydata))
        self.linev.set_visible(self.visible and self.vertOn)
        self.lineh.set_visible(self.visible and self.horizOn)

        self._update()


    def _update(self):
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            else:
                #print "nothing to restore!"
                pass
            
            self.ax.draw_artist(self.linev)
            self.ax.draw_artist(self.lineh)
            self.canvas.blit(self.ax.bbox)
        else:

            self.canvas.draw_idle()

        return False



class Selector(object):
    "Base class for X/YSelector"
    
    def __init__(self, axis, num, tolerance=5, marker='+', color='y', alpha=0.5):
        self.axis = axis
        self.num = num
        self.tol = tolerance
        self.marker = marker
        self.color = color
        self.alpha = alpha

        self.cursor = CamCursor(self.axis, useblit=True)
        self.cursor.lineh.set_color(self.color)
        self.cursor.linev.set_color(self.color)
        self.init_cursor()
                             
        self.activated = False
        self.cursor.visible = False
        
    def init_cursor(self):
        pass
    
    def picker(self, artist, mousevent):
        for overaxis in [mousevent.inaxes]:
            if overaxis is self.axis:
                if not self.activated:
                    distance = self.distance(artist, mousevent)
                    if distance < self.tol: #hitted line
                        self.activated = True
                        #artist.set_alpha(1)
                        self.axis.figure.canvas.draw()
                        self.cursor.visible = True
                        return False, {}
                    
                else: #finish selection
                    self.activated = False
                    self.cursor.visible = False
                    newpos = self.update_position(artist, mousevent)
                    #artist.set_alpha(self.alpha)
                    overaxis.figure.canvas.draw() 
                    newpos['num'] = self.num
                    return True, newpos

        return False, {}

    def update_position(self, artist, mousevent):
        """set position of artist to position given by mouse and return
        new position as dict for use in picker callback"""
        pass

    def update_from_roi(self, roi):
        "set position from roi"
        pass

class YSelector(Selector):
    def __init__(self, axis, pos, num, **kwargs):
        Selector.__init__(self, axis, num, **kwargs)
        self.line = self.axis.axhline(pos,
                                      alpha=self.alpha,
                                      color=self.color,
                                      picker=self.picker)

    def init_cursor(self):
        self.cursor.vertOn = False

    def update_position(self, artist, mousevent):
        pos = round(mousevent.ydata)
        artist.set_ydata([pos] * 2)
        return dict(roiy=pos)

    def update_from_roi(self, roi):
        self.line.set_ydata([roi.gety(self.num)] * 2)
    
    def distance(self, artist, mousevent):
        linepos = artist.get_ydata()[0]
        #sx, sy     = self.axis.transData.xy_tup((0, linepos))
        sx, sy = self.axis.transData.transform((0, linepos))
        return abs(sy - mousevent.y)
            
        
class XSelector(Selector):
    def __init__(self, axis, pos, num, **kwargs):
        Selector.__init__(self, axis, num, **kwargs)
        self.line = self.axis.axvline(pos,
                                      alpha=self.alpha,
                                      color=self.color,
                                      picker=self.picker)
    def init_cursor(self):
        self.cursor.horizOn = False

    def update_position(self, artist, mousevent):
        pos = round(mousevent.xdata)
        artist.set_xdata([pos] * 2)
        return dict(roix=pos)

    def update_from_roi(self, roi):
        self.line.set_xdata([roi.getx(self.num)] * 2)

    def distance(self, artist, mousevent):
        linepos = artist.get_xdata()[0]
        #sx, sy   = self.axis.transData.xy_tup((linepos, 0))
        sx, sy = self.axis.transData.transform((linepos, 0))
        return abs(sx - mousevent.x)


class XYSelector(Selector):
    def __init__(self, axis, pos, num, **kwargs):
        Selector.__init__(self, axis, num, **kwargs)
        self.marker = self.axis.plot([pos[0]], [pos[1]],
                                     marker=self.marker,
                                     markersize=9,
                                     markeredgewidth=1.0,
                                     alpha=self.alpha,
                                     color=self.color,
                                     picker=self.picker)
        self.marker = self.marker[0]
        self.lineV = self.axis.axvline(pos[0],
                                        color=self.color,
                                        alpha=self.alpha)
        self.lineH = self.axis.axhline(pos[1],
                                        color=self.color,
                                        alpha=self.alpha)
        self.set_visible(False)
        
    def update_position(self, artist, mousevent):
        posx = round(mousevent.xdata)
        posy = round(mousevent.ydata)
        self.set_position((posx, posy))
        return dict(x=posx, y=posy)

    def distance(self, artist, mousevent):
        posx = artist.get_xdata()[0]
        posy = artist.get_ydata()[0]
        sx, sy = self.axis.transData.transform((posx, posy))
        return numpy.sqrt((sx - mousevent.x) ** 2 + (sy - mousevent.y) ** 2)

    def set_position(self, point):
        posx, posy = point
        self.marker.set_xdata([posx])
        self.marker.set_ydata([posy])
        self.lineV.set_xdata([posx] * 2)
        self.lineH.set_ydata([posy] * 2)

    def get_position(self):
        posx = self.marker.get_xdata()[0]
        posy = self.marker.get_ydata()[0]
        return (posx, posy)

    position = property(get_position, set_position)

    def set_visible(self, visible):
        self.cross_visible = visible
        self.lineV.set_visible(self.cross_visible)
        self.lineH.set_visible(self.cross_visible)

    def get_visible(self):
        visible = self.lineV._visible
        return visible


class CamSplashScreen(wx.SplashScreen):
    def __init__(self):
        bitmap = wx.Image(os.path.join(settings.bitmappath, 'cam_splash.png'),
                          wx.BITMAP_TYPE_PNG).ConvertToBitmap()
        wx.SplashScreen.__init__(self,
                                 bitmap,
                                 wx.SPLASH_CENTRE_ON_SCREEN | wx.SPLASH_TIMEOUT,
                                 1000, None, - 1)

class ImgPanel(wx.Panel):
    """Class for displaying images and handling fit results, realized
    as wxPanel"""

    ID_Autoscale = wx.NewId()
    
    def __init__(self, parent, panel_id, img, region_of_interest, fitroutine, palette=pylab.cm.jet):

        #Initialize wx stuff 
        wx.Panel.__init__(self,
                          parent,
                          size=(550, 350))

        self.fig = pylab.Figure()
        self.canvas = FigureCanvas(self, -1, self.fig)
        
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

        #setup (and modify) toolbar
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.AddLabelTool(self.ID_Autoscale,
                                  'Autoscale',
                                  wx.Bitmap(os.path.join(settings.bitmappath,
                                                         'autoscale.png'),
                                            wx.BITMAP_TYPE_PNG),
                                  shortHelp = 'Autoscale',
                                  longHelp = 'automatic scaling')
        wx.EVT_TOOL(self, self.ID_Autoscale, self.OnAutoscale)
        self.toolbar.Realize()
        tw, th = self.toolbar.GetSizeTuple()
        fw, fh = self.canvas.GetSizeTuple()
        self.toolbar.SetSize(wx.Size(fw, th))
        self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.toolbar.update()

        #store input data
        self.id = panel_id
        self.rawimg = img
        self.img = img #for the moment
        self.fit = fitroutine
        self.roi = region_of_interest
        
        #initialize attributes
        self.x = range(img.shape[0])
        self.y = range(img.shape[1])
        self.fitpars = fitting.FitParsNoFit()
        self.compensate_offset = False
        self.target_results = None
        self.difference=False
        self.filter=False
        self.filter_scalar=1.0
        self.filter_sigma=1.0
        
        self.reference_diff=False
        self.reference_img=None
        
        self.rotation_angle=0.0
        self.rotation_active=False
        self.bec_residuals=None
        
        #create axes
        #self.axtop   = self.fig.add_axes([0.2, 0.925, 0.6, 0.05]) #top (reload)
        self.aximg = self.fig.add_axes([0.2, 0.3, 0.6, 0.65])    #image
        self.axvprof = self.fig.add_axes([0.05, 0.3, 0.13, 0.65], sharey=self.aximg) #vertical profile
        self.axhprof = self.fig.add_axes([0.2, 0.05, 0.6, 0.23], sharex=self.aximg) #horizontal profile
        #self.axclrb = self.fig.add_axes([0.82, 0.3, 0.05, 0.65])  #colorbar
        #self.axtxt = self.fig.add_axes([0.82, 0.05, 0.15, 0.23])#text
        self.axclrb = self.fig.add_axes([0.82, 0.05, 0.05, 0.23]) #colorbar
        self.axtxt = self.fig.add_axes([0.82, 0.3, 0.05, 0.65])  #text


        #text area
        self.axtxt.set_axis_off()
        self.htxt = self.axtxt.text(0, 0, "Hi",
                                    axes=self.axtxt,
                                    multialignment='left',
                                    transform=self.axtxt.transAxes,
                                    family='monospace',
                                    )
        
        #initialize image area
        
        self.himg = self.aximg.imshow(img, cmap=palette, vmin= - 0.05, vmax=1.5,
                                      origin="upper",
                                      interpolation='bicubic' #'nearest', 'spline36'
                                      )
        
        self.aximg.set_axis_off()
        self.hclrbr = pylab.colorbar(self.himg, self.axclrb)


        #initialize profiles
        self.create_profiles()

        #initialize contour handles
        self.hcontours = []

        #create borders ROI
        self.hrx = []

        self.hrx.append(
            XSelector(
                self.aximg,
                self.roi.xmin,
                num=0))

        self.hrx.append(
            XSelector(
                self.aximg,
                self.roi.xmax,
                num=1))
        
        self.hry = []
        self.hry.append(
            YSelector(
                self.aximg,
                self.roi.ymin,
                num=0
            ))

        self.hry.append(
            YSelector(
                self.aximg,
                self.roi.ymax,
                num=1
            ))
        
        #initialize markers
        self.markers = []
        self.markers.append(
            XYSelector(self.aximg, (440, 295), num=0, marker='x',color='w', alpha=0.8))
        self.markers.append(
            XYSelector(self.aximg, (485, 295), num=1, color='w', alpha=0.8))
        self.markers.append(
            XYSelector(self.aximg, (440, 340), num=2, color='g', alpha=0.8))


        #connect event handler for manipulation of ROI-borders
        self.fig.canvas.mpl_connect('pick_event', self.onpick)

        #set initialize limits for image _and_ profiles
        self.center_roi()

        #TODO: queak
        self.toolbar.push_current()

        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        
        self.fitjobID = 0
        self.Bind(EVT_FIT_COMPLETED, self.OnFitCompleted)
        self.update()

        self.axhprof.autoscale_view(scalex = False, scaley = True)
        self.axvprof.autoscale_view(scalex = True, scaley = False)

    def OnEraseBackground(self, event):
        event.Skip()
        
    def create_profiles(self):
        r = self.roi
        imgroi = self.img[r.y, r.x]
        self.hprof = imgroi.sum(0)
        self.vprof = imgroi.sum(1)

        #self.hprof = ma.array(self.hprof, mask = ~numpy.isfinite(self.hprof))
        #self.vprof = ma.array(self.vprof, mask = ~numpy.isfinite(self.vprof))
        #takex = numpy.isfinite(self.hprof)
        #takey = numpy.isfinite(self.vprof)

        x = r.xrange_clipped(self.img)
        y = r.yrange_clipped(self.img)

        self.hhprof = self.axhprof.plot(x, self.hprof, zorder=10,
                                        #drawstyle = 'steps-mid'
                                        )[0]
        self.hvprof = self.axvprof.plot(self.vprof, y, zorder=10,
                                        #drawstyle = 'steps-mid'
                                        )[0]
        
        self.hhproffit = []
        self.hvproffit = []
        for color, zorder in zip(('r', 'm', 'k:'), (5, 4, 0)):
            self.hhproffit += self.axhprof.plot([], [], color, zorder=zorder)
            self.hvproffit += self.axvprof.plot([], [], color, zorder=zorder)

    def update_profiles(self):
        r = self.roi
        imgroi = self.img[r.y, r.x]

        if self.compensate_offset:
            cimg = imgroi - self.img_background
        else:
            cimg = imgroi

        self.hprof = cimg.filled().sum(0)
        self.vprof = cimg.filled().sum(1)
        
        self.hhprof.set_data(r.xrange_clipped(self.img), self.hprof)
        self.hvprof.set_data(self.vprof, r.yrange_clipped(self.img))
        
        
        
        
        
        
    def mpl_imsave(self,fname, arr, vmin=None, vmax=None, cmap=None, format=None, origin=None, dpi=100):
        """
        Save an array as in image file.
        
        The output formats available depend on the backend being used.
        
        Arguments:
        *fname*:
        A string containing a path to a filename, or a Python file-like object.
        If *format* is *None* and *fname* is a string, the output
        format is deduced from the extension of the filename.
        *arr*:
        An MxN (luminance), MxNx3 (RGB) or MxNx4 (RGBA) array.
        Keyword arguments:
        *vmin*/*vmax*: [ None | scalar ]
        *vmin* and *vmax* set the color scaling for the image by fixing the
        values that map to the colormap color limits. If either *vmin*
        or *vmax* is None, that limit is determined from the *arr*
        min/max value.
        *cmap*:
        cmap is a colors.Colormap instance, eg cm.jet.
        If None, default to the rc image.cmap value.
        *format*:
        One of the file extensions supported by the active
        backend. Most backends support png, pdf, ps, eps and svg.
        *origin*
        [ 'upper' | 'lower' ] Indicates where the [0,0] index of
        the array is in the upper left or lower left corner of
        the axes. Defaults to the rc image.origin value.
        *dpi*
        The DPI to store in the metadata of the file. This does not affect the
        resolution of the output image.
        """
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from matplotlib.figure import Figure
    
        figsize = [x / float(dpi) for x in (arr.shape[1], arr.shape[0])]
        fig = Figure(figsize=figsize, dpi=dpi, frameon=False)
        canvas = FigureCanvas(fig)
        im = fig.figimage(arr, cmap=cmap, vmin=vmin, vmax=vmax, origin=origin)
        fig.savefig(fname, dpi=dpi, format=format, transparent=True)
    
    
    def create_bitmap(self):
        if self.difference and self.imgfit_tot!= None:
            outimg=self.img-self.imgfit_tot
        elif self.compensate_offset:
            outimg=self.img - self.img_background
        else:
            outimg=self.img
        
        
        if self.reference_diff and (self.reference_img!=None):
            outimg-=self.reference_img
        if self.filter:
            outimg=scipy.ndimage.filters.gaussian_filter(outimg,[float(self.filter_sigma),float(self.filter_sigma)])
            #outimg=scipy.signal.wiener(outimg,[float(self.filter_sigma),float(self.filter_sigma)])
            #outimg=scipy.ndimage.filters.median_filter(outimg,[float(self.filter_sigma),float(self.filter_sigma)])
        outimg*=self.filter_scalar
        
        self.calc_residuals()
        
        return outimg
    
    
    def calc_residuals(self):
        if isinstance(self.fitpars, (fitting.FitParsBimodal2d, fitting.FitParsTF2d)):
            
            if self.reference_diff and self.reference_img!=None:
                myimgroi = self.img - self.reference_img
            else:
                myimgroi=self.img
                
            if self.filter:
                myimgroi=scipy.ndimage.filters.gaussian_filter(myimgroi,[float(self.filter_sigma),float(self.filter_sigma)])
            
            r = self.roi
            myimgroi=myimgroi[r.y,r.x]
            
            if isinstance(self.fitpars, fitting.FitParsBimodal2d):
                #imgfitgauss = (self.imgfit[1] - self.img_background)
                imgfitbec   = (self.imgfit[0] - self.imgfit[1])
                mymask = numpy.ma.getmask(numpy.ma.masked_less_equal(imgfitbec, numpy.ma.maximum(imgfitbec)*0.00001))
                myimgfitbec = numpy.ma.masked_array(data=self.imgfit[0], mask=mymask)
            elif isinstance(self.fitpars, fitting.FitParsTF2d):
                imgfitbec = self.imgfit[0]  - self.img_background
                myimgfitbec = self.imgfit[0]
                
            max = numpy.ma.maximum(myimgfitbec)
            min = numpy.ma.minimum(myimgfitbec)
            
            masked_fit = numpy.ma.masked_less_equal(myimgfitbec, (max-min)*0.2+min)
            numpy.ma.set_fill_value(masked_fit, numpy.nan)
            numpy.ma.set_fill_value(myimgroi, numpy.nan)
            self.bec_residuals = numpy.ma.mean(((myimgroi - masked_fit)/imgfitbec)**2)
        else:
            self.bec_residuals=None
        
        self.update_parameters()
        
        return self.bec_residuals
    
    
    
    def save_bitmap(self,filename,cmap):
        min,max=self.himg.get_clim()
        self.mpl_imsave(filename, self.create_bitmap(), cmap=cmap,vmin= min, vmax=max)

    def clear_image_axis_obsolete(self):
        """clean image axis from old contours. for this, clear axis
        completely and readd artist (cursors, markers, ...)"""
       
        self.aximg.clear()
        self.aximg.set_axis_off()
        self.aximg.add_artist(self.himg)
        for marker in self.hrx + self.hry:
            self.aximg.add_artist(marker.cursor.lineh)
            self.aximg.add_artist(marker.cursor.linev)
            self.aximg.add_artist(marker.line)
        for marker in self.markers:
            self.aximg.add_artist(marker.marker)
        
        self.aximg.set_autoscale_on(False)

    def clear_contours(self):
        """remove contour lines from image axis"""
        self.aximg.collections = [h for h in self.aximg.collections if h not in self.hcontours]
        self.hcontours = []

    def update_contours(self):
        try:
            self.clear_contours()
            if self.show_contours:
                xlim = numpy.copy(self.aximg.get_xlim())
                ylim = numpy.copy(self.aximg.get_ylim())

                if isinstance(self.fitpars, fitting.FitParsBimodal2d):
                    imgfitgauss = (self.imgfit[1] - self.img_background)
                    imgfitbec   = (self.imgfit[0] - self.imgfit[1])
                    cs = self.aximg.contour(self.roi.xrange_clipped(self.img),
                               self.roi.yrange_clipped(self.img),
                               imgfitgauss,
                               [numpy.exp(-2)*imgfitgauss.max(),],
                               colors = 'w',
                               alpha = 0.5,
                               linewidths = 1.0,
                               linestyle = 'dashed')

                    cs2 = self.aximg.contour(self.roi.xrange_clipped(self.img),
                               self.roi.yrange_clipped(self.img),
                               imgfitbec,
                               [0.05*imgfitbec.max(),],
                               colors = 'w',
                               alpha = 0.5,
                               linewidths = 1.0)#
                    cs2.collections[0].set_linestyle('dashed')
                    self.hcontours = cs.collections + cs2.collections

                elif isinstance(self.fitpars, fitting.FitParsGauss2d):
                    cs = self.aximg.contour(self.roi.xrange_clipped(self.img),
                               self.roi.yrange_clipped(self.img),
                               self.imgfit[0],
                               [0.37*self.fitpars.OD],
                               colors = 'w',
                               alpha = 0.5,
                               linewidths = 1.0)#
                    self.hcontours = cs.collections

                self.aximg.set_xlim(xlim)
                self.aximg.set_ylim(ylim)

            else:
                self.hcontours = []
                
        except Exception, e:
            print "Error in update contours: ", e
            import traceback
            traceback.print_exc()

    def update(self):
        "calculate profiles and fit results, perform redraw"
        self.invalidate_profiles_fit()
        self.invalidate_parameters()
        self.bec_residuals = None
        self.clear_contours()
        self.create_image()
        self.himg.set_data(self.img)
        self.redraw()
        
        self.fit.rotation_active=self.rotation_active
        self.fit.rotation_angle=self.rotation_angle

        self.fitjobID += 1
        delayedresult.startWorker(consumer = self,
                                  workerFn = self._fit_result_producer,
                                  wargs = (self.fit, #worker arguments
                                           self.img.copy(),
                                           copy.copy(self.roi)), 
                                  cargs = (FitCompletedEvent,), #consumer arguments
                                  ckwargs = {'resultAttr': 'fitresults',
                                             'target': self.target_results}, 
                                  jobID = self.fitjobID)
        
    def _fit_result_producer(self, fit, img, roi):
        if self.difference:
            return fit.do_fit(img, roi, True)
        else:
            return fit.do_fit(img, roi)

    def OnFitCompleted(self, event):
        fitresult = event.fitresults
        jobID = fitresult.getJobID()
        if jobID != self.fitjobID:
            print "Received fit results from outdated fit. Ignoring these results."
            return
        try:
            self.imgfit, self.img_background, self.fitpars, self.imgfit_tot = fitresult.get()
        except Exception, e:
            #if fit failed, create save settings, don't update profiles and
            #contours
            print "Error in multi-threaded fit:", e
            self.imgfit = numpy.empty(shape=(0,0))
            self.img_background = numpy.array([0])
            self.fitpars = fitting.FitParsNoFit()
        else:
            self.update_profiles_fit()
            self.update_contours()

        self.update_profiles()
        self.update_parameters()
        self.update_image()
        self.redraw()

        #tell the world that fit is completed
        if jobID > 1:
            evt = FitResultsEvent(id=0,
                                  source = self.id,
                                  target = event.target,
                                  fitpars = self.fitpars)
            wx.PostEvent(self.canvas, evt)

    def analyze_image(self, img, target):
        """show image and perform fit. store target information for
        later use."""
        self.rawimg = img
        self.target_results = target
        self.update()
        
    def create_image(self):
        """
        apply image filtering, create self.img, (don't update image display data)
        """
        img = self.rawimg
        
        ##compensate for finite optical density
        ODmax = self.fit.imaging_pars.ODmax #TODO: this might fail if
                                            #fit class is changed.
        if ODmax > 0:
            img = numpy.log((1 - numpy.exp(- ODmax)) / (numpy.exp(- img) - numpy.exp(- ODmax)))
            #TODO: remove invalid entries
            #TODO: this should not happen here!
        
        self.img = ma.array(img, mask= ~ numpy.isfinite(img))
        if ODmax > 0:
            self.img.set_fill_value(ODmax)
        else:
            self.img.set_fill_value(3) #which to take?
        
    def redraw(self):   
        self.fig.canvas.draw()
        
    def invalidate_parameters(self):
        self.htxt.set_alpha(0.6)
        self.htxt.set_color('k')
        self.htxt.set_text('fitting...')

    def update_parameters(self):
        s = unicode(self.fitpars)
        if self.bec_residuals!=None:
            s+="\nBEC residuals:\n%.5f"%self.bec_residuals
            
        self.htxt.set_text(s)

        if self.fitpars.valid:
            self.htxt.set_color('k')
            self.htxt.set_alpha(1)
        else:
            self.htxt.set_color('r')
            self.htxt.set_alpha(0.2)
        
    def invalidate_profiles_fit(self):
        for h in self.hhproffit + self.hvproffit:
            h.set_data([], [])

    def update_profiles_fit(self):
        for k, img in enumerate(self.imgfit):
            if self.compensate_offset:
                cimg = img - self.img_background
            else:
                cimg = img

            self.hhproffit[k].set_data(self.roi.xrange_clipped(self.img), cimg.sum(0))
            self.hvproffit[k].set_data(cimg.sum(1), self.roi.yrange_clipped(self.img))

        if self.compensate_offset:
            self.hhproffit[ - 1].set_data([0, 1250], [0, 0])
            self.hvproffit[ - 1].set_data([0, 0], [0, 1040])

        
        
        
    #def do_fit(self):
    #    self.imgfit, self.img_background, self.fitpars = self.fit.do_fit(self.img, self.roi)

    def update_image(self):
        self.himg.set_data(self.create_bitmap())

                
    def set_compensate_offset(self, value=True):
        self.compensate_offset = value
        self.update() #TODO: some optimazation?
        #self.update_profiles()
        #self.invalidate_profiles_fit()
        #self.update_profiles_fit()
        #self.update_image()
        #self.redraw()

    def set_marker_guess(self, value=False):
        self.marker_guess = value
        self.update() 

    def set_roi(self, roi):
        "set roi, update selectors to match roi"
        self.roi = roi

        for selector in self.hrx + self.hry:
            selector.update_from_roi(roi)

        #TODO: update fits?

    def center_roi(self, border=50):
        self.axhprof.set_xlim(self.roi.xmin - border, self.roi.xmax + border)
        self.axvprof.set_ylim(self.roi.ymax + border, self.roi.ymin - border)

    def onpick(self, event):
        need_update = False
        
        "handle pick events. Knows about handling events generated by ROI-picks"
        try:
            self.roi.setx(event.roix, event.num)
        except AttributeError:
            pass
        else:
            need_update = True

        try:
            self.roi.sety(event.roiy, event.num)
        except AttributeError:
            pass
        else:
            need_update = True

        try:
            x, y = event.x, event.y
        except AttributeError:
            pass
        else:
            self.redraw()
        
        if need_update:
            self.update()

    def OnAutoscale(self, event):
        prof = self.vprof
        mi = prof.min()
        ma = prof.max()
        if self.compensate_offset:
            mi = 0
        d = 0.2*(ma - mi)
        self.axvprof.set_xlim(mi-d, ma+d)
        
        prof = self.hprof
        mi = prof.min()
        ma = prof.max()
        if self.compensate_offset:
            mi = 0
        d = 0.2*(ma - mi)
        self.axhprof.set_ylim(mi-d, ma+d)
        
        border = 50
        self.axhprof.set_xlim(self.roi.xmin - border, self.roi.xmax + border)
        self.axvprof.set_ylim(self.roi.ymax + border, self.roi.ymin - border)
        
        self.redraw()
    
##results name

class ImgAppAui(wx.App):

    #Menu: View
    ID_ShowH = wx.NewId()
    ID_ShowV = wx.NewId()
    ID_ShowTools = wx.NewId()
    ID_ShowSavedImages = wx.NewId()
    
    #Menu: Perspectives
    ID_Perspective = wx.NewId()
    ID_PerspectiveFull = wx.NewId()
    ID_PerspectiveH = wx.NewId()
    ID_PerspectiveV = wx.NewId()
    
    #Menu: Fit
    ID_Reload = wx.NewId()
    
    ID_FitShowContours = wx.NewId()
    
    ID_FitHNone = wx.NewId()
    ID_FitHGauss1D = wx.NewId()
    ID_FitHLorentzGauss1D = wx.NewId()
    ID_FitHGaussGauss1D = wx.NewId()
    ID_FitHGauss = wx.NewId()
    ID_FitHGaussSym = wx.NewId()
    ID_FitHGaussBose = wx.NewId()
    ID_FitHTF = wx.NewId()
    ID_FitHTFSlice = wx.NewId()
    ID_FitHBimodal = wx.NewId()
    ID_FitHBimodalSplit = wx.NewId()
    ID_FitHBoseBimodal = wx.NewId()
    ID_FitHBimodalGaussGauss = wx.NewId()
    ID_FitHNGauss = wx.NewId()
    
    ID_FitVNone = wx.NewId()
    ID_FitVGauss1D = wx.NewId()
    ID_FitVLorentzGauss1D = wx.NewId()
    ID_FitVGaussGauss1D = wx.NewId()
    ID_FitVGauss = wx.NewId()
    ID_FitVGaussSym = wx.NewId()
    ID_FitVGaussBose = wx.NewId()
    ID_FitVTF = wx.NewId()
    ID_FitVTFSlice = wx.NewId()
    ID_FitVBimodal = wx.NewId()
    ID_FitVBimodalSplit = wx.NewId()
    ID_FitVBimodalGaussGauss = wx.NewId()
    ID_FitVBoseBimodal = wx.NewId()

    
    ID_CompensateOffsetH = wx.NewId()
    ID_CompensateOffsetV = wx.NewId()
    
    ID_MarkerGuessH = wx.NewId()
    ID_MarkerGuessV = wx.NewId()

    #Menu: Display settings
    ID_V_ContrastHigh = wx.NewId()
    ID_V_ContrastUltra = wx.NewId()
    ID_V_ContrastNormal = wx.NewId()
    ID_V_ContrastLow = wx.NewId()
    
    ID_H_ContrastHigh = wx.NewId()
    ID_H_ContrastUltra = wx.NewId()
    ID_H_ContrastNormal = wx.NewId()
    ID_H_ContrastLow = wx.NewId()
    
    ID_MarkerA = wx.NewId()
    ID_MarkerB = wx.NewId()
    ID_MarkerC = wx.NewId()
    
    #Menu: Results
    ID_ResultsPlot = wx.NewId()
    ID_ResultsSave = wx.NewId()
    ID_ResultsSaveAs = wx.NewId()
    ID_ResultsLoad = wx.NewId()
    ID_ResultsSaveTemplate = wx.NewId()
    ID_ResultsApplyTemplate = wx.NewId()
    ID_ResultsNew = wx.NewId()
    
    #Menu: Settings
    ID_SettingsSave = wx.NewId()
    ID_SettingsLoad = wx.NewId()
    
    #Menu: Autosav
    ID_AutoSaveAbs = wx.NewId()
    ID_AutoSaveRaw = wx.NewId()
    
    #Toolbar:
    ID_SaveImage = wx.NewId()
    ID_Autosave = wx.NewId()
    
    ID_ReloadButton = wx.NewId()
    ID_Autoreload = wx.NewId()
    
    ID_RecordData = wx.NewId()
    
    ID_ImagingChoice = wx.NewId()
    ID_Imaging_AB = wx.NewId()
    ID_ExpansionTimeV = wx.NewId()
    ID_ExpansionTimeH = wx.NewId()
    ID_OpticalDensityMaxV = wx.NewId()
    ID_OpticalDensityMaxH = wx.NewId()
    
    ID_DifferenceV= wx.NewId()
    ID_DifferenceH = wx.NewId()
    ID_FilterV = wx.NewId()
    ID_FilterH = wx.NewId()
    ID_FilterScalarV = wx.NewId()
    ID_FilterScalarH = wx.NewId()
    ID_FilterSigmaV = wx.NewId()
    ID_FilterSigmaH = wx.NewId()
    ID_SaveBitmap = wx.NewId()
    ID_ReferenceV = wx.NewId()
    ID_ReferenceH = wx.NewId()
    ID_CreateReference = wx.NewId()
    
    ID_CMapChoiceV = wx.NewId()
    ID_CMapChoiceH = wx.NewId()
    
    ID_RotateV = wx.NewId()
    ID_RotateH = wx.NewId()
    
    ID_VminV = wx.NewId()
    ID_VminH = wx.NewId()
    
    #ImageTree
    ID_AnalyzeImageButton = wx.NewId()
    ID_ImageTreeReloadSavedImage = wx.NewId()
    ID_ImageTreeRescan = wx.NewId()
    #
    
    #other settings
    imagefilename = os.path.abspath(settings.imagefile)
    try: # used to keep compatibility with older settings
        watchedfiles = [os.path.abspath(f) for f in settings.watchedfiles]
    except NameError:
        watchedfiles = [self.imagefilename]
    
    # TODO: Deprecated?
    rawimg1filename = settings.rawimage1file
    rawimg2filename = settings.rawimage2file
    rawimg3filename = settings.rawimage3file
    autosave_abs = False
    autosave_raw = False

    
    roiV = [ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900)] #list of ROIs for all imaging settings
    roiH = [ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900),ROI(400, 1200, 300, 900)]



    imaging_parlist = [{'V': imagingpars.ImagingParsVertical(), 'H': imagingpars.ImagingParsHorizontal()},
                       {'V': imagingpars.ImagingParsPixelSize(), 'H': imagingpars.ImagingParsPixelSize()},
                       {'V': imagingpars.ImagingParsDark(), 'H': imagingpars.ImagingParsDark()},
                        ]


#-------------------------
  
    imaging_pars = imaging_parlist[0]
    markersHpositions = []
    markersVpositions = []
    


    def select_imaging(self, n):
        """
        select new imaging parameters. sync controls with new settings. 
        keep imaging pars of fitting objects in sync
        """
        
        self.imaging_pars = self.imaging_parlist[n]
        self.sync_imaging_pars_expansion_time(updateControls=True)
        self.sync_imaging_pars_maximum_optical_density(updateControls=True)

        self.H.fit.set_imaging_pars(self.imaging_pars['H'])
        self.H.fitpars.imaging_pars = self.imaging_pars['H']
        
        self.V.fit.set_imaging_pars(self.imaging_pars['V'])
        self.V.fitpars.imaging_pars = self.imaging_pars['V']

    def sync_imaging_pars_expansion_time(self, updateControls=False):
        "keep imaging_pars in sync with controls for expansion time"
        if updateControls:
            self.expansion_time_V.Value = self.imaging_pars['V'].expansion_time
            self.expansion_time_H.Value = self.imaging_pars['H'].expansion_time
        else:
            self.imaging_pars['V'].expansion_time = self.expansion_time_V.Value
            self.imaging_pars['H'].expansion_time = self.expansion_time_H.Value
        
    def sync_imaging_pars_maximum_optical_density(self, updateControls=False):
        if updateControls == False:
            #get values from control
            try:
                ODmaxV = float(self.entry_optical_density_max_V.Value)
            except ValueError:
                ODmaxV = 0
                
            try:
                ODmaxH = float(self.entry_optical_density_max_H.Value)
            except ValueError:
                ODmaxH = 0
            
            self.imaging_pars['V'].ODmax = ODmaxV
            self.imaging_pars['H'].ODmax = ODmaxH
        else:
            #update controls
            ODmaxV = self.imaging_pars['V'].ODmax
            ODmaxH = self.imaging_pars['H'].ODmax
        
        #always update control
        self.entry_optical_density_max_V.Value = "%.1f" % ODmaxV
        self.entry_optical_density_max_H.Value = "%.1f" % ODmaxH    

    def select_roi(self, sel):
        self.H.set_roi(self.roiH[sel])
        self.V.set_roi(self.roiV[sel])
    
    def OnChoiceCMap(self,event):
        if event.Id==self.ID_CMapChoiceV:
            sel=self.cmap_choiceV.GetStringSelection()
            self.imaging_pars['V'].palette=sel
            self.V.himg.set_cmap(pylab.get_cmap(self.imaging_pars['V'].palette))
            self.V.update_image()
            self.V.redraw()
            
        elif event.Id==self.ID_CMapChoiceH:
            sel=self.cmap_choiceH.GetStringSelection()
            self.imaging_pars['H'].palette=sel
            self.H.himg.set_cmap(pylab.get_cmap(self.imaging_pars['H'].palette))
            self.H.update_image()
            self.H.redraw()
        
    def OnChoiceImaging(self, event):
        sel = self.imaging_choice.GetSelection()

        self.select_imaging(sel)
        self.select_roi(sel)

        self.H.center_roi()
        self.V.center_roi()

        self.H.update()
        self.V.update()
        
        self.H.himg.set_cmap(pylab.get_cmap(self.imaging_pars['H'].palette))
        self.V.himg.set_cmap(pylab.get_cmap(self.imaging_pars['V'].palette))
        self.cmap_choiceV.SetStringSelection(self.imaging_pars['V'].palette)
        self.cmap_choiceH.SetStringSelection(self.imaging_pars['H'].palette)

        
    def OnChoice_AB(self, event):
        self.load_image(self.imagefilename)
        
    def OnChangeExpansionTime(self, event):
        id = event.GetId()
        self.sync_imaging_pars_expansion_time()
        if id == self.ID_ExpansionTimeV:
            self.V.update_parameters()
            self.V.redraw()
            self.results.UpdateResults({'V': self.V.fitpars})
        elif id == self.ID_ExpansionTimeH:
            self.H.update_parameters()
            self.H.redraw()
            self.results.UpdateResults({'H': self.H.fitpars})
        
    def OnChangeMaximumOpticalDensity(self, event):
        id = event.Id
        self.sync_imaging_pars_maximum_optical_density()
        #TODO: need reload image
        if id == self.ID_OpticalDensityMaxV:
            self.V.update()
        elif id == self.ID_OpticalDensityMaxH:
            self.H.update()


    def OnChangeRotation(self, event):
        id = event.Id
        #TODO: need reload image
        if id == self.ID_RotateV:
            self.V.rotation_angle = float(self.rotation_V.GetValue())
            if self.V.rotation_angle != 0.0:
                self.V.rotation_active = True
            else:
                self.V.rotation_active = False
            self.V.update()
        elif id == self.ID_RotateH:
            self.H.rotation_angle = float(self.rotation_H.GetValue())
            if self.H.rotation_angle != 0.0:
                self.H.rotation_active = True
            else:
                self.H.rotation_active = False
            self.H.update()
            
            
    def OnChangeVmin(self, event):
        id = event.Id
        #TODO: need reload image
        if id == self.ID_VminV:
            vmin = float(self.vmin_V.GetValue())
            self.V.himg.set_clim(vmin=vmin)
            self.V.redraw()
        elif id == self.ID_VminH:
            vmin = float(self.vmin_H.GetValue())
            self.H.himg.set_clim(vmin=vmin)
            self.H.redraw()


    def OnChangeFilterSigma(self, event):
        id = event.Id
        if id == self.ID_FilterSigmaV:
            self.V.filter_sigma=float(self.entry_filter_sigma_V.GetValue())
            self.V.update_image()
            self.V.redraw()
        elif id == self.ID_FilterSigmaH:
            self.H.filter_sigma=float(self.entry_filter_sigma_H.GetValue())
            self.H.update_image()
            self.H.redraw()
        self.UpdateResiduals()
    
    def UpdateResiduals(self):
        try:
            if self.V.target_results!=None:
                row=self.V.target_results["row"]
            else:
                row=self.results.active_row
            
            if(self.H.bec_residuals!= None):    
                self.results.SetValueNamed(row, "resid H", self.H.bec_residuals)
            
            if(self.V.bec_residuals!= None):    
                self.results.SetValueNamed(row, "resid V", self.V.bec_residuals)
            
            self.results.GetView().Refresh()
        except AttributeError:
            print "error while updating residuals"
            
    def OnChangeFilterScalar(self, event):
        id = event.Id
        if id == self.ID_FilterScalarV:
            self.V.filter_scalar=float(self.entry_filter_scalar_V.GetValue())
            self.V.update_image()
            self.V.redraw()
        elif id == self.ID_FilterScalarH:
            self.H.filter_scalar=float(self.entry_filter_scalar_H.GetValue())
            self.H.update_image()
            self.H.redraw()
        #self.UpdateResiduals()
            
    def OnMenuResultsSave(self, event):
        id = event.GetId()

        filename = self.results.filename
        if id == self.ID_ResultsSaveAs or filename is None:
            filename = self.results_save_filename_ask(self.results)
        
        self.results_save(filename, self.results)

    def results_save_filename_ask(self, results):
        """return filename or None"""
        filename_proposal = time.strftime("%Y%m%d") + '-' + self.results.name 
        savedialog = wx.FileDialog(
            self.frame,
            message="Save results as ...",
            defaultDir=self.get_data_dir(),
            defaultFile=filename_proposal,
            wildcard="CSV file (*.csv)|*.csv|All files (*.*)|*.*",
            style=wx.SAVE | wx.FD_CHANGE_DIR)

        if savedialog.ShowModal() == wx.ID_OK:
            filename = savedialog.GetPath()
        else:
            filename = None

        savedialog.Destroy()
        return filename

    def results_save(self, filename, results, forcesave = False):#False):

        """save results to file, check for overwriting. If forcesave
        is True, don't give up."""

        success = False
        while not success:

            if filename is None:
                if forcesave:
                    filename = self.results_save_filename_ask(results)
                    if filename is None:
                        #continue
                        print "data not saved"
                        break
                else:
                    print "data not saved!"
                    results.filename = None
                    return success

            #check if file exists
            if filename is not None and os.access(filename, os.F_OK):
                #file already exists, ask to overwrite it
                MB = wx.MessageDialog(self.frame,
                                      "File " + filename + 
                                      " already exists. \nDo you want to overwrite it?",
                                      caption="Save File ...",
                                      style=wx.YES_NO | wx.ICON_EXCLAMATION,
                                      )
                answer = MB.ShowModal()
                MB.Destroy()
                if answer == wx.ID_YES:
                    pass
                elif forcesave:
                    #don't overwrite, ask again
                    filename = None
                    continue
                else:
                    print "data file not saved"
                    results.filename = None
                    return success
                
            try:
                results.save_data_csv(filename)
            except StandardError, e:
                import traceback
                traceback.print_exc()
                msg = wx.MessageDialog(self.frame,
                                       'Error saving data to file: "%s"\n%s'\
                                       %(filename, traceback.format_exc(e)),
                                       style = wx.OK | wx.ICON_ERROR,
                                       )
                msg.ShowModal()
                success = False
                filename = None
                results.filename = None
            else:
                results.filename = filename
                success = True

            return success
                
    def OnMenuResultsLoad(self, event):
        loaddialog = wx.FileDialog(
            self.frame,
            message="Load results",
            defaultDir=settings.imagesavepath,
            defaultFile='',
            wildcard="CSV file (*.csv)|*.csv",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)
        if loaddialog.ShowModal() == wx.ID_OK:
            fullpath = loaddialog.GetPath()
            dir, filename = os.path.split(fullpath)
            name, ext = os.path.splitext(filename)
            self.select_measurement(name)
            self.results.load_data_csv(fullpath)

        loaddialog.Destroy()

    def OnMenuResultsSaveTemplate(self, event):
        savedialog = wx.FileDialog(
            self.frame,
            message="Save settings to template file ...",
            defaultDir=settings.templatedir,
            defaultFile='template' + self.results.name,
            wildcard="template file (*.tpl)|*.tpl|All files (*.*)|*.*",
            style=wx.SAVE)

        if savedialog.ShowModal() == wx.ID_OK:
            templatefile = savedialog.GetPath()

            #save metadata to template file
            metadata = self.results.give_metadata()
            output = open(templatefile, 'wb')
            pickle.dump(metadata, output)
            output.close()

        savedialog.Destroy()


    def OnMenuResultsApplyTemplate(self, event):
        loaddialog = wx.FileDialog(
            self.frame,
            message="Load settings from template file ...",
            defaultDir=settings.templatedir,
            defaultFile='',
            wildcard="template file (*.tpl)|*.tpl|All files (*.*)|*.*",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)

        if loaddialog.ShowModal() == wx.ID_OK:
            infile = open(loaddialog.GetPath(), 'rb')
            metadata = pickle.load(infile)
            infile.close()
            self.results.GetView().ApplyTemplate(metadata)

        loaddialog.Destroy()

    def OnMenuResultsNew(self, event):
        """
        Create new measurement, based on actual measurement (copy settings).
        """

        #get metadata from active measurement
        metadata = self.results.give_metadata()

        #TODO: check if saved, save

        #make proposal for measurement
        allnames = self.measurement_combobox.Strings
        name = self.results.name #name of active measurement

        #try to match name like "meas 5"
        pattern = r"(?P<name>[a-zA-Z_]+\s*)(?P<nr>\d+)"
        match = re.match(pattern, name)
        if match is not None:
            #if name matches pattern "meas 5", use "meas 6" as
            #proposal. If this already used, try higher numbers
            gd = match.groupdict()
            prefix = gd['name']
            nr = int(gd['nr'])
            for k in range(100):
                proposal = "%s%d" % (prefix, nr + k)
                if proposal not in allnames:
                    break
        else:
            for k in range(2, 100):
                proposal = "%s %d" % (name.strip(), k)
                if proposal not in allnames:
                    break
        
        #ask for new measurement name, with proposal
        dlg = wx.TextEntryDialog(
            self.frame,
            message="New measurement name:",
            caption="Create new measurement ...",
            defaultValue=proposal)
        if dlg.ShowModal() == wx.ID_OK:
            newname = dlg.Value
        else:
            newname = None
        dlg.Destroy()

        #create new measurement (update combo box)
        #plot results

        if newname:
            self.select_measurement(newname)
            self.results.GetView().ApplyTemplate(metadata)
            self.plot_results()
    
    def OnCreateReferenceFromMask(self,event):
        table=self.gridpanel.pages[self.gridpanel.notebook.GetSelection()].Table
        maskedrows=table.rowmask
        maskedrowsfiles= [table.GetValueNamed(i, 'Filename') for i in range(len(maskedrows)) if maskedrows[i]]
        imagesavedir = self.get_data_dir(subdir = 'images')
        maskedrowsfilespath=[]
        for row in range(len(maskedrowsfiles)):
            if(maskedrowsfiles[row]!=""):
                maskedrowsfilespath.append(os.path.normpath(os.path.join(imagesavedir, maskedrowsfiles[row])))
        
        if len(maskedrowsfilespath)>0:
            mean_sis.create_reference(maskedrowsfilespath)
        
        
    def OnCreateReference(self, event):
        loaddialog = wx.FileDialog(
            self.frame,
            message="Select references",
            defaultDir=settings.imagesavepath,
            defaultFile='',
            wildcard="SIS file (*.sis)|*.sis",
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR | wx.FD_MULTIPLE)
        if loaddialog.ShowModal() == wx.ID_OK:
            fullpath = loaddialog.GetPaths()
            mean_sis.create_reference(fullpath)

        loaddialog.Destroy()    
        
    def OnReferenceSelect(self,event):
        if event.Id==self.ID_ReferenceH:
            if event.IsChecked():
                self.H.reference_diff=True
                _, self.H.reference_img = loadimg(settings.referencefile)
            else:
                self.H.reference_diff=False
            
            self.H.update_image()
            self.H.redraw()
            
        elif event.Id==self.ID_ReferenceV:
            if event.IsChecked():
                self.V.reference_diff=True
                self.V.reference_img, _= loadimg(settings.referencefile)
            else:
                self.V.reference_diff=False
            
            self.V.update_image()
            self.V.redraw()

                   
    def OnDifferenceSelectH(self,event):
        if event.IsChecked():
            self.H.difference=True
        else:
            self.H.difference=False
        self.H.update()
            
    def OnDifferenceSelectV(self,event):
        if event.IsChecked():
            self.V.difference=True
        else:
            self.V.difference=False
        self.V.update()
        
    def OnFilterSelectV(self,event):
        if event.IsChecked():
            self.V.filter=True
        else:
            self.V.filter=False
        self.V.update_image()
        self.V.redraw()
        
        self.UpdateResiduals()
        
    def OnFilterSelectH(self,event):
        if event.IsChecked():
            self.H.filter=True
        else:
            self.H.filter=False
        self.H.update_image()
        self.H.redraw()
        
        self.UpdateResiduals()

    def OnMenuShow(self, event):
        id = event.GetId()

        if id in [self.ID_ShowH]:
            pane = self.mgr.GetPane("H")
            if not pane.IsShown():
                pane.Float().Show()
                self.mgr.Update()

        if id in [self.ID_ShowV]:
            pane = self.mgr.GetPane("V")
            if not pane.IsShown():
                pane.Float().Show()
                self.mgr.Update()

        if id in [self.ID_ShowTools]:
            pane = self.mgr.GetPane('toolbar')
            if not pane.IsShown():
                pane.Show().ToolbarPane()
                self.mgr.Update()

        if id in [self.ID_ShowSavedImages]:
            panetree = self.mgr.GetPane('imagetree')
            paneimgna = self.mgr.GetPane('savedimageH')
            paneimgk = self.mgr.GetPane('savedimageV')
            if not panetree.IsShown():
                panetree.Right().Show()

            if not paneimgna.IsShown():
                paneimgna.Float().Show()

            if not paneimgk.IsShown():
                paneimgk.Float().Show()

            self.mgr.Update()
            
    def OnMenuPerspective(self, event):
        #print self.mgr.SavePerspective()
        id = event.GetId()
        if id in [self.ID_PerspectiveFull]:
            self.mgr.LoadPerspective(self.perspective_full)
        if id in [self.ID_PerspectiveH]:
            self.mgr.LoadPerspective(self.perspective_H)
        if id in [self.ID_PerspectiveV]:
            self.mgr.LoadPerspective(self.perspective_V)

    def OnMenuFit(self, event):
        id = event.GetId()
        self.markersHpositions = [self.H.markers[0].position[0],
                          self.H.markers[0].position[1]]
        self.markersVpositions = [self.V.markers[0].position[0],
                          self.V.markers[0].position[1]]
        # side markers H
        for k, marker in enumerate(self.H.markers):
            if k > 0:
                self.markersHpositions.append(self.H.markers[k].position[0])
                self.markersHpositions.append(self.H.markers[k].position[1])
                # insert automatically symmetric side marker
                # TO DO: check if out of range
                self.markersHpositions.append(2*self.markersHpositions[0]-self.H.markers[k].position[0])
                self.markersHpositions.append(2*self.markersHpositions[1]-self.H.markers[k].position[1])
        # side markers V
        for k, marker in enumerate(self.V.markers):
            if k > 0: 
                self.markersVpositions.append(self.V.markers[k].position[0])
                self.markersVpositions.append(self.V.markers[k].position[1])
                # insert automatically symmetric side marker
                # TO DO: check if out of range
                self.markersVpositions.append(2*self.markersVpositions[0]-self.V.markers[k].position[0])
                self.markersVpositions.append(2*self.markersVpositions[1]-self.V.markers[k].position[1])
        
                
        print 'markers H', self.markersHpositions
        print 'markers V', self.markersVpositions
            

        if id in [self.ID_FitShowContours]:
            self.H.show_contours = event.IsChecked()
            self.V.show_contours = event.IsChecked()
            self.H.update()
            self.V.update()

        if id in [self.ID_FitHNone]:
            self.H.fit = fitting.NoFit(self.imaging_pars['H'])
            self.H.update()
            
        if id in [self.ID_FitHGauss1D]:
            self.H.fit = fitting.Gauss1d(self.imaging_pars['H'])
            self.H.update()

        if id in [self.ID_FitHLorentzGauss1D]:
            self.H.fit = fitting.LorentzGauss1d(self.imaging_pars['H'])
            self.H.update()

        if id in [self.ID_FitHGaussGauss1D]:
            self.H.fit = fitting.GaussGauss1d(
                        self.imaging_pars['H'], self.markersHpositions)
            self.H.update()

        if id in [self.ID_FitHGauss]:
            self.H.fit = fitting.Gauss2d(self.imaging_pars['H'])
            self.H.update()

        if id in [self.ID_FitHGaussSym]:
            self.H.fit = fitting.GaussSym2d(self.imaging_pars['H'])
            self.H.update()

        if id in [self.ID_FitHGaussBose]:
            self.H.fit = fitting.GaussBose2d(self.imaging_pars['H'])
            self.H.update()

        if id in [self.ID_FitHTF]:
            self.H.fit = fitting.ThomasFermi2d(self.imaging_pars['H'])
            self.H.update()     
            
        if id in [self.ID_FitHTFSlice]:
            self.H.fit = fitting.ThomasFermiSlice(self.imaging_pars['H'])
            self.H.update()     

        if id in [self.ID_FitHBimodal]:
            self.H.fit = fitting.Bimodal2d(self.imaging_pars['H'])
            self.H.update()
        
        if id in [self.ID_FitHBimodalSplit]:
            self.H.fit = fitting.Bimodal2dSplit(self.imaging_pars['H'])
            self.H.update()

        if id in [self.ID_FitHBimodalGaussGauss]:
            state = self.menu.FindItemById(self.ID_MarkerGuessH).IsChecked()
            self.H.set_marker_guess(state)
            if self.H.marker_guess:
                self.H.fit = fitting.BimodalGaussGauss2d(
                            self.imaging_pars['H'], self.markersHpositions
                            )
            else:
                self.H.fit = fitting.BimodalGaussGauss2d(
                    self.imaging_pars['H'])
            self.H.update()

        if id in [self.ID_FitHBoseBimodal]:
            self.H.fit = fitting.BoseBimodal2d(self.imaging_pars['H'])
            self.H.update()

        if id in [self.ID_FitHNGauss]:
            self.H.fit = fitting.BimodalGaussGauss2d(
                            self.imaging_pars['H'], self.markersHpositions
                            )
            self.H.fit = fitting.NGauss2d(self.imaging_pars['H'], 
                            self.markersHpositions
                            )
            self.H.update()

        if id in [self.ID_FitVNone]:
            self.V.fit = fitting.NoFit(self.imaging_pars['V'])
            self.V.update()
            
        if id in [self.ID_FitVGauss1D]:
            self.V.fit = fitting.Gauss1d(self.imaging_pars['V'])
            self.V.update()

        if id in [self.ID_FitVLorentzGauss1D]:
            self.V.fit = fitting.LorentzGauss1d(self.imaging_pars['V'], )
            self.V.update()

        if id in [self.ID_FitVGaussGauss1D]:
            self.V.fit = fitting.GaussGauss1d(
                       self.imaging_pars['V'], self.markersVpositions)
            self.V.update()

        if id in [self.ID_FitVGauss]:
            self.V.fit = fitting.Gauss2d(self.imaging_pars['V'])
            self.V.update()

        if id in [self.ID_FitVGaussSym]:
            self.V.fit = fitting.GaussSym2d(self.imaging_pars['V'])
            self.V.update()
            
        if id in [self.ID_FitVGaussBose]:
            self.V.fit = fitting.GaussBose2d(self.imaging_pars['V'])
            self.V.update()

        if id in [self.ID_FitVTF]:
            self.V.fit = fitting.ThomasFermi2d(self.imaging_pars['V'])
            self.V.update()       
        
        if id in [self.ID_FitVTFSlice]:
            self.V.fit = fitting.ThomasFermiSlice(self.imaging_pars['V'])
            self.V.update()   

        if id in [self.ID_FitVBimodal]:
            print 'bimodal V'
            self.V.fit = fitting.Bimodal2d(self.imaging_pars['V'])
            self.V.update()
            
        if id in [self.ID_FitVBimodalSplit]:
            print 'bimodal V - split mxy'
            self.V.fit = fitting.Bimodal2dSplit(self.imaging_pars['V'])
            self.V.update()

        if id in [self.ID_FitVBimodalGaussGauss]:
            state = self.menu.FindItemById(self.ID_MarkerGuessV).IsChecked()
            self.V.set_marker_guess(state)
            if self.V.marker_guess:
                self.V.fit = fitting.BimodalGaussGauss2d(
                            self.imaging_pars['V'], self.markersVpositions
                            )
            else:
                self.V.fit = fitting.BimodalGaussGauss2d(
                    self.imaging_pars['V'])
            self.V.update()
        
        if id in [self.ID_FitVBoseBimodal]:
            self.V.fit = fitting.BoseBimodal2d(self.imaging_pars['V'])
            self.V.update()
        


    def OnMenuCompensate(self, event):
        id = event.GetId()

        if id in [self.ID_CompensateOffsetH]:
            state = self.menu.FindItemById(id).IsChecked()
            self.H.set_compensate_offset(state)

        if id in [self.ID_CompensateOffsetV]:
            if self.menu.FindItemById(id).IsChecked():
                self.V.set_compensate_offset(True)
            else:
                self.V.set_compensate_offset(False)

    def OnMenuMarkerGuess(self, event):
        id = event.GetId()

        if id in [self.ID_MarkerGuessH]:
            state = self.menu.FindItemById(id).IsChecked()
            self.H.set_marker_guess(state)
        if id in [self.ID_MarkerGuessV]:
            if self.menu.FindItemById(id).IsChecked():
                self.V.set_marker_guess(True)
            else:
                self.V.set_marker_guess(False)

    def OnMenuContrast(self, event):
        id = event.GetId()

        V_contrast_dict = {
            self.ID_V_ContrastUltra:  0.2,
            self.ID_V_ContrastHigh:   0.5,
            self.ID_V_ContrastNormal: 1.5,
            self.ID_V_ContrastLow:    2.5,
            }

        H_contrast_dict = {
            self.ID_H_ContrastUltra:  0.2,
            self.ID_H_ContrastHigh:   0.5,
            self.ID_H_ContrastNormal: 1.5,
            self.ID_H_ContrastLow:    2.5,
            }

        cmax = V_contrast_dict.get(id)
        if cmax:
            self.V.himg.set_clim(vmax=cmax)
            self.V.redraw()
            #print id, cmax
            
        cmax = H_contrast_dict.get(id)
        if cmax:
            self.H.himg.set_clim(vmax=cmax)
            self.H.redraw()

    def OnMenuMarker(self, event):
        markerindex = [self.ID_MarkerA, self.ID_MarkerB, self.ID_MarkerC].index(event.Id)
        self.H.markers[markerindex].set_visible(event.IsChecked())
        self.V.markers[markerindex].set_visible(event.IsChecked())
        self.H.redraw()
        self.V.redraw()

    def OnMenuReload(self, event):
        id = event.GetId()

        if id in [self.ID_Reload]:
            self.PostReloadImageEvent()

    def OnMenuAutoSave(self, event):
        id = event.GetId()

        if id in [self.ID_AutoSaveAbs]:
            self.autosave_abs = event.IsChecked()

        if id in [self.ID_AutoSaveRaw]:
            self.autosave_raw = event.IsChecked()

    def OnReloadButton(self, event):
        self.PostReloadImageEvent()

    def OnSaveImage(self, event):
        evt = SaveImageEvent()
        wx.PostEvent(self.frame, evt)

###########
    def OnAutosaveToggle(self, event):
        """Callback for Auto Save Image Checkbox. Stores status to
        self.autosave"""
        if event.IsChecked():
            self.autosave = True
        else:
            self.autosave = False
###############

    def get_data_dir(self, subdir = ''):
        """change current directory to data dir"""
        directory = os.path.join(settings.imagesavepath,
                                 time.strftime("%Y/%Y-%m-%d/"),
                                 subdir)
        if not os.access(directory, os.F_OK):
            try:#try to create dir
                os.makedirs(directory)
            except OSError:
                print "cannot create data dir"
                return os.getcwd()
        return directory
    
    def OnSaveBitmap(self, event):
        imagesavedir = self.get_data_dir(subdir = 'bitmaps')
        if self.V.target_results!=None:
            name=self.V.target_results["name"]
            row=self.V.target_results["row"]
        else:
            row=self.results.active_row
            name=self.results.name
        
        imagesavefilename0 = "%s-%s-%04d_0.png" % (time.strftime("%Y%m%d"),
                                     name,
                                     int(row))
        imagesavefilename1 = "%s-%s-%04d_1.png" % (time.strftime("%Y%m%d"),
                                     name,
                                     int(row))
        
        
        imagesavefilename01 = "%s-%s-%04d.sis" % (time.strftime("%Y%m%d"),
                                     name,
                                     int(row))
        
        imagesavefilenamefull01 = os.path.normpath(os.path.join(imagesavedir, imagesavefilename01))
        
        
        imagesavefilenamefull0 = os.path.normpath(os.path.join(imagesavedir, imagesavefilename0))
        imagesavefilenamefull1 = os.path.normpath(os.path.join(imagesavedir, imagesavefilename1))
        
        #test if file already exists
        if os.access(imagesavefilenamefull0, os.F_OK):
            MB = wx.MessageDialog(self.frame,
                                  "Image file " + imagesavefilename0 + 
                                  " already exists. \nDo you want to overwrite it?",
                                  caption="Save Image File ...",
                                  style=wx.YES_NO | wx.ICON_EXCLAMATION,

                                  )
            answer = MB.ShowModal()
            MB.Destroy()
            if  answer == wx.ID_YES:
                pass
            else:
                print "image file not saved"
                return
        if os.access(imagesavefilenamefull1, os.F_OK):
            MB = wx.MessageDialog(self.frame,
                                  "Image file " + imagesavefilename1 + 
                                  " already exists. \nDo you want to overwrite it?",
                                  caption="Save Image File ...",
                                  style=wx.YES_NO | wx.ICON_EXCLAMATION,

                                  )
            answer = MB.ShowModal()
            MB.Destroy()
            if  answer == wx.ID_YES:
                pass
            else:
                print "image file not saved"
                return

        
        if os.access(imagesavefilenamefull01, os.F_OK):
            MB = wx.MessageDialog(self.frame,
                                  "Image file " + imagesavefilename01 + 
                                  " already exists. \nDo you want to overwrite it?",
                                  caption="Save Image File ...",
                                  style=wx.YES_NO | wx.ICON_EXCLAMATION,

                                  )
            answer = MB.ShowModal()
            MB.Destroy()
            if  answer == wx.ID_YES:
                pass
            else:
                print "image file not saved"
                return
        

        V_bit=self.V.create_bitmap()
        H_bit=self.H.create_bitmap()

        HV_bit=numpy.concatenate((numpy.ma.filled(V_bit,0), numpy.ma.filled(H_bit,0)))
        HV_bit=(HV_bit+1.0)*6553.6

        write_raw_image(imagesavefilenamefull01,HV_bit)
        
        self.V.save_bitmap(imagesavefilenamefull0,pylab.get_cmap(self.imaging_pars['V'].palette))
        self.H.save_bitmap(imagesavefilenamefull1,pylab.get_cmap(self.imaging_pars['H'].palette))
            
    def OnSaveImageEvent(self, event):
        imagesavedir = self.get_data_dir(subdir = 'images')
        imagesavefilename = "%s-%s-%04d.sis" % (time.strftime("%Y%m%d"),
                                             self.results.name,
                                             self.results.active_row)
        rawimg1savefilename = "r1-" + imagesavefilename
        rawimg2savefilename = "r2-" + imagesavefilename
        rawimg3savefilename = "r3-" + imagesavefilename

        imagesavefilenamefull = os.path.normpath(os.path.join(imagesavedir, imagesavefilename))
        rawimg1savefilenamefull = os.path.normpath(os.path.join(imagesavedir, rawimg1savefilename))
        rawimg2savefilenamefull = os.path.normpath(os.path.join(imagesavedir, rawimg2savefilename))
        rawimg3savefilenamefull = os.path.normpath(os.path.join(imagesavedir, rawimg3savefilename))

        #test if file already exists
        if os.access(imagesavefilenamefull, os.F_OK):
            MB = wx.MessageDialog(self.frame,
                                  "Image file " + imagesavefilename + 
                                  " already exists. \nDo you want to overwrite it?",
                                  caption="Save Image File ...",
                                  style=wx.YES_NO | wx.ICON_EXCLAMATION,

                                  )
            answer = MB.ShowModal()
            MB.Destroy()
            if  answer == wx.ID_YES:
                pass
            else:
                print "image file not saved"
                return
                                
        shutil.copy2(self.imagefilename, imagesavefilenamefull)
        statusbar_msg = "image saved as: " + imagesavefilename
        if self.autosave_raw:   # save also raw images
            statusbar_msg = statusbar_msg + "  (+ raw images as: r1-..., r2-..., r3-...)"
            shutil.copy2(self.rawimg1filename, rawimg1savefilenamefull)
            shutil.copy2(self.rawimg2filename, rawimg2savefilenamefull) 
            shutil.copy2(self.rawimg3filename, rawimg3savefilenamefull)

	# add autosave csv
	# self.results_save(os.path.join(self.get_data_dir(), self.results.filename), self.results)

        self.savebutton.SetBackgroundColour(wx.NamedColour("GREEN"))
        self.savebutton.Refresh()

        #self.saved_image_name.SetLabel(imagesavefilename)
        self.statusbar.SetStatusText(statusbar_msg)
        self.UpdateResultsFilename(imagesavefilename)

    def OnInit(self):
        ## create main frame

        splash = CamSplashScreen()
        splash.Show()

        self.SetAppName("Cam")
        
        self.frame = wx.Frame(None,
                              title="Image Display and Fit",
                              size=(1300, 750))

        #set main icon
        icons = wx.IconBundle()
        for icon in ['cam16.png',
                     'cam24.png',
                     'cam32.png',
                     'cam48.png',
                     ]:
                   
            icons.AddIconFromFile(os.path.join(settings.bitmappath, icon),
                                  wx.BITMAP_TYPE_PNG)
        self.frame.SetIcons(icons)

        #create manager
        self.mgr = wx.aui.AuiManager(self.frame,
                                     wx.aui.AUI_MGR_RECTANGLE_HINT 
 | wx.aui.AUI_MGR_ALLOW_FLOATING
                                     )

        self.mgr.SetDockSizeConstraint(0.5, 0.75)
        
        ## create menu bar
        
        #submenus Show
        view_menu = wx.Menu()
        view_menu.Append(self.ID_ShowH, "Show H")
        view_menu.Append(self.ID_ShowV, "Show V")
        view_menu.Append(self.ID_ShowTools, "Show Toolbar")
        view_menu.Append(self.ID_ShowSavedImages, "Browse saved images")
        self.frame.Bind(wx.EVT_MENU, self.OnMenuShow, id=self.ID_ShowH)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuShow, id=self.ID_ShowV)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuShow, id=self.ID_ShowTools)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuShow, id=self.ID_ShowSavedImages)

        #subsubmenu Perspectives
        perspective_menu = wx.Menu()
        perspective_menu.Append(self.ID_PerspectiveFull, "Show all")
        perspective_menu.Append(self.ID_PerspectiveH, "H only")
        perspective_menu.Append(self.ID_PerspectiveV, "V only")
        self.frame.Bind(wx.EVT_MENU, self.OnMenuPerspective, id=self.ID_PerspectiveFull)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuPerspective, id=self.ID_PerspectiveH)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuPerspective, id=self.ID_PerspectiveV)

        view_menu.AppendMenu(self.ID_Perspective, "Perspectives", perspective_menu)

        #submenu Display
        display_menu = wx.Menu()
        #subsubmenu V display
        display_V_menu = wx.Menu()
        display_V_menu.AppendRadioItem(self.ID_V_ContrastHigh, 'high', 'high contrast')
        display_V_menu.AppendRadioItem(self.ID_V_ContrastNormal, 'normal', 'low contrast')
        display_V_menu.AppendRadioItem(self.ID_V_ContrastLow, 'low', 'high contrast')
        display_V_menu.AppendRadioItem(self.ID_V_ContrastUltra, 'ultra', 'ultra high contrast')
        display_menu.AppendMenu(wx.NewId(),
                                "V Color Contrast",
                                display_V_menu)
        #subsubmenu H display
        display_H_menu = wx.Menu()
        display_H_menu.AppendRadioItem(self.ID_H_ContrastHigh, 'high', 'high contrast')
        display_H_menu.AppendRadioItem(self.ID_H_ContrastNormal, 'normal', 'low contrast')
        display_H_menu.AppendRadioItem(self.ID_H_ContrastLow, 'low', 'high contrast')
        display_H_menu.AppendRadioItem(self.ID_H_ContrastUltra, 'ultra', 'ultra high contrast')
        display_menu.AppendMenu(wx.NewId(),
                                "H Color Contrast",
                                display_H_menu)

        ##set defaults
        display_V_menu.Check(self.ID_V_ContrastNormal, True) #TODO
        display_H_menu.Check(self.ID_H_ContrastNormal, True) #TODO

        self.frame.Bind(wx.EVT_MENU_RANGE, self.OnMenuContrast,
                        id=self.ID_V_ContrastHigh,
                        id2=self.ID_H_ContrastLow)

        display_menu.AppendSeparator()
        display_menu.Append(self.ID_MarkerA,
                            'Main marker visible',
                            '',
                            wx.ITEM_CHECK)
        display_menu.Append(self.ID_MarkerB,
                            '2nd marker visible',
                            '',
                            wx.ITEM_CHECK)
        display_menu.Append(self.ID_MarkerC,
                            '3rd marker visible',
                            '',
                            wx.ITEM_CHECK)
        self.frame.Bind(wx.EVT_MENU_RANGE,
                        self.OnMenuMarker,
                        id=self.ID_MarkerA,
                        id2=self.ID_MarkerB)
        self.frame.Bind(wx.EVT_MENU_RANGE,
                        self.OnMenuMarker,
                        id=self.ID_MarkerC)

        #submenu Fit
        fit_menu = wx.Menu()
        fit_menu.Append(self.ID_Reload,
                        "Reload image\tF5",
                        "Reload image from file",
                        )
        fit_menu.AppendSeparator()
        sc = fit_menu.AppendCheckItem(self.ID_FitShowContours,
                        'Show Contours',
                        'Show Contour lines')
        sc.Check()
        fit_menu.AppendSeparator()
        
        #subsubmenu H,1 peak
        fit_H_menu = wx.Menu()
        fit_H_menu.AppendRadioItem(
            self.ID_FitHNone,
            "H no fit",
            "perform no fit for Horizontal",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHGauss1D,
            "H Gauss fit 1D",
            "perform 2 Gauss 1d fits for Horizontal",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHLorentzGauss1D,
            "H Lorentz(x)+Gauss(y) fit 1D",
            "perform two 1d fits for Horizontal",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHGaussGauss1D,
            "H GaussGauss(x)+Gauss(y) fit 1D",
            "perform 1d fits for Horizontal, with 2 gaussian along x",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHGauss,
            "H Gauss fit",
            "perform 2d Gauss fit for Horizontal",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHGaussSym,
            "H symmetric Gauss fit",
            "perform symmetric (sx = sy) 2d Gauss fit for Horizontal",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHGaussBose,
            "H Bose enhanced Gauss fit",
            "perform 2d Bose enhanced Gauss fit for Horizontal",
            )

        fit_H_menu.AppendRadioItem(
            self.ID_FitHTF,
            "H Thomas-Fermi fit",
            "perform 2d Thomas-Fermi fit distribution for Horizontal",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHTFSlice,
            "H Thomas-Fermi fit on a Slice",
            "perform 2d Thomas-Fermi fit distribution for Horizontal",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHBimodal,
            "H bimodal fit",
            "perform 2d fit of Gaussian + Thomas-Fermi distribution for Horizontal",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHBimodalSplit,
            "H bimodal fit - split mxy",
            "perform 2d fit of Gaussian + Thomas-Fermi distribution for Horizontal, with separate centers",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHBoseBimodal,
            "H bose enhanced bimodal fit",
            "perform 2d fit of Bose enhanced Gaussian + Thomas-Fermi distribution for Horizontal",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHBimodalGaussGauss,
            "H bimodal Gauss + Gauss fit",
            "perform 2d fit of Gaussian + Gaussian distribution for Horizontal",
            )
        fit_H_menu.AppendRadioItem(
            self.ID_FitHNGauss,
            "H Gauss xN fit",
            "perform 2d fit with sum of 3/5 Gaussian distribtions for Horizontal",
            )
        fit_menu.AppendMenu(wx.NewId(),
                                "H",
                                fit_H_menu)
        
        fit_menu.Append(self.ID_MarkerGuessH,
                        'H guess from markers',
                        'Guess peak positions from markers',
                        wx.ITEM_CHECK)
        fit_menu.Append(self.ID_CompensateOffsetH,
                        'H compensate background',
                        'Compensate offset of image background',
                        wx.ITEM_CHECK)

        fit_menu.AppendSeparator()
        
        # V
        fit_V_menu = wx.Menu()
        fit_V_menu.AppendRadioItem(
            self.ID_FitVNone,
            "V no fit",
            "perform no fit for Horizontal",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVGauss1D,
            "V Gauss fit 1D",
            "perform 2 Gauss 1d fits for Vertical",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVLorentzGauss1D,
            "V Lorentz(x)+Gauss(y) fit 1D",
            "perform two 1d fits for Vertical",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVGaussGauss1D,
            "V GaussGauss(x)+Gauss(y) fit 1D",
            "perform 1d fits for Vertical, with 2 gaussian along x",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVGauss,
            "V Gauss fit",
            "perform 2d Gauss fit for Vertical",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVGaussSym,
            "V symmetric Gauss fit",
            "perform symmetric (sx = sy) 2d Gauss fit for Vertical",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVGaussBose,
            "V Bose enhanced Gauss fit",
            "perform 2d Bose enhanced Gauss fit for Vertical",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVTF,
            "V Thomas-Fermi fit",
            "perform 2d Thomas-Fermi fit distribution for Vertical",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVTFSlice,
            "V Thomas-Fermi fit on a Slice",
            "perform 2d Thomas-Fermi fit distribution for Vertical",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVBimodal,
            "V bimodal fit",
            "perform 2d fit of Gaussian + Thomas-Fermi distribution for Vertical",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVBimodalSplit,
            "V bimodal fit - split mxy",
            "perform 2d fit of Gaussian + Thomas-Fermi distribution for Vertical, with separate centers",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVBoseBimodal,
            "V bose enhanced bimodal fit",
            "perform 2d fit of Bose enhanced Gaussian + Thomas-Fermi distribution for Horizontal",
            )
        fit_V_menu.AppendRadioItem(
            self.ID_FitVBimodalGaussGauss,
            "V bimodal Gauss + Gauss fit",
            "perform 2d fit of Gaussian + Gaussian distribution for Vertical",
            )

        fit_menu.AppendMenu(wx.NewId(),
                                "V",
                                fit_V_menu)
        
        fit_menu.Append(self.ID_MarkerGuessV,
                        'V guess from markers',
                        'Guess peak positions from markers',
                        wx.ITEM_CHECK)
        fit_menu.Append(self.ID_CompensateOffsetV,
                        'V compensate background',
                        'Compensate offset of image background',
                        wx.ITEM_CHECK)
        
        fit_menu.Check(self.ID_FitHNone, True)
        fit_menu.Check(self.ID_FitVNone, True)
        fit_menu.AppendSeparator()
    
        self.frame.Bind(wx.EVT_MENU,
                        self.OnMenuReload,
                        id=self.ID_Reload)
        self.frame.Bind(wx.EVT_MENU,
                        self.OnMenuFit,
                        id=self.ID_FitShowContours)
        self.frame.Bind(wx.EVT_MENU_RANGE,
                        self.OnMenuFit,
                        id=self.ID_FitHNone, #first item in fit H
                        id2=self.ID_FitVBoseBimodal) #last item in fit V
        self.frame.Bind(wx.EVT_MENU_RANGE,
                        self.OnMenuCompensate,
                        id=self.ID_CompensateOffsetH,
                        id2=self.ID_CompensateOffsetV)
        self.frame.Bind(wx.EVT_MENU_RANGE,
                        self.OnMenuMarkerGuess,
                        id=self.ID_MarkerGuessH,
                        id2=self.ID_MarkerGuessV)
        
        #submenu: Result
        result_menu = wx.Menu()
        result_menu.Append(self.ID_ResultsPlot,
                           "Plot Results",
                           "Plot Results.")
        
        result_menu.Append(self.ID_ResultsSave,
                           "Save", "Save fit results to file")
        result_menu.Append(self.ID_ResultsSaveAs,
                           "Save as ...", "Save fit results to file")
        result_menu.Append(self.ID_ResultsLoad,
                           "Load...", "Load results from file")

        result_menu.Append(self.ID_ResultsSaveTemplate,
                           "Save Template...", "Save settings to template file")
        result_menu.Append(self.ID_ResultsApplyTemplate,
                           "Apply Template", "Apply settings stored in template file")
        result_menu.Append(self.ID_ResultsNew,
                            "New results from template", "")
        
        
        self.frame.Bind(wx.EVT_MENU, self.OnMenuPlotResults, id=self.ID_ResultsPlot)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsSave, id=self.ID_ResultsSave)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsSave, id=self.ID_ResultsSaveAs)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsLoad, id=self.ID_ResultsLoad)

        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsSaveTemplate, id=self.ID_ResultsSaveTemplate)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsApplyTemplate, id=self.ID_ResultsApplyTemplate)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuResultsNew, id=self.ID_ResultsNew)
        
        
        #menu: Autosave
        autosave_menu = wx.Menu()
        autosave_menu.Append(self.ID_AutoSaveAbs, 'Save abs image','',wx.ITEM_CHECK)
        autosave_menu.Append(self.ID_AutoSaveRaw, 'Save raw images','',wx.ITEM_CHECK)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuAutoSave, id=self.ID_AutoSaveAbs)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuAutoSave, id=self.ID_AutoSaveRaw)
        
        #menu: Settings
        settings_menu = wx.Menu()
        settings_menu.Append(self.ID_SettingsSave, "Save", "Save Settings")
        settings_menu.Append(self.ID_SettingsLoad, "Load", "Load Settings")

        self.frame.Bind(wx.EVT_MENU, self.OnMenuSettingsSave, id=self.ID_SettingsSave)
        self.frame.Bind(wx.EVT_MENU, self.OnMenuSettingsLoad, id=self.ID_SettingsLoad)
        
        
        #Menus: finalize
        self.menu = wx.MenuBar()

        self.menu.Append(view_menu, "Show")
        self.menu.Append(fit_menu, "Fit")
        self.menu.Append(display_menu, "Display")
        self.menu.Append(result_menu, "Results")
        self.menu.Append(autosave_menu, "AutoSave")
        self.menu.Append(settings_menu, "Settings")

        self.frame.SetMenuBar(self.menu)



        ## create center panel
        self.gridpanel = FitResultTableGrid.TabbedGridPanel(self.frame)
        
        self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CLOSE, self.OnResultsPageClose, self.gridpanel.notebook)
        self.Bind(wx.aui.EVT_AUINOTEBOOK_PAGE_CLOSED, self.OnResultsPageClosed, self.gridpanel.notebook)
        
        self.results_filename = None

        self.mgr.AddPane(self.gridpanel,
                         wx.aui.AuiPaneInfo().
                         Name("results").Caption("Results")
                         .CenterPane().BestSize(wx.Size(800, 100))
                         )

        ## create and initializer image panels

        #imgH, imgH = loadimg(self.imagefilename)   #####added 14-12-2012
        #imgV, imgV = loadimg(self.imagefilename.replace(".sis","_1.sis"))   #####added 14-12-2012
        imgV, imgH = loadimg(self.imagefilename)   #####original 14-12-2012
        imgV = ma.array(imgV, mask= ~ numpy.isfinite(imgV))
        imgH = ma.array(imgH, mask= ~ numpy.isfinite(imgH))
        
        self.V = ImgPanel(self.frame, 'V', imgV, self.roiV[2], fitting.NoFit(self.imaging_pars['V']), pylab.get_cmap(self.imaging_pars['V'].palette))
        self.H = ImgPanel(self.frame, 'H', imgH, self.roiH[2], fitting.NoFit(self.imaging_pars['H']), pylab.get_cmap(self.imaging_pars['H'].palette))

        #TODO: better: reduce constructor, call later show_image
        
        self.mgr.AddPane(self.V, wx.aui.
                         AuiPaneInfo().Name("V").Caption("Vertical").
                         Top().Position(1).BestSize(wx.Size(640, 480))
                         )
        self.mgr.AddPane(self.H, wx.aui.
                         AuiPaneInfo().Name("H").Caption("Horizontal").
                         Top().Position(2)
                         )
        
        ## create status bar
        
        self.statusbar = StatusBarWx(self.frame)
        self.V.toolbar.set_status_bar(self.statusbar)
        self.H.toolbar.set_status_bar(self.statusbar)
        self.frame.SetStatusBar(self.statusbar)
        
        ## create toolbars top
        self.tb = wx.ToolBar(self.frame, - 1, wx.DefaultPosition,
                             wx.DefaultSize,
                             wx.TB_FLAT | wx.TB_NODIVIDER
                             )
        
        self.tb.SetToolBitmapSize((24, 24)) #TODO: ???

        self.tb2 = wx.ToolBar(self.frame, - 1, wx.DefaultPosition,
                             wx.DefaultSize,
                             wx.TB_FLAT | wx.TB_NODIVIDER
                             )
        
        self.tb2.SetToolBitmapSize((24, 24)) #TODO: ???
        
        self.tb3 = wx.ToolBar(self.frame, - 1, wx.DefaultPosition,
                             wx.DefaultSize,
                             wx.TB_FLAT | wx.TB_NODIVIDER
                             )
        
        self.tb3.SetToolBitmapSize((24, 24)) #TODO: ???
        
        self.tb4 = wx.ToolBar(self.frame, - 1, wx.DefaultPosition,
                             wx.DefaultSize,
                             wx.TB_FLAT | wx.TB_NODIVIDER
                             )
        
        self.tb4.SetToolBitmapSize((24, 24)) #TODO: ???
        
        # measurement
        #self.measurement = wx.TextCtrl(self.tb, -1, "", size = (100,-1))
        self.measurement_combobox = wx.ComboBox(self.tb, - 1, "", size=(100, - 1), style=wx.TE_PROCESS_ENTER)
        self.tb.AddControl(self.measurement_combobox)

        self.frame.Bind(wx.EVT_TEXT_ENTER, self.OnChangeMeasurementName, self.measurement_combobox)
        self.frame.Bind(wx.EVT_COMBOBOX, self.OnChangeMeasurementName, self.measurement_combobox)
        self.tb.AddSeparator()

        # Save Image
        self.savebutton = wx.Button(self.tb, self.ID_SaveImage,
                                    "Save Image", size=(100, 30)
                                    )
        self.savebutton.SetBackgroundStyle(wx.BG_STYLE_COLOUR)
        self.tb.AddControl(self.savebutton)
        self.frame.Bind(wx.EVT_BUTTON, self.OnSaveImage, id=self.ID_SaveImage)

        ## Automatic Save Image
        self.autosave_control = wx.CheckBox(self.tb, self.ID_Autosave,
                                    "Auto Save Image")
        self.tb.AddControl(self.autosave_control)
        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnAutosaveToggle,
                        id=self.ID_Autosave)
        self.autosave = False

        # Reload Button
        self.reloadbutton = wx.Button(self.tb, self.ID_ReloadButton,
                                  "Reload", size=(100, 30)
                                  )
        self.tb.AddControl(self.reloadbutton)
        self.frame.Bind(wx.EVT_BUTTON, self.OnReloadButton, id=self.ID_ReloadButton)

        # Automatic reload
        self.autoreload = wx.CheckBox(self.tb, self.ID_Autoreload,
                                      "Auto Reload")
        self.tb.AddControl(self.autoreload)
        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnAutoreloadClicked,
                        id=self.ID_Autoreload)
        
        # Added when changing filewatch
        self.FileWatcher = filewatch.FileChangeNotifier(self.watchedfiles,
                                                        self.ConcatenateSis)

        # Record data
        self.record_data_button = wx.Button(self.tb, self.ID_RecordData,
                                     "Record Data", size=(100, 30)
                                     )
        self.record_data = True
        self.record_data_button.SetBackgroundStyle(wx.BG_STYLE_COLOUR)
        self.record_data_button.SetBackgroundColour(wx.NamedColour('GREEN'))
        self.tb.AddControl(self.record_data_button)
        self.frame.Bind(wx.EVT_BUTTON,
                        self.OnRecordDataButtonClicked,
                        id=self.ID_RecordData)

        ##save image filename display
        #self.tb.AddSeparator()
        #self.saved_image_name = wx.StaticText(self.tb,
        #                                      -1,
        #                                      'image not saved',
        #                                      size = (200, -1),
        #                                      style = wx.ALIGN_LEFT)
        #self.tb.AddControl(self.saved_image_name)

        # Imaging Parameters:
        # Imaging direction: Choicebox
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                         'Imaging: '))
        choices = []
        for p in self.imaging_parlist:
            choices.append(p['V'].description+" | "+p['H'].description)
        
        self.imaging_choice = wx.Choice(self.tb2,
                                        self.ID_ImagingChoice,
                                        size=(120, - 1),
                                        choices=choices)
        #self.imaging_choice.Select(0)
        #self.select_imaging(0)
        self.tb2.AddControl(self.imaging_choice)
        self.Bind(wx.EVT_CHOICE, self.OnChoiceImaging, self.imaging_choice)



        choices=["cam0-cam1","cam0-cam0","cam1-cam1"]
        self.imaging_AB =wx.Choice(self.tb2,self.ID_Imaging_AB,size=(80,-1),choices=choices)
        self.tb2.AddControl(self.imaging_AB)
        self.Bind(wx.EVT_CHOICE,self.OnChoice_AB,self.imaging_AB)



        # Expansion times
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " Expansion time V: "))
        self.expansion_time_V = wx.SpinCtrl(self.tb2,
                                            self.ID_ExpansionTimeV,
                                            min=0,
                                            max=500,
                                            initial=10,
                                            name="expansion time V",
                                            size=(50, - 1))
        self.tb2.AddControl(self.expansion_time_V)

        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " H: "))
        self.expansion_time_H = wx.SpinCtrl(self.tb2,
                                             self.ID_ExpansionTimeH,
                                             min=0,
                                             max=500,
                                             initial=10,
                                             name="expansion time H",
                                             size=(50, - 1))
        self.tb2.AddControl(self.expansion_time_H)

        self.Bind(wx.EVT_SPINCTRL, self.OnChangeExpansionTime, self.expansion_time_V)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeExpansionTime, self.expansion_time_V)
        self.Bind(wx.EVT_SPINCTRL, self.OnChangeExpansionTime, self.expansion_time_H)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeExpansionTime, self.expansion_time_H)                                           
        #maximum optical density entries
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " Max. OD V: "))
        self.entry_optical_density_max_V = wx.TextCtrl(self.tb2,
                                                 self.ID_OpticalDensityMaxV,
                                                 '0.0',
                                                 size=(50, - 1),
                                                 name="max. OD V",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb2.AddControl(self.entry_optical_density_max_V)
        
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " H: "))
        self.entry_optical_density_max_H = wx.TextCtrl(self.tb2,
                                                 self.ID_OpticalDensityMaxH,
                                                 '0.0',
                                                 size=(50, - 1),
                                                 name="max. OD H",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb2.AddControl(self.entry_optical_density_max_H)
        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeMaximumOpticalDensity, self.entry_optical_density_max_H)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeMaximumOpticalDensity, self.entry_optical_density_max_V)
                


        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " Rot V: "))
        self.rotation_V = wx.TextCtrl(self.tb2,
                                                 self.ID_RotateV,
                                                 '0.0',
                                                 size=(50, - 1),
                                                 name="Rotation V",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb2.AddControl(self.rotation_V)
        
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " H: "))
        self.rotation_H = wx.TextCtrl(self.tb2,
                                                 self.ID_RotateH,
                                                 '0.0',
                                                 size=(50, - 1),
                                                 name="Rotation H",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb2.AddControl(self.rotation_H)
        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeRotation, self.rotation_H)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeRotation, self.rotation_V)





        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " Vmin V: "))
        self.vmin_V = wx.TextCtrl(self.tb2,
                                                 self.ID_VminV,
                                                 '-0.05',
                                                 size=(50, - 1),
                                                 name="Vmin V",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb2.AddControl(self.vmin_V)
        
        self.tb2.AddControl(wx.StaticText(self.tb2,
                                          -1,
                                          " H: "))
        self.vmin_H = wx.TextCtrl(self.tb2,
                                                 self.ID_VminH,
                                                 '-0.05',
                                                 size=(50, - 1),
                                                 name="Vmin H",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb2.AddControl(self.vmin_H)
        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeVmin, self.vmin_H)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeVmin, self.vmin_V)










        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          "V residuals:"))
        self.differenceV = wx.CheckBox(self.tb3, self.ID_DifferenceV,
                                      "")
        
        self.tb3.AddControl(self.differenceV)
        
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " filter:"))
        self.filterV = wx.CheckBox(self.tb3, self.ID_FilterV,
                                      "")
        self.tb3.AddControl(self.filterV)
        
        
        
        
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " filter sigma: "))
        self.entry_filter_sigma_V = wx.TextCtrl(self.tb3,
                                                 self.ID_FilterSigmaV,
                                                 str(self.V.filter_sigma),
                                                 size=(50, - 1),
                                                 name="sigma V",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb3.AddControl(self.entry_filter_sigma_V)
        
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " scale: "))
        self.entry_filter_scalar_V = wx.TextCtrl(self.tb3,
                                                 self.ID_FilterScalarV,
                                                 str(self.V.filter_scalar),
                                                 size=(50, - 1),
                                                 name="scale V",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb3.AddControl(self.entry_filter_scalar_V)
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " reference: "))
        self.referenceV = wx.CheckBox(self.tb3, self.ID_ReferenceV,
                                      "")
        self.tb3.AddControl(self.referenceV)
        
        
        
        self.tb3.AddSeparator()
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " H residuals: "))
        self.differenceH = wx.CheckBox(self.tb3, self.ID_DifferenceH,
                                      "")

        self.tb3.AddControl(self.differenceH)
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " filter: "))
        self.filterH = wx.CheckBox(self.tb3, self.ID_FilterH,
                                      "")
        self.tb3.AddControl(self.filterH)
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " filter sigma: "))
        self.entry_filter_sigma_H = wx.TextCtrl(self.tb3,
                                                 self.ID_FilterSigmaH,
                                                 str(self.H.filter_sigma),
                                                 size=(50, - 1),
                                                 name="sigma H",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb3.AddControl(self.entry_filter_sigma_H)
        
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " scale: "))
        self.entry_filter_scalar_H = wx.TextCtrl(self.tb3,
                                                 self.ID_FilterScalarH,
                                                 str(self.H.filter_scalar),
                                                 size=(50, - 1),
                                                 name="scale H",
                                                 style=wx.TE_PROCESS_ENTER,
                                                 )
        self.tb3.AddControl(self.entry_filter_scalar_H)
        
        self.tb3.AddControl(wx.StaticText(self.tb3,
                                          -1,
                                          " reference: "))
        self.referenceH = wx.CheckBox(self.tb3, self.ID_ReferenceH,
                                      "")
        self.tb3.AddControl(self.referenceH)
        
        
        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnDifferenceSelectV,
                        id=self.ID_DifferenceV)
        
        
        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnDifferenceSelectH,
                        id=self.ID_DifferenceH)

        
        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnFilterSelectV,
                        id=self.ID_FilterV)

        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnFilterSelectH,
                        id=self.ID_FilterH)


        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnReferenceSelect,
                        id=self.ID_ReferenceV)


        self.frame.Bind(wx.EVT_CHECKBOX,
                        self.OnReferenceSelect,
                        id=self.ID_ReferenceH)

        

        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeFilterSigma, self.entry_filter_sigma_H)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeFilterSigma, self.entry_filter_sigma_V)
        
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeFilterScalar, self.entry_filter_scalar_H)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnChangeFilterScalar, self.entry_filter_scalar_V)
                



        
        self.savebitmapbutton = wx.Button(self.tb4, self.ID_SaveBitmap,
                                    "Save Bitmaps", size=(100, 30)
                                    )
        self.savebitmapbutton.SetBackgroundStyle(wx.BG_STYLE_COLOUR)
        self.tb4.AddControl(self.savebitmapbutton)
        self.frame.Bind(wx.EVT_BUTTON, self.OnSaveBitmap, id=self.ID_SaveBitmap)

        
        self.tb4.AddControl(wx.StaticText(self.tb4,
                                          -1,
                                         ' Colors V: '))
        
        self.maps=[map for map in pylab.cm.datad if not map.endswith("_r")]
        self.maps.sort()
        
        self.cmap_choiceV = wx.Choice(self.tb4,
                                        self.ID_CMapChoiceV,
                                        size=(75, - 1),
                                        choices=self.maps)
        self.tb4.AddControl(self.cmap_choiceV)
        self.Bind(wx.EVT_CHOICE, self.OnChoiceCMap, self.cmap_choiceV)
        
        self.tb4.AddControl(wx.StaticText(self.tb4,
                                          -1,
                                         ' H: '))
        self.cmap_choiceH = wx.Choice(self.tb4,
                                        self.ID_CMapChoiceH,
                                        size=(75, - 1),
                                        choices=self.maps)
        self.tb4.AddControl(self.cmap_choiceH)
        self.Bind(wx.EVT_CHOICE, self.OnChoiceCMap, self.cmap_choiceH)

        self.createreference = wx.Button(self.tb4, self.ID_CreateReference,
                                  "Create reference", size=(100, 30)
                                  )
        self.tb4.AddControl(self.createreference)
        self.frame.Bind(wx.EVT_BUTTON, self.OnCreateReference, id=self.ID_CreateReference)



        #finalize toolbars
        self.tb.Realize()
        self.tb2.Realize()
        self.tb3.Realize()
        self.tb4.Realize()
        
        self.mgr.AddPane(self.tb,
                         wx.aui.AuiPaneInfo().
                         Name("toolbar").
                         ToolbarPane().Top().Row(1).Position(1))


        self.mgr.AddPane(self.tb2,
                         wx.aui.AuiPaneInfo().
                         Name("toolbar2").
                         ToolbarPane().Top().Row(1).Position(2))
        self.mgr.AddPane(self.tb3,
                         wx.aui.AuiPaneInfo().
                         Name("toolbar3").
                         ToolbarPane().Top().Row(2).Position(1))
        self.mgr.AddPane(self.tb4,
                         wx.aui.AuiPaneInfo().
                         Name("toolbar4").
                         ToolbarPane().Top().Row(2).Position(2))

        #ImageTree
        self.savedimagestree = ImageTree.TreeModelSync(os.path.join(settings.imagesavepath, '')) #TODO
        self.imagetreepanel = ImageTree.ImageTreePanel(self.frame,
                                             self.savedimagestree)
        self.selectedsavedimage = None
        self.mgr.AddPane(self.imagetreepanel,
                         wx.aui.AuiPaneInfo().
                         Name("imagetree").
                         Caption("Saved Images").
                         Right().BestSize(wx.Size(200, 400))
                         )
        self.imagetreepanel.Bind(wx.EVT_TREE_ITEM_ACTIVATED, self.OnActivateImageTree)
        self.imagetreepanel.Bind(wx.EVT_TREE_ITEM_RIGHT_CLICK, self.OnRightClickImageTree)
        

        #ImagePanels for saved images
        self.savedimageV = ImagePanel.CamImageMarkersPanel(self.frame)
        self.savedimageH = ImagePanel.CamImageMarkersPanel(self.frame)

        for panel in [self.savedimageV, self.savedimageH]:
            panel.toolbar.AddSeparator()
            panel.toolbar.AddLabelTool(self.ID_AnalyzeImageButton,
                                       'Analyze',
                                       wx.Bitmap(os.path.join(settings.bitmappath,
                                                              'cam24.png'),
                                                 wx.BITMAP_TYPE_PNG),
                                       shortHelp = 'Analyze',
                                       longHelp = 'load image in cam for analyzing')
            panel.toolbar.Realize()
            panel.Bind(wx.EVT_TOOL, self.OnAnalyzeImage, id=self.ID_AnalyzeImageButton)

        self.mgr.AddPane(self.savedimageV,
                         wx.aui.AuiPaneInfo().
                         Name("savedimageV").
                         Caption("saved image of V").
                         Float().Dockable(False).
                         BestSize(wx.Size(400, 300))
                         )
        self.mgr.AddPane(self.savedimageH,
                         wx.aui.AuiPaneInfo().
                         Name("savedimageH").
                         Caption("saved image of H").
                         Float().Dockable(False).
                         BestSize(wx.Size(400, 300))
                         )
        
        ## create some perspectives
        self.mgr.GetPane("V").Show()
        self.mgr.GetPane("H").Show()
        self.mgr.GetPane("results").Show()
        self.mgr.GetPane("Toolbar").Show()
        self.mgr.GetPane("savedimageH").Hide()
        self.mgr.GetPane("savedimageV").Hide()
        self.perspective_full = self.mgr.SavePerspective()

        self.mgr.GetPane("V").Hide()
        self.mgr.GetPane("H").Show().Left()
        self.perspective_H = self.mgr.SavePerspective()

        self.mgr.GetPane("H").Hide()
        self.mgr.GetPane("V").Show().Left()
        self.perspective_V = self.mgr.SavePerspective()

        self.mgr.LoadPerspective(self.perspective_full)

        ## Reload Event
        self.frame.Bind(EVT_RELOAD_IMAGE, self.OnReloadEvent)

        ## Save Image Event
        self.frame.Bind(EVT_SAVE_IMAGE, self.OnSaveImageEvent)

        ##Fit update event
        self.frame.Bind(EVT_FIT_RESULTS, self.OnFitResultsEvent)
        
        #tweak scaling
        self.H.axhprof.set_xlim(self.H.roi.xmin - 100, self.H.roi.xmax + 100)
        self.H.toolbar.push_current()

        ##initialize








        ######default settings 2013-12-03######################################
        fit_menu.Check(self.ID_FitHGauss, True)
        fit_menu.Check(self.ID_FitVGauss, True)
        fit_menu.Check(self.ID_CompensateOffsetH, True)
        fit_menu.Check(self.ID_CompensateOffsetV, True)
        
        self.H.fit = fitting.Gauss2d(self.imaging_pars['H'])
        self.H.update()
        self.V.fit = fitting.Gauss2d(self.imaging_pars['V'])
        self.V.update()
        self.H.set_compensate_offset(True)
        self.V.set_compensate_offset(True)
        self.H.OnAutoscale(None)
        self.V.OnAutoscale(None)
        
        self.autoreload.SetValue(True)
        self.OnAutoreloadClicked(self.autoreload)

        
        self.plotpanels = dict()
        self.ding = ding.Ding() #CLASS DING TO BE DONE
        self.InitSettings()

        self.mgr.Update()
        self.frame.Show(True)

        return True


############# added 17-12-2012
    def ConcatenateSis(self):
        #self.imagefilename
        
        img0 = read(self.imagefilename.replace(".sis","_0.sis"))
        img1 = read(self.imagefilename.replace(".sis","_1.sis"))
        h, w = img0.shape

        img=numpy.concatenate((img0[h/2:], img1[h/2:]))

        img.shape=img0.shape

        #img=numpy.concatenate([img0, img1])
                
        write_raw_image(self.imagefilename, img)
        self.CreateReloadEvent()
        
##############


    def InitSettings(self):
        """Initialize all values and settings which is not yet possible to do 
        in OnInit since some of the controls do not yet exist"""
        self.H.show_contours = True
        self.V.show_contours = True
        
        self.cmap_choiceV.SetStringSelection(self.imaging_pars['V'].palette)
        self.cmap_choiceH.SetStringSelection(self.imaging_pars['H'].palette)
        self.imaging_choice.Select(0)
        self.select_imaging(0)
        self.imaging_AB.Select(0)

        
        self.measurements = []
        self.select_measurement("results")
        self.sync_imaging_pars_expansion_time()
        imgV, imgH = loadimg(self.imagefilename)
        ############# added 17-12-2012
#        self.ConcatenateSis()
#        self.FileWatcher = filewatch.FileChangeNotifier(self.imagefilename.replace(".sis", "_0.sis"),
#                                                        self.ConcatenateSis)
        
#        self.fileconcatenatethread0 = filewatch.FileChangeNotifier(self.imagefilename.replace(".sis","_0.sis"),
#                                                    callback=self.ConcatenateSis)
##        self.fileconcatenatethread1 = filewatch.FileChangeNotifier(self.imagefilename.replace(".sis","_1.sis"),
#                                                    callback=self.ConcatenateSis)
#        self.fileconcatenatethread0.start()
        #self.fileconcatenatethread1.start()
        #############
        
        
        
    def UpdateResults(self):
        self.results.UpdateResults({'V': self.V.fitpars, 'H': self.H.fitpars})

    def UpdateResultsFilename(self, filename):
        """Update entry filename in results (after saving or reloading from file)"""
        head, tail = os.path.split(filename)
        self.results.UpdateFilename(tail)

    def OnReloadEvent(self, event):
        """Reload image. load image (default location if
        event.filename == None), store results into event.target"""

        if event.filename is not None:
            full_path = self.savedimagestree.find_file(event.filename)
            if full_path:
                filename = full_path
            else:
                print "can't find image!"
                return
        else:
            filename = self.imagefilename
            
        self.results.activate() #activate results table on reload
        self.load_image(filename, event.target)
        self.savebutton.SetBackgroundColour(wx.NamedColour("RED"))
        self.savebutton.Refresh() #necessary?
        
        #self.saved_image_name.SetLabel('image not saved') #TODO: ??
        
        self.statusbar.StatusText = 'image not saved'

        if (self.autosave_abs or self.autosave_raw) and event.target is None:
            wx.PostEvent(self.frame, SaveImageEvent())
        event.Skip()

    def load_image(self, filename, target = None):

        """load image, analyze image. Store results in target. Note
        results are created asynchronously in parallel threads,
        received by FitResultEvent."""
        #imgH, imgH = loadimg(filename)   #####added 14-12-2012
        #imgV, imgV = loadimg(filename.replace(".sis","_1.sis"))   #####added 14-12-2012
        imgV, imgH = loadimg(filename)   #####original 14-12-2012

        imgV_old=imgV
        imgH_old=imgH
        
        if self.imaging_AB.GetSelection()==0:
            imgV=imgV_old
            imgH=imgH_old
        if self.imaging_AB.GetSelection()==1:
            imgV=imgV_old
            imgH=imgV_old
        if self.imaging_AB.GetSelection()==2:
            imgV=imgH_old
            imgH=imgH_old
        
        #loading new image also means appending row
        if self.record_data and target is None:
            self.results.AppendRows(1) #TODO: this should be done somewhere else

        self.V.analyze_image(imgV, target)
        self.H.analyze_image(imgH, target)

        
    def OnFitResultsEvent(self, event):
        species = event.source
        target = event.target
        fitpars = event.fitpars
        if target is None:
            self.results.UpdateResults({species: fitpars})
        else:
            if target['name'] != self.results.name:
                print "target measurement '%s' does not match active measurement '%s'"%(target['name'], self.results.name)
            else:
                self.results.UpdateResults(data = {species: fitpars},
                                           row = target['row'])

        event.Skip()
        self.UpdateResiduals()



    def OnAutoreloadClicked(self, event):
        if event.IsChecked():
            #activate automatic reload
            self.FileWatcher.setEnabled(True)
#            self.filewatchthread = filewatch.FileChangeNotifier(self.imagefilename,    ############commented on 01-02-2013
#                                                callback=self.CreateReloadEvent)
#            self.filewatchthread.start()
        else:
            self.FileWatcher.setEnabled(False)
#            self.filewatchthread.keeprunning = False
#            self.filewatchthread.join()
            

    def OnRecordDataButtonClicked(self, event):
        if self.record_data:
            self.record_data_button.SetBackgroundColour(wx.NamedColour("RED"))
            #self.record_data_button.SetForegroundColour(wx.NamedColour("RED"))
            self.record_data = False
            self.results.record = False
        else:
            self.record_data_button.SetBackgroundColour(wx.NamedColour("GREEN"))
            #self.record_data_button.SetForegroundColour(wx.NamedColour("GREEN"))
            self.record_data = True
            self.results.record = True

        self.record_data_button.Refresh()
        self.gridpanel.Refresh()

    def PostReloadImageEvent(self):
        """post an ReloadImageEvent with empty arguments: load default
        image, append results to active measurement""" 
        wx.PostEvent(self.frame, ReloadImageEvent(filename = None, target = None))

    def AcousticSignal(self):
        self.ding.play()

    def CreateReloadEvent(self):
        self.AcousticSignal()
        self.PostReloadImageEvent()

    def OnMenuSettingsSave(self, event):
        import shelve

        S = shelve.open(os.path.join(settings.basedir,
                                            'settings/settings_cam'))
        try:
            #Window layout
            S['perspective'] = self.mgr.SavePerspective()

            #regions of interest
            S['roiV'] = self.roiV
            S['roiH'] = self.roiH

            #imaging pars
            S['imaging_selection'] = self.imaging_choice.GetSelection()

            #markers
            markersVpositions = []
            for k, marker in enumerate(self.V.markers):
                markersVpositions.append(marker.position)

            markersHpositions = []
            for k, marker in enumerate(self.H.markers):
                markersHpositions.append(marker.position)

            S['markersVpositions'] = markersVpositions
            S['markersHpositions'] = markersHpositions

            S['imaging_parlist'] = self.imaging_parlist
            S['imaging_choice'] = self.imaging_choice.Selection
            
                
        finally:
            S.close()
        
    def OnMenuSettingsLoad(self, event):
        import shelve
        #with shelve.open('settings') as settings:
        S = shelve.open(os.path.join(settings.basedir,
                                            'settings/settings_cam'))
        try:
            #perspectives
            self.mgr.LoadPerspective(S['perspective'])

            #region of interests
            self.roiV = S['roiV']
            self.roiH = S['roiH']

            #imaging settings
            imaging_sel = S['imaging_selection']
            self.imaging_choice.SetSelection(imaging_sel)
            self.select_imaging(imaging_sel)
            self.select_roi(imaging_sel)

            #markers
            mVpos = S['markersVpositions']
            for k, pos in enumerate(mVpos):
                self.V.markers[k].position = pos

            mHpos = S['markersHpositions']
            for k, pos in enumerate(mHpos):
                self.H.markers[k].position = pos
                
            self.imaging_parlist = S['imaging_parlist']
            imagingsel = S['imaging_choice'] 
            
            self.imaging_choice.Selection = imagingsel
            self.select_roi(imagingsel)
            self.imaging_parlist = S['imaging_parlist']
            
            self.select_imaging(imagingsel)
             
        finally:
            S.close()
        
        self.V.center_roi()
        self.H.center_roi()
            
        self.V.update()
        self.H.update()

    

    #def reload(self, mplevent):
    #    imgV, imgH = loadimg(self.filename)
    #    if mplevent.inaxes == self.H.axtop:
    #        self.H.show_image(imgH)
    #        self.V.show_image(imgV)
    #    elif mplevent.inaxes == self.V.axtop:
    #        self.V.show_image(imgV)
    #        self.H.show_image(imgH)
    #    else:
    #        self.V.show_image(imgV)
    #        self.H.show_image(imgH)
    #    #TODO: share following code with OnReloadEvent ?
    #    self.results.AppendRows()
    #    self.UpdateResults()
    #    self.savebutton.SetBackgroundColour(wx.HmedColour("RED"))


    def OnChangeMeasurementName(self, event):
        name = event.GetString()
        self.select_measurement(name)

    def new_measurement(self, name):
        """create notebook panel, add name to
         combo box list and set combo box
         value
        """
        self.gridpanel.addpage(name)
        self.measurement_combobox.Append(name)
        self.measurement_combobox.Value = name
        self.measurements.append(name)
        
    def delete_measurement(self, name):
        #remove from combobox
        index = self.measurements.index(name)
        self.measurement_combobox.Delete(index)
        #remove from internal list
        self.measurements.remove(name)

        #remove plotpanel
        panelname = name + ' plot' #see plot_results
        plotpane = self.mgr.GetPane(panelname)
        if plotpane.IsOk():
            datapanel = self.plotpanels.pop(panelname)
            self.mgr.DetachPane(datapanel)
            #del datapanel
        
    def select_measurement(self, name):
        """select measurement. if measurement does not exist, add new
        page, update control.
        """
        if self.measurement_combobox.FindString(name) == wx.NOT_FOUND:
            self.new_measurement(name)
            
        index = self.measurements.index(name)
        self.measurement_combobox.SetSelection(index)
        self.gridpanel.notebook.SetSelection(index)

        #deactivate Table
        try:
            self.results.activate(False) #will fail if self.results is
                                         #not yet initialized
        except AttributeError:
            pass 

        #change short reference 'results' and 
        self.results = self.gridpanel.pages[index].Table
        self.results.name = name
        #self.results.activate() #activation done on first reload

        #mark active page, first remove mark
        for i in range(self.gridpanel.notebook.GetPageCount()):
            self.gridpanel.notebook.SetPageBitmap(i, wx.NullBitmap)
        self.gridpanel.notebook.SetPageBitmap(index,
                                              wx.Bitmap(
                                                  os.path.join(settings.bitmappath,
                                                               'star-small.png'),
                                                  wx.BITMAP_TYPE_PNG)
                                              )
        self.gridpanel.Refresh()

    def create_plotpanel(self, measurement_name, panelname):
        """create plotpanel, add to AUI Manager, establish connections
        for automatic updates, store it."""

        plotpanel = DataPanel.DataPlotPanel(self.frame, measurement_name)
        self.mgr.AddPane(plotpanel,
                         wx.aui.AuiPaneInfo()
                         .Name(panelname)
                         .Caption(measurement_name)
                         .BestSize(wx.Size(400, 300))
                         .Float()
                         .Dockable(False)
                         )
        self.mgr.Update()

        #plotpanel.set_table(self.results)
        self.results.add_observer(plotpanel)
        return plotpanel

    def OnMenuPlotResults(self, event):
        """Handle selection of 'plot results' menu. Calls L{create_plotpanel}
        """
        self.plot_results()

    def plot_results(self):
        """Create plot panel
        if necessary, show if hidden, update plotpanel.
        """

        measurement_name = self.results.name
        panelname = measurement_name + ' plot'
        
        #try to get plot panel (AUIPaneInfo) associated with current measurement
        plotpane = self.mgr.GetPane(panelname)
        if not plotpane.IsOk(): #plot panel does not yet exist
            self.plotpanels[panelname] = self.create_plotpanel(measurement_name, panelname)
            plotpane = self.mgr.GetPane(panelname)

        if not plotpane.IsShown():    
            plotpane.Show()
            self.mgr.Update()
            self.plotpanels[panelname].refresh()
            
        #TODO: plotpanel skips updates if not visible
        #self.plotpanels[name].update()
        self.plotpanels[panelname].draw()
        

    def OnActivateImageTree(self, event):
        index = self.imagetreepanel.treectrl.GetIndexOfItem(event.Item)
        filename = self.savedimagestree.GetItemFile(index)

        if len(index)<3:
            #not on image
            return

        if self.savedimagestree.check_day(index[0]):
            print "refresh day"
            self.imagetreepanel.treectrl.RefreshItems()

        self.displayedsavedimage = filename

        #imgV, imgV = ImagePanel.loadimg(filename)       ######added 14-12-2012
        #imgH, imgH = ImagePanel.loadimg(filename.replace(".sis","_1.sis"))       ######added 14-12-2012
        imgV, imgH = ImagePanel.loadimg(filename)       ######original 14-12-2012        
        self.savedimageV.show_image(imgV, description=filename)
        self.savedimageH.show_image(imgH, description=filename)

        paneimgna = self.mgr.GetPane('savedimageH')
        paneimgk = self.mgr.GetPane('savedimageV')
        if not paneimgna.IsShown() and not paneimgk.IsShown():
            paneimgna.Float().Show()
            paneimgk.Float().Show()
            self.mgr.Update()

        next = self.imagetreepanel.treectrl.GetNextSibling(event.Item)
        if next.IsOk():
            self.imagetreepanel.treectrl.SelectItem(next)

    def OnRightClickImageTree(self, event):
        index = self.imagetreepanel.treectrl.GetIndexOfItem(event.Item)
        if len(index)>=3:
            filename = self.savedimagestree.GetItemFile(index)
            self.selectedsavedimage = filename
        else:
            self.selectedsavedimage = None

        menu = wx.Menu()
        menu.Append(self.ID_ImageTreeRescan, 'Rescan')
        if len(index)>=1:
            menu.Append(-1, "Day Menu")
        if len(index)>=2:
            menu.Append(-1, "Measurement Menu")
        if len(index)>=3:
            menu.Append(self.ID_ImageTreeReloadSavedImage, "Reload Image")
        self.imagetreepanel.Bind(wx.EVT_MENU, self.OnReloadSavedImage, id = self.ID_ImageTreeReloadSavedImage)
        self.imagetreepanel.Bind(wx.EVT_MENU, self.OnRescanImageTree, id = self.ID_ImageTreeRescan)
        self.imagetreepanel.PopupMenu(menu)
        menu.Destroy()

    def OnAnalyzeImage(self, event):
        """Load and analyze images. Bound to button in saved image panels"""
        self.load_image(self.displayedsavedimage)
        self.UpdateResultsFilename(self.displayedsavedimage)
        self.statusbar.StatusText = 'reload: ' + self.displayedsavedimage

    def OnReloadSavedImage(self, event):
        """Load and analyze images. Bound to popup menu in image tree"""
        self.load_image(self.selectedsavedimage)
        self.UpdateResultsFilename(self.selectedsavedimage)
        self.statusbar.StatusText = 'reload: ' + self.selectedsavedimage

    def OnRescanImageTree(self, event):
        print "rescan files", self.imagetreepanel.treemodel.root
        self.imagetreepanel.treemodel.createfiletree()
        self.savedimagestree = self.imagetreepanel.treemodel
        self.imagetreepanel.treectrl.RefreshItems()
        
    def OnResultsPageClose(self, event):
        page = self.gridpanel.notebook.GetPage(event.Selection)
        results = page.Table
        need_save = results.modified
        name = results.name
        
        if need_save:
            dlg = wx.MessageDialog(self.frame,
                               "Data changes are not yet saved.\n"+\
                                "Do you want to save data?\n"+\
                                "(otherwise they are lost, forever!)",
                               caption = "Close data window",
                               style = wx.YES_NO | wx.CANCEL | wx.ICON_EXCLAMATION,
                               )
            answer = dlg.ShowModal()
            dlg.Destroy()
            if answer == wx.ID_YES:
                success = self.results_save(filename = results.filename, results = results, forcesave = True)
                if success:
                    self.delete_measurement(name)
                else:
                    event.Veto()
            elif answer == wx.ID_NO:
                print "don't save data, it's your choice"
                self.delete_measurement(name)
            else: #CANCEL
                event.Veto()
        else:
            self.delete_measurement(name)
        
    def OnResultsPageClosed(self, event):
        self.select_measurement("results")
        print "OnResultsPageClosed"
        
"""
def test_img_from_disk():

    d = shelve.open('sis.pkl')

    print d.keys()
    img1 = d['img1']
    img2 = d['img2']
    bkg = d['img3']
    d.close()

    img1 = (img1 - bkg).astype(Float)
    img2 = (img2 - bkg).astype(Float)

    den = - log(img1 / img2)

    imgV = den[0:1040, :]
    imgH = den[1040:2080, :]

    img = ma.array(imgH, mask= ~ isfinite(imgH))
    img.set_fill_value(0)

    #img = ma.masked_inside(imgH, -0.1, 3)


    roi1 = ROI(450, 650, 350, 480)
    roi2 = ROI(450, 650, 350, 481)

    imgdisp = ImgDisp(img, roi1, fitting.Fitting(), fig=1)
    #imgdisp2 = ImgDisp(ma.array(imgV, mask = ~isfinite(imgV)), roi2, Fitting(), fig=2)


    show()
"""

def run_cam():
    gui = ImgAppAui(redirect=False)
    gui.MainLoop()
    return gui

"""
def profile():
    import cProfile
    cProfile.run('run()', stats)
"""    
if __name__ == '__main__':
    gui = run_cam()
    
