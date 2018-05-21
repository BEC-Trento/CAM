import os
from sys import path

#basedir = r"c:/fit/siscam"
#basedir = "/home/gabriele/Desktop/CAMERA/CAM/"
basedir = path[0]


### Full path to image file
#imagefile = os.path.join(basedir, 'img/test.sis')
#imagefile = os.path.join('/home/gabriele/sis-fish/', 'test.sis')
imagefile = os.path.join(basedir, 'test.sis') # the combined .sis that will be saved

### Full path to incoming images files
watchedfiles = ['test_0.sis', 'test_1.sis'] # files that are monitored for autoreload

#are they now deprecated?
rawimage1file = imagefile
rawimage2file = imagefile
rawimage3file = imagefile

referencefile = os.path.join(basedir, 'img/reference.sis')

#where to save images
#imagesavepath = r"./img"
#imagesavepath = os.path.join(basedir, 'img')

imagesavepath = '/backup/img'

#icons, etc.
bitmappath = os.path.join(basedir, 'bitmaps')

#directory to store template files
templatedir = os.path.join(basedir, 'templates')

##acquire
useTheta = True
#useTheta = False
useBluefox = True
#useBluefox = False
useSony = False

#settings for Theta Systems cam
#configfile = r"c:\WinSIS6\py\config.ini"
configfile = r"./WinSIS6/py/config.ini"
#usePseudoCam = False
usePseudoCam = True
