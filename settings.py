### TEMPLATE FILE FOR CAM SETTINGS

import os

#basedir = r"c:/fit/siscam"
#basedir = "/home/gabriele/Desktop/CAMERA/CAM/"
basedir = "/home/arturo/Documents/CAM"

### Full path to image file
#imagefile = os.path.join(basedir, 'img/test.sis')

#Imaging PC writing on host
imagefile = os.path.join(basedir, 'test.sis')

#Reading from imaging PC
#imagefile = os.path.join('sftp://bec1@10.194.32.22/home/bec1/sis-fish', 'test.sis')

### Full path to incoming image-files
watchedfiles = [os.path.join(basedir, 'test_0.sis'),
                os.path.join(basedir, 'test_1.sis')]

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
