import os

#basedir = r"c:/fit/siscam"
basedir = "/home/gabriele/Desktop/CAMERA/CAM/"

#full path to image file

#imagefile = os.path.join(basedir, 'img/test.sis')

#WHEN USING IMAGING-PC ON WINDOWS
imagefile = os.path.join('/home/gabriele/sis-fish/', 'test.sis')

#WHEN USING IMAGING-PC ON UBUNTU
#imagefile = os.path.join('sftp://bec1@10.194.32.22/home/bec1/sis-fish', 'test.sis')


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
