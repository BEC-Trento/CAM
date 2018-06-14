#!/usr/bin/python
#-*- coding: latin-1 -*-
"""Watch if a file is changed. Used for automatic reloading."""

#import threading
#import os
#import time
#import sys

#
# def dosomething():
#    print "filechange!"
#    sys.stdout.flush()
#
#
# class FileChangeNotifier(threading.Thread):
#
#    def __init__(self, filename, callback=dosomething, delay=2):
#        threading.Thread.__init__(self)
#        self.filename = filename
#        self.callback = callback
#        self.delay = delay
#
#        self.s = os.stat(self.filename)
#        self.keeprunning = True
#
#    def run(self):
#        while self.keeprunning:
#            s = os.stat(self.filename)
#            if s.st_mtime > self.s.st_mtime:
#                self.s = os.stat(self.filename)
#                time.sleep(self.delay)
#                # horrible hack to check if transmission is finished
#                while(s.st_mtime != self.s.st_mtime):
#                    time.sleep(self.delay)
#
#                self.callback()
#                time.sleep(1)
#
#            time.sleep(0.2)


from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import time
import os

class FileChangeNotifier():
    """
    Glue class between cam.py and watchdog. It will launch callback function
    only when ALL files have been update once
    usage:
        files: files to be watched (must be in the same folder!)
        callback: function to call when files are modified
        delay: delay before function call
        enabled: toggles if the callback is called or not
    """

    def __init__(self, files, callback, delay=0.5, enabled=False):

        self.callback = callback
        self.delay = delay
        self.monitoredfiles = files

        self.observer = Observer()
        self.event_handler = FileChangeNotifier.MyHandler(
            callback,
            self.monitoredfiles,
            enabled,
            delay)
        self.setEnabled(enabled)

        self.observer.schedule(
            self.event_handler, os.path.dirname(self.monitoredfiles[0]))
        print "Observing the folder: "+os.path.dirname(self.monitoredfiles[0])
        self.observer.start()

    def setEnabled(self, enabled):
        """
        Enabled =
        True: callback is called
        False: callback is never called (but observer is running)
        """
        self.event_handler.enabled = enabled

    def Stop(self):
        """"
        Stops the watchdog. It can't be restarted after, another one must
        be created
        """
        self.observer.stop()

    class MyHandler(FileSystemEventHandler):

        def __init__(self, callback, files, enabled=False, delay=0.):
            self.callback = callback
            self.enabled = enabled
            self.delay = delay
            self.monitoredfiles = files
            self.hasBeenUpdated = [False]*len(files)
#            print files

        def on_moved(self, event):
#            print event
            try:
                i = self.monitoredfiles.index(event.dest_path)
                self.hasBeenUpdated[i] = True
#                print str(i)+" had changed"
                if all(self.hasBeenUpdated):
                    self.hasBeenUpdated = [False]*len(self.hasBeenUpdated)
                    if self.enabled:
                        time.sleep(self.delay)
                        self.callback()
#                        print "gotcha"

            except ValueError:
                pass


#%% FOR TESTING
if __name__ == '__main__':
    def dosmt():
        print "me"

    d = FileChangeNotifier(os.getcwd(), ['test_0.sis', 'test_1.sis'], dosmt)
    d.setEnabled(True)
    time.sleep(2)

    with open(os.path.join(os.getcwd(), 'test_0.sis'), 'w+b') as fid:
        fid.write(b'adhksaljdh')
    print "1"

    time.sleep(1)

    with open(os.path.join(os.getcwd(), 'test_1.sis'), 'w+b') as fid:
        fid.write(b'adhksaljdh')
