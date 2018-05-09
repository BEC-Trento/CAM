#!/usr/bin/python
#-*- coding: latin-1 -*-
"""Watch if a file is changed. Used for automatic reloading."""

import threading
import os
import time
import sys

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


class MyHandler(FileSystemEventHandler):

    def __init__(self, callback, files, enabled=False):
        self.callback = callback
        self.enabled = enabled
        self.monitoredfiles = files
        print files

    def on_modified(self, event):
        if (event.src_path in self.monitoredfiles) and self.enabled:
            print "gotcha"
#            self.callback()
#        print event


class FileChangeNotifier():
    """
    Glue class between cam.py and watchdog
    usage:
        folder: folder where to watch
        files: files to be watched
        callback: function to call when files are modified
        delay: delay before function call
        enabled: toggles if the callback is called or not
    """

    def __init__(self, folder, files, callback, delay=0.5, enabled=False):

        self.callback = callback
        self.delay = delay
        self.monitoredfiles = [os.path.join(folder, f) for f in files]

        self.observer = Observer()
        self.event_handler = MyHandler(None, self.monitoredfiles, enabled)
        self.setEnabled(enabled)

        self.observer.schedule(self.event_handler, folder)
        self.observer.start()

    def setEnabled(self, enabled):
        self.event_handler.enabled = enabled


if __name__ == '__main__':
    d = FileChangeNotifier(os.getcwd(), ['test.sis'], lambda x: None)
    d.setEnabled(True)
