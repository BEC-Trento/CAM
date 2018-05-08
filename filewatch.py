#!/usr/bin/python
#-*- coding: latin-1 -*-
"""Watch if a file is changed. Used for automatic reloading."""

import threading
import os
import time
import sys

def dosomething():
    print "filechange!"
    sys.stdout.flush()

class FileChangeNotifier(threading.Thread):

    def __init__(self, filename, callback = dosomething, delay = 2):
        threading.Thread.__init__(self)
        self.filename = filename
        self.callback = callback
        self.delay = delay

        self.s = os.stat(self.filename)
        self.keeprunning = True

    def run(self):
        while self.keeprunning:
            s = os.stat(self.filename)
            if s.st_mtime > self.s.st_mtime:
                self.s = os.stat(self.filename)
                time.sleep(self.delay)
                # horrible hack to check if transmission is finished
                while(s.st_mtime != self.s.st_mtime):
                    time.sleep(self.delay)
                    
                self.callback()
                time.sleep(1)
              
                
            time.sleep(0.2)

    
if __name__ == '__main__':
    d = FileChangeNotifier('test')
    d.start()
    
    
        
