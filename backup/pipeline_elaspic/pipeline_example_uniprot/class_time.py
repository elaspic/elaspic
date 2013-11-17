# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 10:25:09 2013

@author: niklas
"""

from time import time


class runTime():
    
    def __init__(self):
        self.start = time()
    
    def __call__(self, mode=''):
        t = time() - self.start
        if mode == '':
            return t
        elif mode == 'human':
            return self.secondsToStr(t)
        else:
            print 'Unregognised format!'
            return
    
    def secondsToStr(self, t):
        """
        credit goes to:
            http://stackoverflow.com/questions/1557571/how-to-get-time-of-a-python-program-execution/1557906#1557906
        """
        return "%d:%02d:%02d.%03d" % \
            reduce(lambda ll,b : divmod(ll[0],b) + ll[1:],
                [(t*1000,),1000,60,60])
    