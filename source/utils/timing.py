from __future__ import division
#----------------------------------------#
import time
#----------------------------------------#
class Timer:
    def __init__(self):
        self.start_time = 0
        self.stop_time = 0
        self.duration = 0

    def __enter__(self):
        self.start_time = time.time()
        self.duration = 0

    def __exit__(self, *kwargs):
        self.stop_time = time.time()
        self.duration = self.stop_time - self.start_time
        print self.__str__()
        
    def __str__(self):
        return 'Timer duration: {} s.'.format(self.duration)
