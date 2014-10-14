# flake tools for flagr
from __future__ import print_function, division, generators
import time, itertools, math, sys
import numpy as np


def next_neighbour(nn='diag'):
    """
    Creates the next-neighbours vector. First creates all combinations of
    (0,1,-1) and then all permutations, keeping only the unique vector.
    pads a zero for the fourth dimension at the end.
    Returns index-type
    """
    if nn == 'cross':
        return np.array(((1, 0, 0), (0, 1, 0), (0, 0, 1),
                        (-1, 0, 0), (0, 0, -1), (0, -1, 0)))
    elif nn == 'real':
        return np.array(((0,0,-1),(1,0,-1),(0,1,-1),
        (1,0,0),(1,-1,0),(0,-1,0),(-1,0,0),(-1,1,0),(0,1,0),
        (0,0,1),(-1,0,1),(0,-1,1)))
    else:
        types = itertools.combinations_with_replacement([1,-1,0],3)
        perms = []
        for i in types:
            t = set(itertools.permutations(i))
            while t:
                perms.append(t.pop())
        return perms[:-1]
        #return np.pad(perms[:-1], ((0,0),(0,1)), mode='constant')

def vector(i, j, k, regular=True, shift=0):
    """
    Build crystal as i*a+j*b+k*c with lattice vectors:
    a = [2, 0, 0], b = [1, sqrt(3), 0], c = [1, 1/sqrt(3), 2/sqrt(3)]
    Regular refers to FCC-close-packing order A-B-C, where not
    regular refers to inverted order C-B-A, needed for twin-planes.
    returns coordinate-type
    """
    ABCABC = np.array(( 2*i + j + (k+shift)%3,
                        j*np.sqrt(3) + ((k+shift)%3)/np.sqrt(3),
                         k*2./np.sqrt(3), 0))
    CBACBA = np.array(( 2*i + j - (k+1+shift)%3,
                        j*np.sqrt(3) - ((k+1+shift)%3)/np.sqrt(3),
                         k*2./np.sqrt(3), 0))
    if regular:
        return ABCABC
    else:
        return CBACBA

class RedirectStdoutTo:
    """
    Class for redirecting the output to a chosen destination (e.g. a file)
    The exit argument automatically reverts the standard output on exit, see
    _with_ statement.
    """
    def __init__(self, out_new):
        self.out_new = out_new
    def __enter__(self):
        self.out_old = sys.stdout
        sys.stdout = self.out_new
    def __exit__(self, *args):
        sys.stdout = self.out_old

#progress bar:
# http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console

class ProgressBar():
    DEFAULT_BAR_LENGTH = float(65)

    def __init__(self, end, start=0):
        self.end    = end
        self.start  = start
        self._barLength = ProgressBar.DEFAULT_BAR_LENGTH

        self.setLevel(self.start)
        self._plotted = False

    def setLevel(self, level, initial=False):
        self._level = level
        if level < self.start:  self._level = self.start
        if level > self.end:    self._level = self.end

        self._ratio = float(self._level - self.start) / float(self.end - self.start)
        self._levelChars = int(self._ratio * self._barLength)

    def plotProgress(self):
        sys.stdout.write("\r  %3i%% [%s%s]" %(
            int(self._ratio * 100.0),
            '=' * int(self._levelChars),
            ' ' * int(self._barLength - self._levelChars),
        ))
        self._plotted = True

    def setAndPlot(self, level):
        oldChars = self._levelChars
        self.setLevel(level)
        if (not self._plotted) or (oldChars != self._levelChars):
            self.plotProgress()

    def __del__(self):
        sys.stdout.write("\n")

if __name__ == "__main__":
    import time
    count = 5
    print( "starting things:")

    pb = ProgressBar(count)

    curProgress = 0
    #pb.plotProgress()
    while curProgress <= count:
        pb.setAndPlot(curProgress)
        curProgress += 1
        time.sleep(1)
    del pb

    print( "done")
