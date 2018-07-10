#!/usr/bin/env python

class utils:
    def __init__(self):
        pass

    def dynamicProfileBins(self,hits,coord='x',relError=0.1):
        minPixels = max(1,1/relError/relError)
        index = 0 if coord=='x' else 1
        xmin=min([h[index] for h in hits])
        xmax=max([h[index] for h in hits])
        x=int(xmin)
        xedges=[x]
        integral=0
        while x<xmax:
            if integral<minPixels:
                integral += sum([int(h[index])==int(x) for h in hits])
            else:
                xedges.append(x)
                integral=0
            x+=1
        xedges.append(int(xmax))
        return xedges

    def rotate_around_point(self, hit, dir, pivot, inverse=False):
        x,y = hit[:-1]
        ox, oy = pivot
        cos,sin = dir
        if inverse: cos = -1*cos
        qx = ox + cos * (x - ox) + sin * (y - oy)
        qy = oy - sin * (x - ox) + cos * (y - oy)
        return qx, qy

