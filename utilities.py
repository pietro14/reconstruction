#!/usr/bin/env python
import subprocess

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

    def gen_rand_limit(self, x1, x2, y1, y2, maxx=2048, maxy=2048):
        import random
        # generate x, y O(1)
        # --x
        left = random.randrange(0, x1)
        right = random.randrange(x2+1, maxx)
        withinx = random.randrange(x1, x2+1)
        # adjust probability of a point outside the box columns
        # a point outside has probability (1/(maxx-w)) v.s. a point inside has 1/w
        # the same is true for rows. adjupx/y adjust for this probability 
        w = abs(x2-x1)
        h = abs(y2-y1)
        adjpx = ((maxx - w)/w/2)
        x = random.choice([left, right] * adjpx + [withinx])
        # --y
        top = random.randrange(0, y1)
        bottom = random.randrange(y2+1, maxy)
        withiny = random.randrange(y1, y2+1)
        if x == left or x == right:
            adjpy = ((maxy- h)/h/2)
            y = random.choice([top, bottom] * adjpy + [withiny])
        else:
            y = random.choice([top, bottom])
        return x, y 

    def get_git_revision_hash(self):
        return subprocess.check_output(['git', 'rev-parse', 'HEAD'])
