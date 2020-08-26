#where the coordinate transform magic happens
#2 * arctan(0.015 / (2 * 0.018)) in degrees
#need to know focal length or sensor is ambiguous


#[x][y][z][1]


import math
import numpy as np
def worldToPixelCoords(T, fov, h, w, f, worldxyz1): 
    #
    TW = np.matmul(T,worldxyz1)
    pixelSize = 2 * f * (math.tan(fov/2))
    x = (f/TW[2]) * TW[0]
    y = (f/TW[2]) * TW(1)

    x = x / pixelSize
    y = y / pixelSize

    x += w/2
    y += -h/2
    return
#worldxyz1 = [x, y, z, 1.0] in world coords

#1. multiply worldxyz1 by T. T is the camera transform. It puts image in the correct viewing position of the camera (which is always at the origin looking down -Z)
#2. project result from (1) onto focal plane at distance f to get film coordinates  x = Xf / Z, y = Yf / Z where X Y Z are from (1). See lecture slides about similar triangle derivation
#3. divide by pixel size. Pixel size is how wide one pixel is on the focal plane (basically the size of a pixel in world coordinates). since you know f you can figure this out from the vertical FOV and the image size h.
#4. multiply y by -1 and add w/2 to x and h/2 to y.
# 
# 
# re-adjust from center to top left 0,0
#
#rather than do it in all matrix like the slides,
