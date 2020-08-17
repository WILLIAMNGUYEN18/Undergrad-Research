import numpy as np
import math
#where the coordinate transform magic happens
#2 * arctan(0.015 / (2 * 0.018)) in degrees
#need to know focal length or sensor is ambiguous
def worldToPixelCoords(T, fov, h, w, f, worldxyz1):
    #matrix multiply worldxyz1 by T where T is the camera transform
    cameraxyz1 = np.matmul(worldxyz1,T)

    #project result from (1) onto focal plane at distance f to get film coordinates
    #x = Xf / Z, y = Yf / Z where X Y Z are from (1)

    #actually don't think we need these
    #X = x * (z/f)
    #Y = y * (z/f)

    #Are we only generating this for a single point?
    #Or are we generating this for an entire matrix?

    #assuming worldxyz1 has 4 indexes
    x1 = worldxyz1[0] * f / worldxyz1[2]
    y1 = worldxyz1[1] * f / worldxyz1[2]


    #find pixel size from right triangle of f and h (a^2 + b^2 = c^2)
    pxAngle = fov / 2
    pxSize = math.sqrt(pow(f,2) + pow(h,2))






    

    





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
