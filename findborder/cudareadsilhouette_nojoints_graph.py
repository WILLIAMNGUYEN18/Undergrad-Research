import cv2
import numpy as np
import math
import sys
import struct

#takes a png file with an alpha channel and plots it to a ply file scaled to the subject's known height.
def cudaPlotSilhouette(png, flength, hsense):

    #4 channel RGBA image represented as a numpy array
    #img[r,c] = [b, g, r, a], img[r,c,3] = alpha channel
    img = cv2.imread(png, -1)
    plypoints = []

    # Used to loop through neighbors, counterclockwise
    #row,column, (Y,X)
    #east, south, west, north
    bd = [[0,1],  [-1,0], [0, -1], [1,0], [-1,1], [-1,-1], [1, -1], [1,1]]
    point_map = {}
    first = True
    
    print("Generating Graph")
    #iterate through the entire map to create a full graph. We can use find_border_children and add them to the map
    #this will allow us to find a comprehensive graph/mapping that we can iterate through
    #let's still store first point as a starting point for the loop
    for r in range(0, img.shape[0]):
        for c in range(0, img.shape[1]):
            currP = (r,c)
            #if alpha value 1 (border)
            #find all borders around this one and map it
            if img[r,c,3] != 0 and is_border_pixel(img,currP, bd):
                
                if(first):
                    firstPoint = currP
                    first = False
                children = find_border_children(img, currP, bd)
                point_map[currP] = children
    if firstPoint is None:
        print("No pixel values with value 1")
        exit()

    seen_points = set()
    seen_points.add(firstPoint)
    currPoint = firstPoint
    point_list = []


    print("Exploring graph")
    #Calvin implementation with direction reversed
    #his was supposed to be clockwise but it went reverse
    while True:
        next = None
        currChildren = point_map[currPoint]
        for i, child in enumerate(currChildren):
            if child not in seen_points:
                next = child
                break
        if next is None:
            if points_in_range(firstPoint, currPoint, bd):
                break
            if len(point_list) == 1:
                print("first_border is the only pixel with alpha = 1 in range")
                exit()
            currPoint = point_list.pop()   
        else:
            point_list.append(next)
            seen_points.add(next)
            currPoint = next

    print( "border pixel found: " + str(currPoint))
    print("point list: " + str(point_list))
    #cv2.imshow('image', img)
    fileName = "C:\\Users\\Brilliance\\Desktop\\Projects\\Undergrad-Research\\findborder\\result.png"
    cv2.imwrite(fileName,img)
    #pixel size in meters
    pixelsize = hsense / float(img.shape[1])
    centerx = img.shape[1] / 2.0
    centery = img.shape[0] / 2.0

    print(hsense)
    print(img.shape[1])

    # Compute "normal" to each point, store in an array. Each point has two neighbors, normal for point j is the average of the two lines j - i and j - k for consecutive points i, j, k
    size = len(point_list)
    print(size)
    for j in range(0, size):

        # Computing indices in a circular fashion
        i = j-1
        k = j+1
        if j == 0:
            i = size-1
        elif j == size - 1:
            k = 0
        #print str(j) + " , " + str(k)
        # Compute normal
        v_i = np.array([point_list[i][1], point_list[i][0], 0])
        v_j = np.array([point_list[j][1], point_list[j][0], 0])
        v_k = np.array([point_list[k][1], point_list[k][0], 0])

        ji_norm = (v_j - v_i) / np.linalg.norm(v_j - v_i)
        jk_norm = (v_j - v_k) / np.linalg.norm(v_j - v_k)

        avg_norm = (ji_norm + jk_norm)
        length = np.linalg.norm(avg_norm)
        if length != 0.0:
            avg_norm /= np.linalg.norm(avg_norm)
        else:	#straight line, so just make it perpindicular
            avg_norm = v_j - v_i
            a0 = avg_norm[0]
            avg_norm[0] = -avg_norm[1]
            avg_norm[1] = a0
            avg_norm = avg_norm / np.linalg.norm(avg_norm)
        #need to do an inside outside test
        inoutpix = np.array(point_list[j]) + 3 * np.array([avg_norm[1], avg_norm[0]])
        #flip normal if its pointing inside
        inoutpixtrunc = np.rint(inoutpix).astype(np.int32)
        if img[inoutpixtrunc[0], inoutpixtrunc[1], 3] != 0:
            avg_norm*= -1

        #print(avg_norm)
        y = avg_norm[1] * -1.0
        x = avg_norm[0]

        normal = (x,y)
        #print(normal)

        #point_list is an array of tuples (r,c) so the y coordinate is first
        #x, y, z, nx, ny, nz. coordinates are metric in m
        plypoints.append([(point_list[j][1] - centerx) * pixelsize, -1 * (point_list[j][0] - centery) * pixelsize, -flength ,normal[0], normal[1], 0.0])

    # #append joint positions
    # dct = open(sys.argv[7] + '/deepcut_imgcoords', 'wb')
    # for j in range(0, joints.shape[1]):
    #     #1/27/2020 dump the joints to a text file so the C++ routine can read it
    #     dct.write(struct.pack('i', joints[0][j]))
    #     dct.write(struct.pack('i', joints[1][j]))
    #     plypoints.append(np.array([(joints[0][j] - centerx) * pixelsize, -1 * (joints[1][j] - centery) * pixelsize, -flength, 0.0, 0.0, 0.0]))
        
    plypoints.append([0.0,0.0,0.0, 0.0, 0.0, 0.0])

    #WRITE TO PLY FILE IN ASCII, don't worry about the shit below
    pointsheader = "ply\nformat ascii 1.0\nelement vertex " + str(len(plypoints)) +"\nproperty float x\nproperty float y\nproperty float z\nproperty float nx\nproperty float ny\nproperty float nz\n"
    edgeheader = "element edge " + str(len(plypoints)-1) + "\nproperty int vertex1\nproperty int vertex2\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nend_header\n"
    plyheader = pointsheader + "end_header\n"  #edgeheader
    plyfile = open("C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\silhouette_4direction.ply", 'w')
    plyfile.write(plyheader)
    for v in plypoints:
        line = ""
        for f in v:
            line += str(f) + " "
        line = line[:-1]
        line = line + '\n'
        plyfile.write(line)
        #write edge v1 v2 r g b
        #http://paulbourke.net/dataformats/ply/
    return plypoints


# Find neighboring points that are also border pixels
def find_border_children(img, point, bd):
    border_pixels = []
    # Append them in the order starting from above, going counterclockwise.
    for i in range(0,8):
        r = point[0] + bd[i][0]
        c = point[1] + bd[i][1]
        child = (r,c)

        if img[r,c,3] != 0 and is_border_pixel(img, child, bd):
            border_pixels.append(child)

    return border_pixels


def is_border_pixel(img, point, bd):

    # Check if there is a pixel surrounding current point with alpha value = 0
    for i in range(0,8):
        r = point[0] + bd[i][0]
        c = point[1] + bd[i][1]

        # Alpha = 0 so this pixel is a border pixel
        if img[r,c,3] == 0:
            return True

    # No neighbors with alpha = 0
    return False


# Check whether p2 is 1 away from p1
def points_in_range(p1, p2, bd):
  
    for i in range(0,8):
        r = p2[0] + bd[i][0]
        c = p2[1] + bd[i][1]
        if r == p1[0] and c == p1[1]:
            return True

    return False


##takes a png file with an alpha channel and plots it to a ply file scaled to the subject's known height.
#def cudaPlotSilhouette(png, flength, hsense):
if __name__ == "__main__":
    cudaPlotSilhouette("C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\mask0002.png", 0.018, 0.0225)
