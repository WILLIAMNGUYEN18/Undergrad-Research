import cv2
import numpy as np
import math
import sys
import struct

#takes a png file with an alpha channel and plots it to a ply file scaled to the subject's known height.
def cudaPlotSilhouette(png):

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


    print(point_map)
    #print(len(point_map.keys()))
    #print(len(point_map.values()))
    print(point_list)
    print("point_map size: " + str(len(point_map)))
    print("point_list size: " + str(len(point_list)))
    print("first pixel: " + str(firstPoint))
    print("last pixel: " + str(currPoint))
    fileName = "C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\result_graph.png"
    cv2.imwrite(fileName,img)

# Find neighboring points that are also border pixels
def find_border_children(img, point, bd):
    #print("Finding borders for point: (" + str(point[0]) + "," + str(point[1]) + ")")
    border_pixels = []
    # Append them in the order starting from above, going counterclockwise.
    for i in range(0,4):
        r = point[0] + bd[i][0]
        c = point[1] + bd[i][1]
        child = (r,c)

        if img[r,c,3] != 0 and is_border_pixel(img, child, bd):
            #print("Adding Border Pixel: (" + str(r) + "," + str(c) + ")" )
            border_pixels.append(child)

    return border_pixels


def is_border_pixel(img, point, bd):

    # bd 0 --> 7, [0] is row, [1] is column.
    # point is current point
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

#C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\mask0002.png
#C:\\Users\\Brilliance\\Desktop\\Projects\\Undergrad-Research\\findborder\\mask0002.png
if __name__ == "__main__":
    cudaPlotSilhouette("C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\mask0002.png")
