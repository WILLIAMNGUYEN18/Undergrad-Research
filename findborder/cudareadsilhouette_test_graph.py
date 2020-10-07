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
    #east, southeast, south, southwest, west, northwest, north, northeast
    bd = [[0,1], [-1,1], [-1,0], [-1,-1], [0, -1], [1,-1], [1,0], [1,1]]
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
    maxpoint_list = []
    maxlistlength = len(point_list)
    maxpoint = firstPoint
    pop_list = []
    seen_list = []


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
            if(maxlistlength < len(point_list)):
                maxlistlength = len(point_list)
                maxpoint_list = point_list.copy()
                maxpoint = currPoint    
            currPoint = point_list.pop()
            pop_list.append(currPoint)   
        else:
            point_list.append(next)
            seen_points.add(next)
            currPoint = next
            seen_list.append(next)


    #print(point_map)
    #print(len(point_map.keys()))
    #print(len(point_map.values()))
    #print(point_list)
    print("point_map size: " + str(len(point_map)))
    print("point_list size: " + str(len(point_list)))
    print("first pixel: " + str(firstPoint))
    print("last pixel: " + str(currPoint))
    print(maxpoint_list)
    print("maxpoint: " + str(maxpoint))
    print("maxlist length: " + str(maxlistlength))
    fileName = "C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\result_graph.png"
    cv2.imwrite(fileName,img)
    outlinefileName = "C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\resultgraph_outline.png"
    outlineImage = np.zeros(shape=[img.shape[0],img.shape[1], 3], dtype=np.uint8)

    for p in seen_list:
        outlineImage[p[0],p[1],0] = 0
        outlineImage[p[0],p[1],1] = 250
        outlineImage[p[0],p[1],2] = 0

    print("seen_list length: " + str(len(seen_list)))
    print("pop_list length: " + str(len(pop_list)))
    cv2.imwrite(outlinefileName, outlineImage)
    outline_pop_filename = "C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\resultgraph_outlinepop.png"
    for p in pop_list:
        outlineImage[p[0],p[1],0] = 0
        outlineImage[p[0],p[1],1] = 0
        outlineImage[p[0],p[1],2] = 250

    cv2.imwrite(outline_pop_filename, outlineImage)
    outline_maxlist_filename = "C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\resultgraph_outlinemaxlist.png"
    
    n = 0
    for p in maxpoint_list:
        n = n + 1
        variance = (255 / len(maxpoint_list)) * n
        outlineImage[p[0],p[1],0] = variance
        outlineImage[p[0],p[1],1] = 255 - variance
        outlineImage[p[0],p[1],2] = 250

    cv2.imwrite(outline_maxlist_filename, outlineImage)

# Find neighboring points that are also border pixels
def find_border_children(img, point, bd):
    #print("Finding borders for point: (" + str(point[0]) + "," + str(point[1]) + ")")
    border_pixels = []
    # Append them in the order starting from above, going counterclockwise.
    for i in range(0,8):
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
