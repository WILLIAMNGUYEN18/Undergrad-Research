import cv2
import numpy as np
import math
import sys
import struct

#takes a png file with an alpha channel and plots it to a ply file scaled to the subject's known height.
def cudaPlotSilhouette(png1, png2):

    #4 channel RGBA image represented as a numpy array
    #img[r,c] = [b, g, r, a], img[r,c,3] = alpha channel
    img1 = cv2.imread(png1, -1)
    img2 = cv2.imread(png2, -1)

    plypoints = []

    # Used to loop through neighbors, counterclockwise
    #row,column, (Y,X)
    #east, southeast, south, southwest, west, northwest, north, northeast
    bd = [[0,1], [-1,1], [-1,0], [-1,-1], [0, -1], [1,-1], [1,0], [1,1]]
    png1_list = []
    png2_list = []
    
    print("Generating Lists")
    #iterate through the entire map to create a full graph. We can use find_border_children and add them to the map
    #this will allow us to find a comprehensive graph/mapping that we can iterate through
    #let's still store first point as a starting point for the loop
    for r in range(0, img1.shape[0]):
        for c in range(0, img1.shape[1]):
            currP = (r,c)
            #if alpha value 1 (border)
            #find all borders around this one and map it
            if img1[r,c,3] != 0 and is_border_pixel(img1,currP, bd):
                png1_list.append(currP)

    for r in range(0, img2.shape[0]):
        for c in range(0, img2.shape[1]):
            currP = (r,c)
            #if alpha value 1 (border)
            #find all borders around this one and map it
            if img2[r,c,3] != 0 and is_border_pixel(img2,currP, bd):
                png2_list.append(currP)

    print("Comparing:")
    print(Diff(png1_list,png2_list))
    #list(set(li1)-set(li2))
    #list(set(li2)-set(li1))
    print("Exclusive to mask 1")
    print(list(set(png1_list) - set(png2_list)))
    print("Exclusive to mask 2")
    print(list(set(png2_list) - set(png1_list)))

    fileName = "C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\result_graph.png"
    cv2.imwrite(fileName,img1)

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

def Diff(li1, li2): 
    return (list(list(set(li1)-set(li2)) + list(set(li2)-set(li1))))     

#C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\mask0002.png
#C:\\Users\\Brilliance\\Desktop\\Projects\\Undergrad-Research\\findborder\\mask0002.png
if __name__ == "__main__":
    cudaPlotSilhouette("C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\mask0001.png","C:\\Users\\willi\\Desktop\\Programming\\Undergrad-Research\\findborder\\mask0002.png")
