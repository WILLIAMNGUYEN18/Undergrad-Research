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
    #print(str(img))
    plypoints = []
    #traverse array until you find a border pixel, i.e. pixel where any one of the 8 neighbors has an alpha value of 0
    #repeat until you get back to where you started:
        #find neighbor that is also a border pixel
        #append neighbor as array [r,c] to array plypoints
        #handle edge cases / dead ends, but there should be no holes and only one path around the subject

    # Find first pixel with alpha value 1; our starting point
    # Start from top left, going down every coloumn
    first_border = None
    for r in range(0, img.shape[0]):
        for c in range(0, img.shape[1]):

            if img[r,c,3] != 0:
                first_border = (r,c)
                break

        # Exit after finding the first border pixel
        if first_border is not None:
            break
    print("first border pix: " + str(first_border))
    # Case where there is no pixel with alpha value 1
    if first_border is None:
        print("No pixel with alpha value 1!")
        exit()
  
    # Used to loop through neighbors, counterclockwise
    #row,column, (Y,X)
    #east, southeast, south, southwest, west, northwest, north, northeast
    bd = [[0,1], [-1,1], [-1,0], [-1,-1], [0, -1], [1,-1], [1,0], [1,1]]


    point_map = {}
    seen_points = set()
    seen_points.add(first_border)
    point_list = []
    point_list.append(first_border)
    pop_list = []
    seen_list = []
    streak_pop = []

    curr = first_border

    nonStreakPop = True

    #Loop until we find first point again.
    while True:
        children = find_border_children(img, curr, bd)

        point_map[curr] = children

        next = None
        #print("looping through children:" + str(children))
        #it loops through all children not in seen_points, and then iterates through it, so it might be going in reverse order
        for i, child in enumerate(children):
            if child not in seen_points:
                #print(str(i) + ": " + str(child))
                # Found next point
                next = child
        
        
        # Case where there are no children; we go backwards
        if next is None:

            # Check if we are in range of starting pixel. If we are, then we are done
            if points_in_range(first_border, curr, bd):
                # This line appends the start point again, which i assume should not happen
                # point_list.append(first_border)
                break

            # Edge case where first_border is the only pixel with alpha = 1 in range
            if len(point_list) == 1:
                print("first_border is the only pixel with alpha = 1 in range")
                exit()

            # go back and get rid of last element
            curr = point_list.pop()
            if(nonStreakPop):
                streak_pop.append(curr)
                nonStreakPop = False
            pop_list.append(curr)

        else:
            # Found next point
            # Append and move on to next set of children
            point_list.append(next)
            seen_points.add(next)
            seen_list.append(next)
            curr = next
            if(nonStreakPop == False):
                nonStreakPop = True


    print( "border pixel found: " + str(curr))
    #print(seen_list)
    print(len(pop_list))
    print(len(point_list))
    print(len(streak_pop))
    print(point_map)
    num = 1
    for p in seen_list:
        img = cv2.circle(img, (p[1],p[0]), 20, (0,255,255), 1)
        num = num + 1
    #cv2.imshow('image', img)
    fileName = "C:\\Users\\Brilliance\\Desktop\\Projects\\Undergrad-Research\\findborder\\result_circles.png"
    cv2.imwrite(fileName,img)

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

if __name__ == "__main__":
    cudaPlotSilhouette("C:\\Users\\Brilliance\\Desktop\\Projects\\Undergrad-Research\\findborder\\mask0002.png")
