import numpy as np

def computeAvgVector(v_list):
    total = np.array([0,0,0])
    n = len(v_list)
    for i in range(0, n):
        total = np.add(total, v_list[i])

    total *= 1.0/n

    return total / np.linalg.norm(total)

def recomputeSurfaceNormals(points3d, facelist):
    #points3d is flat
    points3d = np.reshape(points3d, (len(points3d) / 3, 3))
    facelist = facelist.split('\n')[:-1] #extra whitespace at end for some reason
    newnormalslist = []
    for i in range(0, len(points3d)):
        newnormalslist.append([])

    for i in range (0, len(facelist)):
        facei = facelist[i].split()[1:] #get rid of extra 3 at the front
        facei = [int(facei[0]), int(facei[1]), int(facei[2])]
	v1 = points3d[facei[0]] - points3d[facei[1]]
        v2 = points3d[facei[2]] - points3d[facei[1]]
        v3 = np.cross(v2, v1)
        v3 /= np.linalg.norm(v3)

        newnormalslist[facei[0]].append(v3)
        newnormalslist[facei[1]].append(v3)
        newnormalslist[facei[2]].append(v3)
    

    newnormals = []
    for i in range(0,  len(newnormalslist)):
        newnormals.append(computeAvgVector(newnormalslist[i]))
    
    
    return newnormals

