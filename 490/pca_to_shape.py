import cPickle as pickle
import numpy as np
import math
import struct
import sys
import os
sys.path.insert(0, '../')
import surfacenormals

def writePCAtoPly(plyfile, points, header, facelist):
	plyfile.write(header)
        normals = surfacenormals.recomputeSurfaceNormals(points, facelist)
	for i in range(0, len(points) / 3):
		line = str(points[i*3]) + ' ' + str(points[i*3 + 1]) + ' ' + str(points[i*3 + 2]) + ' ' + str(normals[i][0]) + ' ' + str(normals[i][1]) + ' ' + str(normals[i][2]) + '\n'
		plyfile.write(line)
	plyfile.write(facelist)

def readTemplatePlyFile(infile):
	headerEnd = False
	vertNum = 0
	headerstring = ""
	facelist = ""
	while not headerEnd:
		line = infile.readline()
		#print line
		if line == "end_header\n":
			headerEnd = True
		if "element vertex" in line:
			vertNum = int(line.split()[2])
			#print vertNum
		#if "float nx" not in line and "float ny" not in line and "float nz" not in line:
		headerstring += line

	vertices = np.zeros(vertNum * 3)
	for i in range(0, vertNum):
		coord = infile.readline()
		values = coord.split()
		point = [float(values[0]), float(values[1]), float(values[2])]
		vertices[i*3] = point[0]
		vertices[i*3 + 1] = point[1]
		vertices[i*3 + 2] = point[2]
	
	face = infile.readline()
	while len(face) > 0:
		facelist += face
		face = infile.readline()
	return vertices, headerstring, facelist

if __name__ == "__main__":
  
    gender = sys.argv[4]
    for i in range (0, int(sys.argv[3])):
        #meanshape = open('../pcainputs/markus_UCSF_pcadata/mean_male_11_18_17.ply', 'r')
        
        meanshape = open('../pcainputs/nov2018_bennett_8020/ply/male_avg_8020split_ascii.ply', 'r')
        pcatofeatures = np.load("../pcainputs/nov2018_bennett_8020/malep_to_f_50.npy")
        featurestopca = np.load("../pcainputs/nov2018_bennett_8020/malef_to_p_50.npy")
        mpca = np.load("../pcainputs/nov2018_bennett_8020/mpca8020_50_np.npy")
	sigmas = np.load("../pcainputs/nov2018_bennett_8020/msigma_50_np.npy")
	if gender == "1":
	    meanshape = open('../pcainputs/nov2018_bennett_8020/ply/female_avg_8020split_ascii.ply', 'r')
            pcatofeatures = np.load("../pcainputs/nov2018_bennett_8020/femalep_to_f_50.npy")
            featurestopca = np.load("../pcainputs/nov2018_bennett_8020/femalef_to_p_50.npy")
            mpca = np.load("../pcainputs/nov2018_bennett_8020/fpca8020_50_np.npy")
            sigmas = np.load("../pcainputs/nov2018_bennett_8020/fsigma_50_np.npy")

	w0 , headerstring, facelist = readTemplatePlyFile(meanshape)
        #mpcapkl = pickle.load(mpcafile)
        #mpca = np.transpose(mpcapkl.components_)
        #ftop_pfat = np.load("../pcainputs/nov2018_bennett_8020/malef_to_p_pfat.npy")
        
        heightm = np.float64(sys.argv[1])
        weightkgcbrt = np.float64(sys.argv[2]) ** (1. / 3.0)
        f = np.array([[weightkgcbrt], [heightm], [1]])
        test = featurestopca.dot(f)
        s = np.random.normal(np.zeros(test.shape), sigmas.reshape(test.shape))
        test += s * 0.5	#perturbs by 3 standard deviations
        result = np.dot(mpca, test) + w0.reshape((w0.shape[0], 1))
        result = result.flatten()
            
        test_append1 = np.append(test.flatten(), 1)
        test_append1 = test_append1.flatten().reshape((test_append1.shape[0] , 1))
        print pcatofeatures.shape
        predfeatures = pcatofeatures.dot(test_append1)
        
        fmass = pow(predfeatures[2][0],3)
        lmass = pow(predfeatures[3][0],3)
	armleanmass = pow(predfeatures[5][0],3)
	legleanmass = pow(predfeatures[6][0],3)

        print "weight = " + str(pow(predfeatures[0][0], 3)) + "kg"
        print "height = " + str(predfeatures[1][0]) + "m"
        print "fat mass = " + str(fmass) + "kg"
        print "lean mass = " + str(lmass) + "kg"
        #print "percent fat = " + str(predfeatures[4][0]) + "%"
        print "Visceral fat = " + str(pow(predfeatures[4][0],3)) + "kg"
        print "armlean = " + str(armleanmass) + "kg"
        print "leg lean = " + str(legleanmass) + "kg"
        
        print "pfat: " + str(fmass / (fmass + lmass))
        print "appLMI: " + str((armleanmass + legleanmass) / pow(heightm, 2))
	#cube root weight
        #height m
        #fat mass
        #lean mass
        #pfat
        # visceral fat
        # arm lean
        # leg lean

	#WILL'S CHANGES

	#string file to be made into text file
	comp_pred = str(predfeatures[1][0])
	comp_pred += "\n" + str(pow(predfeatures[0][0], 3))
	comp_pred += "\n" + str(fmass)
	comp_pred += "\n" + str(lmass)
	comp_pred += "\n" + str(pow(predfeatures[4][0],3))
	comp_pred += "\n" + str(armleanmass)
	comp_pred += "\n" + str(legleanmass)
        comp_pred += "\n" + str(fmass / (fmass + lmass))
	
	print("comp_pred file:")	
	print(comp_pred)
	print("")

	#directory name
	dirnam = "F_"
	if gender == '0':
		dirnam = "M_"
	dirnam += "h" + str('%.3f'%(predfeatures[1][0])) + "_w" + str('%.3f'%(pow(predfeatures[0][0], 3)))
	#Using i from 0 to sysargs[3] (Sample #) here for the current iteration	
	dirnam += "_" + str(i) + "_" + str('%.3f'%(fmass / (fmass + lmass)))
	
	print("directory name:")	
	print(dirnam)
	print("")
	#Rename ply file

		#os.rename('a.txt', 'b.kml')
		#where is the ply file?	
		
	#Create directory
		#need to define a path?
		#path is probably ../599g1_training/male_synth?
		#pathing from current directory of validation? 
	path = "../599g1_training/males_synth"	
	if gender == '1':	
		path = "../599g1_training/females_synth"
	
	path += "/" + dirnam
	
	directory = path
	#need to add info to mkdir?
	#unsure how to name it
	
	if not os.path.exists(directory):
		os.makedirs(directory)
	

	#send comp_pred to directory
	#save path
	save_path = directory
	completeName = os.path.join(save_path, "comp_pred.txt")
	comp_file = open(completeName, "w")
	comp_file.write(comp_pred)
	comp_file.close()	
        
	#change title of bodyresult
	plyname = os.path.join(save_path, dirnam + ".ply")
        bodyresult = open( plyname, 'w')
        writePCAtoPly(bodyresult, result, headerstring, facelist)
        #print "pca weights:"
        #print test

	
	
	#send PLY to directory
		
	
	

        #HEY WILL
        #mkdir output directory with format string
        #make txt file with all stats in it
        #OUTPUT ALL SHIT INTO THIS DIR
	#output comp_pred
	

        
#don't use this
        if len(sys.argv) > 5:
            print "meep"
            pfatin = np.float64(sys.argv[4])
            print pfatin
            fmassincbrt = (pfatin / 100.0 * pow(weightkgcbrt,3)) ** (1. / 3.0)
            f2 = np.array([[weightkgcbrt], [heightm], [fmassincbrt], [1]])
            testpfat = ftop_pfat.dot(f2)
            resultpfat = np.dot(mpca, testpfat) + w0.reshape((w0.shape[0], 1))
            resultpfat = resultpfat.flatten()
            pfatvis = open('generated_' + str(i) + 'pfat.ply', 'w')
            writePCAtoPly(pfatvis, resultpfat, headerstring, facelist)
            
            testpfat_append1 = np.append(testpfat.flatten(), 1)
            testpfat_append1 = testpfat_append1.flatten().reshape((testpfat_append1.shape[0] , 1))
            predfeatures = pcatofeatures.dot(testpfat_append1)
            fmass = pow(predfeatures[2][0],3)
            lmass = pow(predfeatures[3][0],3)
            print "weight = " + str(pow(predfeatures[0][0], 3)) + "kg"
            print "height = " + str(predfeatures[1][0]) + "m"
            print "fat mass = " + str(pow(predfeatures[2][0],3)) + "kg"
            print "lean mass = " + str(pow(predfeatures[3][0],3)) + "kg"
            #print "percent fat = " + str(predfeatures[4][0]) + "%"
            print "Visceral fat = " + str(pow(predfeatures[4][0],3)) + "kg"
            print "armlean = " + str(pow(predfeatures[5][0],3)) + "kg"
            print "leg lean = " + str(pow(predfeatures[6][0],3)) + "kg"
            print "pfat: " + str(fmass / (fmass + lmass))
            print testpfat
            
