import cPickle as pickle
import numpy as np
import math
import struct
import sys

def writePCAtoPly(plyfile, points, header, facelist):
	plyfile.write(header)
	for i in range(0, len(points) / 3):
		line = str(points[i*3]) + ' ' + str(points[i*3 + 1]) + ' ' + str(points[i*3 + 2]) + '\n'
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
		if "float nx" not in line and "float ny" not in line and "float nz" not in line:
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
		facelist += face + '\n'
		face = infile.readline()
	return vertices, headerstring, facelist

if __name__ == "__main__":
    for i in range (0, int(sys.argv[3])):
        meanshape = open('../pcainputs/markus_UCSF_pcadata/mean_female_11_18_17.ply', 'r')
        w0 , headerstring, facelist = readTemplatePlyFile(meanshape)
        mpcafile = open("../pcainputs/markus_UCSF_pcadata/PCA_female_11_18_17.pkl")
        mpcapkl = pickle.load(mpcafile)
        mpca = np.transpose(mpcapkl.components_)
        featurestopca = np.load("../pcainputs/markus_UCSF_pcadata/femalef_to_p.npy")
        ftop_pfat = np.load("../pcainputs/markus_UCSF_pcadata/femalef_to_p_pfat.npy")
        
        heightm = np.float64(sys.argv[1])
        weightkgcbrt = np.float64(sys.argv[2]) ** (1. / 3.0)
        f = np.array([[weightkgcbrt], [heightm], [1]])
        test = featurestopca.dot(f)
        s = np.random.normal(np.zeros(test.shape), np.sqrt(np.sqrt(mpcapkl.explained_variance_).reshape(test.shape)))
        #test += s * 2
        result = np.dot(mpca, test) + w0.reshape((w0.shape[0], 1))
        result = result.flatten()
            
        pcatofeatures = np.load("../pcainputs/markus_UCSF_pcadata/femalep_to_f.npy")
        test_append1 = np.append(test.flatten(), 1)
        test_append1 = test_append1.flatten().reshape((test_append1.shape[0] , 1))
        print pcatofeatures.shape
        predfeatures = pcatofeatures.dot(test_append1)
        
        fmass = pow(predfeatures[2][0],3)
        lmass = pow(predfeatures[3][0],3)
        print "weight = " + str(pow(predfeatures[0][0], 3)) + "kg"
        print "height = " + str(predfeatures[1][0]) + "m"
        print "fat mass = " + str(fmass) + "kg"
        print "lean mass = " + str(lmass) + "kg"
        #print "percent fat = " + str(predfeatures[4][0]) + "%"
        print "Visceral fat = " + str(pow(predfeatures[4][0],3)) + "kg"
        print "armlean = " + str(pow(predfeatures[5][0],3)) + "kg"
        print "leg lean = " + str(pow(predfeatures[6][0],3)) + "kg"
        
        print "pfat: " + str(fmass / (fmass + lmass))
        #cube root weight
        #height m
        #fat mass
        #lean mass
        #pfat
        # visceral fat
        # arm lean
        # leg lean
            
        bodyresult = open('generated_' + str(i) + '.ply', 'w')
        writePCAtoPly(bodyresult, result, headerstring, facelist)
        print "pca weights:"
        print test
        
        if len(sys.argv) > 4:
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
            
