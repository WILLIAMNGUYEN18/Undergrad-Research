import sys
import os
import subprocess
import numpy as np
from subprocess import Popen
import glob

input_dir = sys.argv[1]

jpgs_lower = glob.glob(input_dir + "/*.jpg")
jpgs = glob.glob(input_dir + "/*.JPG")
tifs = glob.glob(input_dir + "/*.tif")

imglist = sorted(jpgs_lower + jpgs + tifs)

dxacsvfile = open(sys.argv[2], 'r')

pcomps = int(sys.argv[3])

lines = dxacsvfile.readlines()
labels = lines.pop(0).split(',')

rows = len(lines)
cols = len(labels)
dxacsv = []
subj_ids = []
#populate arrays with dxa information and subj id for easy indexing
for r in range(0, rows):
    dxacsv.append(lines[r].split(','))    
    subj_ids.append(dxacsv[r][1])

#read in camera params, should be in a txt file.
#focal_length_m, sensorw_m, sensorh_m

cameraparams = open(input_dir + "/camera.txt")
cameralines = cameraparams.readlines()
flength = float(cameralines[0])
hfov = float(cameralines[1])
vfov = float(cameralines[2])

processes = []

for img in imglist:
    #find index of subject in spreadsheet
    subject = img[-20:-11]
    print img
    print subject
    idx = subj_ids.index(subject)
    
    subj_wt = float(dxacsv[idx][9])
    subj_height = float(dxacsv[idx][10]) / 100.0
    gender = 0
    genderprefix = 'males/m'
    if dxacsv[idx][8] != 'M':
        gender = 1
        genderprefix = 'females/f'

    resultdir = "results/OBESITYWEEK2018_RESULTS/" + genderprefix + str(subject) + "_results/d=" + str(pcomps)
    stdoutfile=None
    stderrfile=None
    try:
        os.mkdir(resultdir)
    except OSError as e:
        print(e)

    try:
        stdoutfile = open(resultdir + "/stdout.txt", 'w')
        stderrfile = open(resultdir + "/stderr.txt", 'w')
    except IOError as e:
        print(e)

    #run an instance of the selfie pipeline 

    pargs = " ".join(["./runselfie.sh", str(img), str(gender), str(subj_height), str(subj_wt), str(flength), str(hfov), str(vfov), str(pcomps)])
    print pargs
    processes.append(subprocess.Popen(pargs, stdout=stdoutfile, stderr=stderrfile, shell=True))

for p in processes:
    p.wait()
