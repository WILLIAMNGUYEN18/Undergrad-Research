import sys
import os
import subprocess
import numpy as np
from subprocess import Popen
import glob
import fnmatch
from os import path

overwrite = True

print "OVERWRITE FLAG: " + str(overwrite)

input_dir = sys.argv[1]
pcomps = int(sys.argv[2])
baseresdir = sys.argv[3]
devicei = int(sys.argv[4])
print baseresdir
split = input_dir.split('/')[-2]
print split

plys = sorted(glob.glob(input_dir + "/*.ply"))

dxacsvfile = open("./dxa/SUQ4_DXA_200421.csv", 'r')

device = 'fit3d'
if devicei == 1:
    device = 'styku'
if devicei == 2:
    device = 'ss'
if devicei == 3:
    device = 'super'
if devicei == 4:
    device = 'tpose'
if devicei == 5:
    device = 'naked'
lines = dxacsvfile.readlines()
labels = lines.pop(0).split(',')

rows = len(lines)
cols = len(labels)
dxacsv = []
subj_ids = []
sid_idx = labels.index("SUBJECTID")
#populate arrays with dxa information and subj id for easy indexing
for r in range(0, rows):
    dxacsv.append(lines[r].split(','))    
    subj_ids.append(dxacsv[r][sid_idx])
#print subj_ids

#read in camera params, should be in a txt file.
#focal_length_m, sensorw_m, sensorh_m


processes = []

for mesh in plys:
    #find index of subject in spreadsheet
    subject = mesh.split('/')[-1][0:-4]
    #print mesh
    #print subject
    sid = subject.split('_')[0]
    idx = subj_ids.index(sid)
    gender = dxacsv[idx][labels.index("SEX")]
    if gender == 'M':
        gender = 0
    else:
        gender = 1
    
    try:
        os.mkdir(baseresdir + "/males/" + device)
        os.mkdir(baseresdir +  "/females/" + device)
    except OSError as e:
        print(e)

    genderprefix = 'males'
    if gender == 1:
        genderprefix = 'females'    
        
    subj_wt = float(dxacsv[idx][labels.index("WBTOT_MASS")]) / 1000.0
    subj_height = float(dxacsv[idx][labels.index("HEIGHT")]) / 100.0
    
    subjdir = baseresdir + "/" + genderprefix + '/'  + device + "/" + split + "/" + genderprefix[0] + str(subject)
    print subjdir
    resultdir = subjdir + "/d=" + str(pcomps)
    stdoutfile=None
    stderrfile=None
    try:
        os.mkdir(subjdir)
        os.mkdir(resultdir)
    except OSError as e:
        print(e)

    if not overwrite and path.exists(resultdir + "/result.ply") and os.stat(resultdir + "/result.ply").st_size != 0: #if empty file, overwrite
        print "result exists for " + subject
        continue

    try:
        stdoutfile = open(resultdir + "/stdout.txt", 'w')
        stderrfile = open(resultdir + "/stderr.txt", 'w')
    except IOError as e:
        print(e)
        
    print subject
    #run an instance of the selfie pipeline 
    os.chdir("./C")
    #lda = 0.001 #original
    lda = 0.01
    #lda = 0.1   
    #lda = 0.0 #for training set members only
    #tol = 0.01
    tol = 0.1
    #tol = 0.0 #training only
    pargs = " ".join(["./pcamatch", "../" + mesh, str(lda) ,str(gender), str(subj_height), str(subj_wt) , "../" + resultdir, str(tol), device])
    print pargs
    processes.append(subprocess.Popen(pargs, stdout=stdoutfile, stderr=stderrfile, shell=True))
    os.chdir("../")
    while len(processes) >= 16:
        processes[0].wait()
        processes.pop(0)

for p in processes:
    p.wait()
