import sys
import os
import subprocess
import numpy as np
from subprocess import Popen
import glob
import fnmatch
from os import path


input_dir = sys.argv[1]
print("input dir: " + input_dir)
plys = sorted(glob.glob(input_dir + "/*.ply"))
#print("ply files: " + plys)

#O:\unix\projects\grail\iytian\separatepcperdevice_102920\dxa\SUQ4_DXA_200421.csv
#location will be in separatepcperdevice_102920 so we can do relative path
#dxacsvfile = open(".\\dxa\\SUQ4_DXA_200421.csv", 'r')
#dxacsvfile = open("O:\unix\projects\grail\iytian\separatepcperdevice_102920\dxa\SUQ4_DXA_200421.csv", 'r')
dxacsvfile = open("./dxa/SUQ4_DXA_200421.csv", 'r')


lines = dxacsvfile.readlines()
labels = lines.pop(0).split(',')
#print("Labels: " + labels)
rows = len(lines)
cols = len(labels)

sid_idx = labels.index("SUBJECTID")

subject_ids = []
dxacsv = []

for r in range(0, rows):
    dxacsv.append(lines[r].split(','))
    subject_ids.append(dxacsv[r][sid_idx])

#initialize all the arrays outside of the loop
percentFat = []
leanMass = []
fatMass = []
viscFat = []
legLean = []
armLean = []
trunkLean = []
trunkFat = []
legFat = []
armFat = []
FMI = []
FFMI = []

for mesh in plys:
    print("mesh")
    print(mesh)
    #parse for subject id with -1 being last token(filename) [0:-4 drops the .ply]
    #need '\\' for windows and '/' for linux
    subject = mesh.split('\\')[-1][0:-4]
    print("subject")
    print(subject)

    sid = subject.split('_')[0]
    print("sid: " + sid)

    #find index of subject in spreadsheet
    idx = subject_ids.index(sid)
    print("index in dxacsv: " + str(idx))    
    subj_fatMass = float(dxacsv[idx][labels.index("")])