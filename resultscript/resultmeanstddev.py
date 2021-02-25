import sys
import os
import subprocess
import numpy as np
from subprocess import Popen
import glob
import fnmatch
from os import path
from datetime import date
from scipy import stats

def standarddeviation(value):
    print(str(np.format_float_positional(np.std(np.array(value)), 3)))
def mean(value):
    print(str(np.format_float_positional(np.mean(np.array(value)), 3)))
def max(value):
    print(str(np.format_float_positional(np.max(np.array(value)), 3)))
def min(value):
    print(str(np.format_float_positional(np.min(np.array(value)), 3)))
def pval(train,test):
    print(str(stats.ttest_ind(train, test)[1]))

def scan_age(birthdate, scandate):
    barray = birthdate.split('/')
    sarray = scandate.split('/')
    print(str(barray))
    print(str(sarray))
    #year/month/day
    bday = date(int(barray[2]), int(barray[1]), int(barray[0]))
    sday = date(int(sarray[2]), int(sarray[1]), int(sarray[0]))
    return (sday - bday).year 


dir_train = sys.argv[1]
dir_bootstrap = sys.argv[2]
dir_test = sys.argv[3]
print("train dir: " + dir_train)
print("bootstrap dir: " + dir_bootstrap)
print("test dir: "+ dir_test)
plys = sorted(glob.glob(dir_train + "/*.ply") + glob.glob(dir_bootstrap + "/*.ply"))

#print("ply files: " + plys)

#O:\unix\projects\grail\iytian\separatepcperdevice_102920\dxa\SUQ4_DXA_200421.csv
#location will be in separatepcperdevice_102920 so we can do relative path
#dxacsvfile = open(".\\dxa\\SUQ4_DXA_200421.csv", 'r')
#dxacsvfile = open("O:\unix\projects\grail\iytian\separatepcperdevice_102920\dxa\SUQ4_DXA_200421.csv", 'r')
dxacsvfile = open("./SUQ4_DXA_200421.csv", 'r')


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
age = []
height = []
mass = []
BMI = []

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
    subject = mesh.split('/')[-1][0:-4]
    print("subject")
    print(subject)

    sid = subject.split('_')[0]
    print("sid: " + sid)

    #find index of subject in spreadsheet
    idx = subject_ids.index(sid)
    print("index in dxacsv: " + str(idx))    
    #subj_fatMass = float(dxacsv[idx][labels.index("")])
    #print("Appending to percentFat: " + dxacsv[idx][labels.index("WBTOT_PFAT")])


    height.append(float(dxacsv[idx][labels.index("HEIGHT")]) / 100.0)
    mass.append(float(dxacsv[idx][labels.index("WBTOT_MASS")]) / 1000.0)
    age.append(scan_age(dxacsv[idx][labels.index("BIRTHDATE")], dxacsv[idx][labels.index("SCAN_DATE")]))
    BMI.append(float(dxacsv[idx][labels.index("BMI")]))

    percentFat.append(float(dxacsv[idx][labels.index("WBTOT_PFAT")]))
    leanMass.append(float(dxacsv[idx][labels.index("WBTOT_LEAN")]) / 1000.0)
    fatMass.append(float(dxacsv[idx][labels.index("WBTOT_FAT")]) / 1000.0)
    viscFat.append(float(dxacsv[idx][labels.index("VFAT_MASS")]) / 1000.0)
    legLean.append(float(dxacsv[idx][labels.index("LEG_LEAN")]) / 1000.0)
    armLean.append(float(dxacsv[idx][labels.index("ARM_LEAN")]) / 1000.0)
    trunkLean.append(float(dxacsv[idx][labels.index("TRUNK_LEAN")]) / 1000.0)
    trunkFat.append(float(dxacsv[idx][labels.index("TRUNK_FAT")]) / 1000.0)
    legFat.append(float(dxacsv[idx][labels.index("LEG_FAT")]) / 1000.0)
    armFat.append(float(dxacsv[idx][labels.index("ARM_FAT")]) / 1000.0)
    currHeight = float(dxacsv[idx][labels.index("HEIGHT")]) / 100.0
    fMass = float(dxacsv[idx][labels.index("WBTOT_FAT")]) / 1000.0
    lMass = float(dxacsv[idx][labels.index("WBTOT_LEAN")]) / 1000.0
    FMI.append(fMass / (currHeight * currHeight))
    FFMI.append(lMass / (currHeight * currHeight))


    
plys_test = sorted(glob.glob(dir_test + "/*.ply"))

age_test = []
height_test = []
mass_test = []
BMI_test = []

percentFat_test = []
leanMass_test = []
fatMass_test = []
viscFat_test = []
legLean_test = []
armLean_test = []
trunkLean_test = []
trunkFat_test = []
legFat_test = []
armFat_test = []
FMI_test = []
FFMI_test = []
for mesh in plys_test:
    print("mesh")
    print(mesh)
    subject = mesh.split('/')[-1][0:-4]
    print("subject")
    print(subject)
    sid = subject.split('_')[0]
    print("sid: " + sid)
    idx = subject_ids.index(sid)
    print("index in dxacsv: " + str(idx))
    height_test.append(float(dxacsv[idx][labels.index("HEIGHT")]) / 100.0)
    mass_test.append(float(dxacsv[idx][labels.index("WBTOT_MASS")]) / 1000.0)
    age_test.append(scan_age(dxacsv[idx][labels.index("BIRTHDATE")], dxacsv[idx][labels.index("SCAN_DATE")]))
    BMI_test.append(float(dxacsv[idx][labels.index("BMI")]))

    percentFat_test.append(float(dxacsv[idx][labels.index("WBTOT_PFAT")]))
    leanMass_test.append(float(dxacsv[idx][labels.index("WBTOT_LEAN")]) / 1000.0)
    fatMass_test.append(float(dxacsv[idx][labels.index("WBTOT_FAT")]) / 1000.0)
    viscFat_test.append(float(dxacsv[idx][labels.index("VFAT_MASS")]) / 1000.0)
    legLean_test.append(float(dxacsv[idx][labels.index("LEG_LEAN")]) / 1000.0)
    armLean_test.append(float(dxacsv[idx][labels.index("ARM_LEAN")]) / 1000.0)
    trunkLean_test.append(float(dxacsv[idx][labels.index("TRUNK_LEAN")]) / 1000.0)
    trunkFat_test.append(float(dxacsv[idx][labels.index("TRUNK_FAT")]) / 1000.0)
    legFat_test.append(float(dxacsv[idx][labels.index("LEG_FAT")]) / 1000.0)
    armFat_test.append(float(dxacsv[idx][labels.index("ARM_FAT")]) / 1000.0)
    currHeight_test = float(dxacsv[idx][labels.index("HEIGHT")]) / 100.0
    fMass_test = float(dxacsv[idx][labels.index("WBTOT_FAT")]) / 1000.0
    lMass_test = float(dxacsv[idx][labels.index("WBTOT_LEAN")]) / 1000.0
    FMI_test.append(fMass / (currHeight * currHeight))
    FFMI_test.append(lMass / (currHeight * currHeight))



#np.format_float_positional(np.mean())
print("Train")
print("Mean")
mean(age)
mean(height)
mean(mass)
mean(BMI)

mean(percentFat)
mean(leanMass)
mean(fatMass)
mean(viscFat)
mean(legLean)
mean(armLean)
mean(trunkLean)
mean(trunkFat)
mean(legFat)
mean(armFat)
mean(FMI)
mean(FFMI)
print("\nStandard Deviation")
standarddeviation(age)
standarddeviation(height)
standarddeviation(mass)
standarddeviation(BMI)

standarddeviation(percentFat)
standarddeviation(leanMass)
standarddeviation(fatMass)
standarddeviation(viscFat)
standarddeviation(legLean)
standarddeviation(armLean)
standarddeviation(trunkLean)
standarddeviation(trunkFat)
standarddeviation(legFat)
standarddeviation(armFat)
standarddeviation(FMI)
standarddeviation(FFMI)
print("\nMax")
max(age)
max(height)
max(mass)
max(BMI)

max(percentFat)
max(leanMass)
max(fatMass)
max(viscFat)
max(legLean)
max(armLean)
max(trunkLean)
max(trunkFat)
max(legFat)
max(armFat)
max(FMI)
max(FFMI)

print("\nMin")
min(age)
min(height)
min(mass)
min(BMI)

min(percentFat)
min(leanMass)
min(fatMass)
min(viscFat)
min(legLean)
min(armLean)
min(trunkLean)
min(trunkFat)
min(legFat)
min(armFat)
min(FMI)
min(FFMI)



print("\n\nTest")
print("\nMean")
mean(age_test)
mean(height_test)
mean(mass_test)
mean(BMI_test)

mean(percentFat_test)
mean(leanMass_test)
mean(fatMass_test)
mean(viscFat_test)
mean(legLean_test)
mean(armLean_test)
mean(trunkLean_test)
mean(trunkFat_test)
mean(legFat_test)
mean(armFat_test)
mean(FMI_test)
mean(FFMI_test)
print("\nStandard Deviation")
standarddeviation(age_test)
standarddeviation(height_test)
standarddeviation(mass_test)
standarddeviation(BMI_test)

standarddeviation(percentFat_test)
standarddeviation(leanMass_test)
standarddeviation(fatMass_test)
standarddeviation(viscFat_test)
standarddeviation(legLean_test)
standarddeviation(armLean_test)
standarddeviation(trunkLean_test)
standarddeviation(trunkFat_test)
standarddeviation(legFat_test)
standarddeviation(armFat_test)
standarddeviation(FMI_test)
standarddeviation(FFMI_test)
print("\nMax")
max(age_test)
max(height_test)
max(mass_test)
max(BMI_test)

max(percentFat_test)
max(leanMass_test)
max(fatMass_test)
max(viscFat_test)
max(legLean_test)
max(armLean_test)
max(trunkLean_test)
max(trunkFat_test)
max(legFat_test)
max(armFat_test)
max(FMI_test)
max(FFMI_test)

print("\nMin")
min(age_test)
min(height_test)
min(mass_test)
min(BMI_test)

min(percentFat_test)
min(leanMass_test)
min(fatMass_test)
min(viscFat_test)
min(legLean_test)
min(armLean_test)
min(trunkLean_test)
min(trunkFat_test)
min(legFat_test)
min(armFat_test)
min(FMI_test)
min(FFMI_test)

print("\npValue")
pval(age,age_test)
pval(height,height_test)
pval(mass,mass_test)
pval(BMI,BMI_test)

pval(percentFat,percentFat_test)
pval(leanMass,leanMass_test)
pval(fatMass,fatMass_test)
pval(viscFat,viscFat_test)
pval(legLean, legLean_test)
pval(armLean, armLean_test)
pval(trunkLean, trunkLean_test)
pval(trunkFat, trunkFat_test)
pval(legFat, legFat_test)
pval(armFat, armFat_test)
pval(FMI, FMI_test)
pval(FFMI, FFMI_test)





