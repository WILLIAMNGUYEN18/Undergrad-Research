import sys
import os
import numpy as np
import sklearn
import sklearn.metrics
import scipy.stats as stat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob

#strictly better version of computeR2 script from earlier. Cleaner because arrays are run through in a loop and you can pick any two metrics to compare

dxacsv = open("./dxa/SUQ4_DXA_200421.csv", 'r')
resultdir = glob.glob(sys.argv[1] + '/*')

#resultdir = glob.glob(sys.argv[1] + '/fit3d/test/*') + glob.glob(sys.argv[1] + '/styku/test/*') + glob.glob(sys.argv[1] + '/ss/test/*') #hack to look at all devices at once

idxA, idxB = int(sys.argv[2]), int(sys.argv[3]) #each item has 4 entries, pick which ones to plot
#0 is ground truth, 1 is PCA fit prediction, 2 is markered ganger fit projected prediction, 3 is refined ganger fit prediction

if sys.argv[4] == 'm':
    gender = "male"
elif sys.argv[4] == 'f':
    gender = "female"
else:
    print "BAD GENDER STRING"

print gender
comps = sys.argv[5]

printlabels = False

#column names
lines = dxacsv.readlines()
labels = lines.pop(0).split(',')
for l in range(0, len(labels)):
    labels[l] = labels[l].replace('"','')

sids = [l.split(',')[0] for l in lines]

percentfat, lmassp, fmassp, vfatp, legleanp, armleanp, trunkleanp, trunkfatp, legfatp, armfatp, FMIp, FFMIp = ([] for i in range(12))

#for each folder, pull the txt files representing the composition predictions ganger projection, the refined ganger projection, and the PCA prediction

#Compare them to each other or to the DXA file

for f in resultdir:
    subjid = f.split('/')[-1].split('_')[0][1:]
    print subjid
    line1 = lines[sids.index(subjid)].split(',')
    #print line1
    pctxt = open(f + '/d=-1/selfiereturnlinear.txt', 'r').readlines()
    
    compstr = '_' + comps
    if int(comps) < 0:
        compstr = ''
    #gangertxt = open(f + '/gangerprojectpredict.txt', 'r').readlines()
    refinedgangertxt = open(f + '/gangerrefinedprojectpredict' + compstr + '.txt', 'r').readlines()
    gangertxt = refinedgangertxt #hack this for now, index 2 doesnt do anything because not all of these have markers
    
    ht1 = float(line1[labels.index("HEIGHT")]) / 100.0
    wt1 = float(line1[labels.index("WBTOT_MASS")]) / 1000.0
    
    percentfat.append((float(line1[labels.index("WBTOT_PFAT")]), float(pctxt[2]) / wt1 * 100.0, float(gangertxt[8]) * 100.0, float(refinedgangertxt[3]) / wt1 * 100.0))

    lmass1 = float(line1[labels.index("WBTOT_LEAN")]) / 1000.0
    lmassp.append((lmass1, float(pctxt[1]), float(gangertxt[4]), float(refinedgangertxt[4])))

    fmass1 = float(line1[labels.index("WBTOT_FAT")]) / 1000.0
    fmassp.append((fmass1, float(pctxt[2]), float(gangertxt[3]), float(refinedgangertxt[3])))
    
    vfatp.append((float(line1[labels.index("VFAT_MASS")]) / 1000.0, float(pctxt[3]), float(gangertxt[5]), float(refinedgangertxt[5])))

    legleanp.append((float(line1[labels.index("LEG_LEAN")]) / 1000.0, float(pctxt[4]), float(gangertxt[7]), float(refinedgangertxt[7])))

    armleanp.append((float(line1[labels.index("ARM_LEAN")]) / 1000.0, float(pctxt[5]), float(gangertxt[6]), float(refinedgangertxt[6])))

    trunkfatp.append((float(line1[labels.index("TRUNK_FAT")]) / 1000.0, float(pctxt[9]), float(gangertxt[9]), float(refinedgangertxt[9])))

    trunkleanp.append((float(line1[labels.index("TRUNK_LEAN")]) / 1000.0, float(pctxt[10]), float(gangertxt[10]), float(refinedgangertxt[10])))
   
    armfatp.append((float(line1[labels.index("ARM_FAT")]) / 1000.0, float(pctxt[11]), float(gangertxt[11]), float(refinedgangertxt[11])))

    legfatp.append((float(line1[labels.index("LEG_FAT")]) / 1000.0, float(pctxt[12]), float(gangertxt[12]), float(refinedgangertxt[12])))


    fmi1 = fmass1 / (ht1 * ht1)

    ffmi1 = lmass1 / (ht1 * ht1)

    FMIp.append((fmi1,  float(pctxt[7]), float(gangertxt[8]), float(refinedgangertxt[3]) / (ht1 * ht1)))

    FFMIp.append((ffmi1,  float(pctxt[8]), float(gangertxt[9]), float(refinedgangertxt[4]) / (ht1 * ht1)))

def computeCV(pairs, A, B):
    sd = 0
    x = 0
    for j in pairs:
        sd_j = abs(j[A] - j[B]) / 2.0
        sd += pow(sd_j,2)
        x_j = (j[A] + j[B]) / 2.0
        x += x_j
    sd /= float(len(pairs))
    sd = np.sqrt(sd)
    x /= float(len(pairs))
    cv = sd / x
    return cv, sd

allarrays = [percentfat, lmassp, fmassp, vfatp, legleanp, armleanp, trunkleanp, trunkfatp, legfatp, armfatp, FMIp, FFMIp]
names = ["pfat", "leanmass", "fatmass", "visceralfat", "leglean", "armlean", "trunklean", "trunkfat", "legfat", "armfat", "FMI", "FFMI"]

ai = 0
sigs = 3
for a in allarrays:
    cv, sd = computeCV(a, idxA, idxB)
    array1 = [t[idxA] for t in a]
    array2 = [t[idxB] for t in a]
    if printlabels:
        print names[ai] + " %CV: " + str(round(cv * 100.0, sigs)) + "%" 
        print names[ai] + " RMSE: " + str(round(sd, sigs))
        print names[ai] + "R2: " + str(round(sklearn.metrics.r2_score(array1, array2), sigs))
        print "p-value " + str(round(stat.ttest_rel(array1, array2)[1], sigs)) + '\n'
    else:
        #print names[ai]
        #print round(cv * 100.0, sigs) 
        print round(sd, sigs)
        #print round(sklearn.metrics.r2_score(array1, array2), sigs)
        #print round(stat.ttest_rel(array1, array2)[1], sigs)

    #plot
    plt.clf()
    plt.scatter(array1, array2)
    plt.plot([0, max(array1) + 1],[0, max(array1)+ 1 ]) #y=x line
    plt.savefig('./refinedGangerCharts/' + names[ai] + '_' + gender + '_' + str(idxA) + str(idxB) + ".png")
    
    ai +=1
