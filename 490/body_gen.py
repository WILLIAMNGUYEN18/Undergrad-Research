import cPickle as pickle
import numpy as np
import math
import struct
import sys
import os 
import subprocess
from subprocess import Popen
from subprocess import PIPE


#https://stackoverflow.com/questions/7974849/how-can-i-make-one-python-file-run-another
#
#can import file, execfile('file.py') or os.system('python file.py')
#file.name_of_your_func()
#from subprocess import Popen
#Popen('python filename.py')
#os.system("python myOtherScript.py arg1 arg2 arg3")  
#Using os you can make calls directly to your terminal. If you want to be even more specific you can 
#concatenate your input string with local variables, ie.
#command = 'python myOtherScript.py ' + sys.argv[1] + ' ' + sys.argv[2]
#os.system(command)

print 'Number of arguments: ', len (sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
#os.path.dirname(os.path.abspath(__file__))


#sample count
n = 4

g = 0

#male loop?

mheight = np.linspace(1.6,1.9,17)


print ("starting male generations")

for height in mheight:
	for weight in range(65, 140, 1):
		#is sysargs going to be off by one since pca_to_shape may be included?
		command = 'python pca_to_shape.py ' + str(height) + ' ' + str(weight) 
		command += ' ' + str(n) + ' ' + str(g)
		print (command)
                #stdoutfile = open("stdout.txt", 'w')
                #stderrfile = open("stderr.txt", 'w')
		#jumbled mess piping to terminal
		print("Starting Popen for: " + str(height) + ' ' + str(weight) + ' ' + str(n) + ' ' + str(g))
		p = subprocess.Popen(command, shell=True)	
		#p.wait()
		print("Finished Popen for: " + str(height) + ' ' + str(weight) + ' ' + str(n) + ' ' + str(g))		




		

print ("starting female generations")
fheight = np.linspace(1.55,1.75, 15)
g = 1
#female loop?
for height in fheight:
	for weight in range(55, 140, 1):
		command = 'python pca_to_shape.py ' + str(height) + ' ' + str(weight) 
		command += ' ' + str(n) + ' ' + str(g)
		print (command)
		#jumbled mess piping to terminal
		print("Starting Popen for: " + str(height) + ' ' + str(weight) + ' ' + str(n) + ' ' + str(g))
		p = subprocess.Popen(command, shell=True)	
		print("Finished Popen for: " + str(height) + ' ' + str(weight) + ' ' + str(n) + ' ' + str(g))	
                #p.wait()


"""
    #run an instance of the selfie pipeline 

    pargs = " ".join(["./runselfie.sh", str(img), str(gender), str(subj_height), str(subj_wt), str(flength), str(hfov), str(vfov), str(pcomps)])
    print pargs
    processes.append(subprocess.Popen(pargs, stdout=stdoutfile, stderr=stderrfile, shell=True))

for p in processes:
    p.wait()

		#print("Waiting un previous subprocess")
                #p.wait()
		# of examples = 10
		# gender [0 (male) or 1 (female)]
		#https://stackoverflow.com/questions/3781851/run-a-python-script-from-another-python-script-passing-in-args
		#Popen('python pca_to_shape.py')
"""		
