# -*- coding: utf-8 -*-
"""BodyStatsNet_MSE.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1kAVQwALHaYfmCKVc3VRhn5ChpfBv6HIx
"""

#install pytorch
import os
from wheel.pep425tags import get_abbr_impl, get_impl_ver, get_abi_tag
platform = '{}{}-{}'.format(get_abbr_impl(), get_impl_ver(), get_abi_tag())

accelerator = 'cu80' if os.path.exists('/opt/bin/nvidia-smi') else 'cpu'
!pip install http://download.pytorch.org/whl/{accelerator}/torch-0.4.1-{platform}-linux_x86_64.whl torchvision

import torch
print('Version', torch.__version__)
print('CUDA enabled:', torch.cuda.is_available())

!apt-get -qq install -y libsm6 libxext6 && pip install -q -U opencv-python

import cv2

!pip install -U scikit-learn

# Connect Google Drive to Colab
# Load the Drive helper and mount
from google.colab import drive

# This will prompt for authorization.
drive.mount('/fdrive')
# Create a directory and mount Google Drive using that directory.

#setting up path
gender = "male"
#gender = "female"
BASE_PATH = '/fdrive/My Drive/Research/bodystats/'
if not os.path.exists(BASE_PATH):
  os.makedirs(BASE_PATH)
DATA_PATH = BASE_PATH

os.chdir('/fdrive/My Drive/Research/bodystats')
!ls
os.chdir('/content')

#imports
import torch
import torch.nn as nn
from torchvision import datasets
from torchvision import transforms

import math
import numpy as np
import os
import torch.nn.functional as F
import torch.optim as optim
import glob
import sys
import sklearn
import sklearn.metrics
import matplotlib
import matplotlib.pyplot as plt
import math
import copy

#print("CURRENT DIRECTORY")
#print(os.getcwd())
#!ls
sys.path.append(BASE_PATH)
#!ls

#print("CURRENT DIRECTORY")
#print(os.getcwd())
print(BASE_PATH)
#BREAKING HERE, CAN'T FIND PT_UTIL
import pt_util

"""DEFINING THE NEURAL NETWORK"""

LOSSFN = 'mse'
#Can keep core structure
from torchvision import models

#Can have same structure as BodyCompNet
#Need to have nn.linear instead of nn.Module
#https://stackoverflow.com/questions/45869131/all-tensorflow-outputs-are-nan
class BodyStatsNet(nn.Module):
    bestacc = -1.0

    def __init__(self):
        super(BodyStatsNet, self).__init__()
        meep = 1
        #if LOSSFN == 'crosse':
        #  meep = 200
        #80 input to 128 node 2nd layer  
        #self.fc1 = nn.Linear(80,32)
        #self.drop = nn.Dropout()
        
        #128 node 2nd layer to 200 node output, used for classification in this case
        
        #potentially reduce middle layers due to overfitting
        #https://stats.stackexchange.com/questions/181/how-to-choose-the-number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw
        self.fc2 = nn.Linear(80, 64)
        self.fc3 = nn.Linear(64, meep)

        
    def forward(self, x):
        # TODO define the forward pass
        #this is the method to define how layers pass through each other
        #input goes into fully connected layer 1
        #then that output goes into fully connected layer 2
        #the result of that is returned
        #x = F.relu(self.fc1(x))
        #x = F.relu(self.fc2(x))
        
        #x = F.relu(self.fc1(x))
        
        #x = F.relu(self.fc1(x))
        #x = self.drop(x)
        x = self.fc2(x)
        x = F.relu(self.fc3(x))
        #Hidden layers have too many features that effectively memorize the training set
        #Overfitting
        #need to add dropout layers that effectively remove weights, removing features to reduce overfitting
       #Isaac's example had a resnet, and the last layer of resnet is a fully connected layer
        return x
        
    def loss(self, prediction, label, reduction='elementwise_mean'):
        loss_val = None

        #classification is done in the loss function
        #200 bins, 1 if correctly classified, 0 if not, show loss
        #if LOSSFN == 'crosse':
        #    loss_val = F.cross_entropy(prediction, label.squeeze(), reduction=reduction)
        if LOSSFN == 'mse':
            loss_val = F.mse_loss(prediction.squeeze(), label.squeeze())

        return loss_val

    def save_model(self, file_path, num_to_keep=1):
        pt_util.save(self, file_path, num_to_keep)
        
    def save_best_model(self, accuracy, file_path, num_to_keep=1):
        # TODO save the model if it is the best
        if accuracy > BodyStatsNet.bestacc:
          pt_util.save(self, file_path, num_to_keep)
        
    def load_model(self, file_path):
        pt_util.restore(self, file_path)

    def load_last_model(self, dir_path):
        return pt_util.restore_latest(self, dir_path)

"""Data Loader:
Load in Body Statistics


> GetItem == People with their list of stats? OR people's DXA

> labels = % fat

> LEN: # of people (number of PCA values from CSV)
"""

#Will's Data Loader
class BodyStatsDataset(torch.utils.data.Dataset):
    def __init__(self, datadir, transform=None):
      # TODO Implement data loading.
      #defin
      self.transform = transform
      self.datadir = datadir
      dxacsv = open(datadir)
         
      #print("LINES:")
      self.lines = dxacsv.readlines()
      #print(self.lines)
      
      
      #print("LABELS:")
      self.labels = self.lines.pop(0).split(',')
      for l in range(0, len(self.labels)):
        self.labels[l] = self.labels[l].replace('"','')
      #Adjust spreadsheet to include all PCA vectors
      #slice manually at 80
      self.labels[-1] = "PC80"
      #print(self.labels)
      
      #print("SUBJECT IDs:")
      #self.sids = [l.split(',')[self.labels.index("SubjectID")] for l in self.lines]
      #print(self.sids)
      

      #self.filelist = glob.glob(datadir + '/*')
    def __len__(self):
      #number of subfolders = number of subjects
      #return len(self.filelist)
      # of rows -1
      #print(len(self.lines))
      return len(self.lines)
        
        
    #index ith index
    #pulls out ith object in vector
    def __getitem__(self, idx):
      PCAvec = []
      lbls = 0
      #print("SIZE")
      #print(len(self.lines))
      
      #print("FLAG 1")
      #print(idx)
      #dxacsv = open(datadir)
      #NEED TO HARD DODGE if 0?
      #if idx == 0:
      #    print("What the fuck")
      #else:    
      #print("DXA LINE")
      #dxaline = self.lines[idx]
      #not split by commas, it's a string atm
      #print("FLAG 2")
      dxaline = self.lines[idx].split(',')
      #print(dxaline)
      
      
      #temp = dxaline[-1]

      #temp = temp.rstrip()
      #print("LAST CHARACTER OF DXALINE")
      #print(str(temp))
      #print(len(temp))
      #dxaline[-1] = temp

      #print("LABELS1:")
      #print(self.labels)

      #print("FLAG 3")
      #print("INDEX OF WBTOT_FAT")
      #print(self.labels.index("WBTOT_FAT"))
      #WBTOT_FAT has quotes around it in later iteration?
      #for l in self.labels:
      #    l = l.rstrip()
      #fmass = (float(dxaline[self.labels.index('\"WBTOT_FAT\"')])/1000.0)
      #lmass = (float(dxaline[self.labels.index('\"WBTOT_LEAN\"')])/1000.0)
      #print("LABELS2:")
      #print(self.labels)
      fmass = (float(dxaline[self.labels.index('WBTOT_FAT')])/1000.0)
      lmass = (float(dxaline[self.labels.index('WBTOT_LEAN')])/1000.0)
      #sca(float(dxaline[labels.index("WBTOT_FAT")])/1000.0)
      pfat = (fmass / (fmass + lmass)) * 100.0
      #print("PFAT:" + str(pfat))
      #lbls = pfat
      #print("FLAG 4")


      for i in range(1,81,1):
          currPCA = "PC" + str(i)
          #print(str(currPCA))
          #print(str(float(dxaline[self.labels.index(currPCA)])))
          #print("Type of PCAval " + str(type(float(dxaline[self.labels.index(currPCA)]))))
          PCAvec.append(float(dxaline[self.labels.index(currPCA)]))
      
      #print("FLAG 5")
      #pfatcalc = pfat * 100.0
      #pfatcalc = torch.tensor(pfat * 100.0, dtype=torch.float)
      #pfatcalc = torch.tensor(pfat, dtype=torch.float)
      pfatcalc = torch.tensor( 100.0 * min(round( (4.0 * pfat)), 200.0))
      #print("PFAT: " + str(pfat))
      #print("tensor of PFAT: " + str(pfatcalc))

      #PCAvec = torch.tensor(PCAvec, dtype = torch.float)
      lbls = pfatcalc
      PCAvecs = torch.tensor(PCAvec, dtype = torch.float)
      #print("RETURN LABEL")
      #print(str(lbls))
      #print(PCAvecs)
      #print("FLAG 6")
      #CHECK TYPES OF PCAvec, PCAvec values, and lbls
      #print("Type of PCAvec" + str(type(PCAvecs)))
      #print("Type of lbls" + str(type(lbls)))
      return (PCAvecs, lbls)

print(os.path.exists(DATA_PATH + "shapeup_q2_8020split_dxa_3d_"+ gender +"_train.csv"))
print(os.path.exists(DATA_PATH + "dxa_pcaweights_"+ gender +"_valid.csv"))
      
data_train = BodyStatsDataset(DATA_PATH + "shapeup_q2_8020split_dxa_3d_"+ gender +"_train.csv")

data_test = BodyStatsDataset(DATA_PATH + "dxa_pcaweights_"+ gender +"_valid.csv")

"""Train and Test can be the same here as the internal mechanics of the neural net are unchanged (effectively a black box)."""

import time

def train(model, device, train_loader, optimizer, epoch, log_interval):
    model.train()
    lastloss = -1.0
    sse = 0
    gt = []
    prediction = []
    
    for batch_idx, (data, label) in enumerate(train_loader):
        #print(data.shape)
        #print(label.shape)
        data, label = data.to(device), label.to(device)
        #label = label.long()
        optimizer.zero_grad()
        output = model(data.float())
        #print(data.shape)
        #print(output.shape)
        #print(label.shape)
        loss = model.loss(output, label)
        #print("Train Batch Predictions: " + str(output.tolist()))
        loss.backward()
        optimizer.step()
        lastloss = loss
        
        #pred = output.max(1)[1]     #indices
        #li = 0
        #lablist = label.squeeze().tolist()
        #for p in pred.squeeze().tolist():
        #    l = lablist[li]
        #    p_pfat = p / 4.0
        #    l_pfat = l / 4.0
        #    gt.append(l_pfat)
        #    prediction.append(p_pfat)
        #    sse += math.pow(p_pfat - l_pfat, 2)
        #    li += 1
        
        
        #if mse
        #output.item() instead of max to print out the singular output that contains all values
        #
        #p instead of p /4.0
        #l instead of l /4.0
        pred = output     #indices
        li = 0
        lablist = label.squeeze().tolist()
        print("Output of Pred and Output during TRAIN")
        #print("output.")
        
        for p in pred.squeeze().tolist():
            l = lablist[li]
            p_pfat = p
            l_pfat = l
            #print("p_fat: ")
            #print(p_pfat)
            #print("l: ")
            #print(l_pfat)
            gt.append(l_pfat)
            prediction.append(p_pfat)
            sse += math.pow(p_pfat - l_pfat, 2)
            li += 1
         #if LOSSFN == 'mse':
         #    pred = output
         #    lablist = label.squeeze().tolist()
         #    li = 0
             #for p in pred.squeeze().tolist():
         #         l = lablist[li]
         #         regr_results.append((p, l))
         #         gt.append(l)
         #         prediction.append(p)
         #         li += 1    
        
        print("Printing GT and Prediction")
        print(gt)
        print(prediction)
        if batch_idx % log_interval == 0:
            print('{} Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                time.ctime(time.time()),
                epoch, batch_idx * len(data), len(train_loader.dataset),
                100. * batch_idx / len(train_loader), loss.item()))
            
    print("epoch: " + str(epoch))        
    rmse = math.sqrt(sse / len(train_loader.dataset))
    r2val = sklearn.metrics.r2_score(gt, prediction)
    print("Train R2: " + str(r2val))
    print("Train RMSE: " + str(rmse))
    plt.clf()
    plt.scatter(gt, prediction)
    return lastloss
            
def test(model, device, test_loader, log_interval=None):
    model.eval()
    test_loss = 0
    correct = 0
    
    sse = 0
    
    regr_results = []
    
    gt = []
    prediction = []
    
    with torch.no_grad():
        for batch_idx, (data, label) in enumerate(test_loader):
            data, label = data.to(device), label.to(device)
            #print(label)
            #print(data.shape)
            #label = label.long()
            
            output = model(data.float())
            
            #print(label.size(0))
            #print(output.size(0))
            test_loss_on = model.loss(output, label, reduction='sum').item()
            test_loss += test_loss_on
            
            #print(LOSSFN)
            pred = None
            correct_mask = None
            num_correct = None
            #cross entropy prediction
            if LOSSFN == 'crosse':
                pred = output.max(1)[1]     #indices           
                correct_mask = pred.eq(label.view_as(pred))
                num_correct = correct_mask.sum().item()
                correct += num_correct
                li = 0
                lablist = label.squeeze().tolist()
                for p in pred.squeeze().tolist():
                  l = lablist[li]
                  
                  p_pfat = p / 4.0
                  l_pfat = l / 4.0
                  regr_results.append((p_pfat, l_pfat))
                  gt.append(l_pfat)
                  prediction.append(p_pfat)
                  sse += math.pow(p_pfat - l_pfat, 2)
                  li += 1
            
            #mse
            if LOSSFN == 'mse':
                #https://pytorch.org/docs/stable/tensors.html
                #Use torch.Tensor.item() to get a Python number from a tensor containing a single value:
                #only one element tensors can be converted to Python scalars
                #so output has more than one value currently even though output is 1?
                print("Printing Output")
                print(output)
                #pred = output.item()
                pred = output
                print("Output/pred is type:")
                print(type(output))
                lablist = label.squeeze().tolist()
                li = 0
                for p in pred.squeeze().tolist():
                  l = lablist[li]
                  regr_results.append((p, l))
                  gt.append(l)
                  #VALUES HERE ARE NaN
                  #need to check if pred is correct
                  #output.item?
                  prediction.append(p)
                  li += 1
            
            if log_interval is not None and batch_idx % log_interval == 0:
                print('{} Test: [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                    time.ctime(time.time()),
                    batch_idx * len(data), len(test_loader.dataset),
                    100. * batch_idx / len(test_loader), test_loss_on))
 
    if LOSSFN == 'crosse':
        test_loss /= len(test_loader.dataset)
        test_accuracy = 100. * correct / len(test_loader.dataset)
        print('\nTest set: Average loss: {:.4f}, Accuracy: {}/{} ({:.0f}%)\n'.format(
            test_loss, correct, len(test_loader.dataset), test_accuracy))
        
        rmse = math.sqrt(sse / len(test_loader.dataset))
        print('RMSE: {:.4f}\n'.format(rmse))
        
        #for i in range(0, len(regr_results)):
        #  print('Test {}: Predicted: {:.4f}, Actual: {:.4f}\n'.format(
        #    i, regr_results[i][0], regr_results[i][1]))
        #print("\nPredicted:")
        #for i in range(0, len(regr_results)):
        #  print(regr_results[i][0])
        #print("\nActual:")
        #for i in range(0, len(regr_results)):
        #  print(regr_results[i][1])
        
        r2val = sklearn.metrics.r2_score(gt, prediction)
        print("R2: " + str(r2val))
        plt.clf()
        plt.scatter(gt, prediction)
        
        return test_loss, rmse
    if LOSSFN == 'mse':
        rmse = math.sqrt(test_loss)
        print('\nTest set: Average loss: {:.4f}, RMSE: {:.4f}\n'.format(
            test_loss, rmse))
        
        #print out predicted vs actual %fat
        for i in range(0, len(regr_results)):
          #print("fuck")
          #issue on this line. Printing NaN or INF
          print('Test {}: Predicted: {:.4f}, Actual: {:.4f}\n'.format(
            i, regr_results[i][0], regr_results[i][1]))
        print("Printing GT and Pred for TEST")
        print(gt)
        print(prediction)
        r2val = sklearn.metrics.r2_score(gt, prediction)
        print("R2: " + str(r2val))
        plt.clf()
        plt.scatter(gt, prediction)
        
        return test_loss, rmse

# Play around with these constants, you may find a better setting.
REALT_BATCH_SIZE = 147 #17 #36 #255
SYNTHT_BATCH_SIZE = 102
TEST_BATCH_SIZE = 38 #31 for male 39 for female
EPOCHS = 100 #50
LEARNING_RATE = 5e-7
MOMENTUM = 0.9
USE_CUDA = True
PRINT_INTERVAL = 5
WEIGHT_DECAY = 5e-7
LOG_PATH = DATA_PATH + '_log.pkl'
CHECKPOINT_PATH = DATA_PATH + '/checkpoints'


#https://stackoverflow.com/questions/45869131/all-tensorflow-outputs-are-nan
#lower learning rate to deal with alternating negatives that lead to NaN


# Now the actual training code
use_cuda = USE_CUDA and torch.cuda.is_available()

device = torch.device("cuda" if use_cuda else "cpu") #not here if you get cuda error
print('Using device', device)
import multiprocessing
print('num cpus:', multiprocessing.cpu_count())

kwargs = {'num_workers': multiprocessing.cpu_count(),
          'pin_memory': True} if use_cuda else {}

#train_loader = torch.utils.data.DataLoader(data_train_synth, batch_size=SYNTHT_BATCH_SIZE,
#                                           shuffle=True, **kwargs)
realtrain_loader = torch.utils.data.DataLoader(data_train, batch_size=REALT_BATCH_SIZE,
                                           shuffle=True, **kwargs)
test_loader = torch.utils.data.DataLoader(data_test, batch_size=TEST_BATCH_SIZE,
                                          shuffle=False, **kwargs)

#synthtrain80_loader = torch.utils.data.DataLoader(data_train80_synth, batch_size=SYNTHT_BATCH_SIZE,
#                                           shuffle=True, **kwargs)
#synthtest20_loader = torch.utils.data.DataLoader(data_test20_synth, batch_size=TEST_BATCH_SIZE,
#                                          shuffle=False, **kwargs)

model = BodyStatsNet().to(device)
optimizer = optim.SGD(model.parameters(), lr=LEARNING_RATE, momentum=MOMENTUM, weight_decay=WEIGHT_DECAY)
start_epoch = model.load_last_model(DATA_PATH + 'checkpoints')
print("Start: " + str(start_epoch))
#start_epoch = 0

# You may want to define another default for your log data depending on how you save it.
log_data = pt_util.read_log(LOG_PATH, None)

#matrix somewhere with NaN or INF
#https://stackoverflow.com/questions/49013049/runtimewarning-invalid-value-encountered-in-reduce
losstest, acctest = test(model, device, test_loader, True)

train_loss = []
test_loss = []
test_acc = []
epochcount = []

if log_data != None:
  train_loss = log_data[0]
  test_loss = log_data[1]
  test_acc = log_data[2]
  epochcount = log_data[3]
  
currepoch = 0

print(epochcount)
try:
    for epoch in range(start_epoch, start_epoch + EPOCHS + 1):
        losstrain = train(model, device, realtrain_loader, optimizer, epoch, PRINT_INTERVAL) #change train loader here
        losstest, acctest = test(model, device, test_loader) #change t est loader here
        model.save_best_model(acctest, CHECKPOINT_PATH + '/%03d.pt' % epoch)
        train_loss.append(losstrain)
        test_loss.append(losstest)
        test_acc.append(acctest)
        epochcount.append(epoch)
        currepoch = epoch
        print(epochcount)

    #dont train on real photos
    ##train again on the real photos
    #for epoch in range(start_epoch + EPOCHS + 1, start_epoch + 3 * EPOCHS + 1):
    #    losstrain = train(model, device, realtrain_loader, optimizer, epoch, REALT_BATCH_SIZE)
    #    losstest, acctest = test(model, device, test_loader)
    #    model.save_best_model(acctest, CHECKPOINT_PATH + '/%03d.pt' % epoch)
    #    train_loss.append(losstrain)
    #    test_loss.append(losstest)
    #    test_acc.append(acctest)
    #    epochcount.append(epoch)
    #    currepoch = epoch
    #    print(epochcount)
        
except KeyboardInterrupt as ke:
    print('Interrupted')
except:
    import traceback
    traceback.print_exc()
finally:
    # Always save the most recent model, but don't delete any existing ones.
    model.save_model(DATA_PATH + 'checkpoints/%03d.pt' % (currepoch), 0)     
    pt_util.write_log(LOG_PATH, [train_loss, test_loss, test_acc, epochcount])
    
    #plot shit
    pt_util.plot(epochcount, train_loss, "Training Loss", 'Epoch', 'Train Loss')
    pt_util.plot(epochcount, test_loss, "Test Loss", 'Epoch', 'Test Loss')
    pt_util.plot(epochcount, test_acc, "Test Accuracy", 'Epoch', 'Test Accuracy')

#plot shit
    pt_util.plot(epochcount[:-10], train_loss[:-10], "Training Loss", 'Epoch', 'Train Loss')
    pt_util.plot(epochcount[:-10], test_loss[:-10], "Test Loss", 'Epoch', 'Test Loss')
    pt_util.plot(epochcount[:-10], test_acc[:-10], "Test Accuracy", 'Epoch', 'Test Accuracy')

"""#cross entropy is classification
          last layer is 200 bins
          
          each line is a row
          split by commas,
          access column by label
          
          
          
          For every person in training set (subjectID)
          Rows in CSV - 1
          Feature is PCA vector (first 80 pca vectors)
          index of "PC1" - "PC80"
          Input is size 80, pulled out of CSV for each person
          Label = some derivative of fat (pfat)
          ROUGHLYselfiecomput float(dxaline label.index "wbtotalfat" /1000)
          Fully connected network
          nn.linear
          80 first, 200 end
          
          
          FOR TEST SET:
          no vectors at end of valid
          CSVs
          2 training
          2 validation
          2 wonky
          Identified by MeshID
          
          Open up DXA file, open up weights
          for each person in DXA file, look up meshID in PCA weights
          
          OR
          
          can sort both and combine them
          
          fora data loader, return a vector tensor and label (scalar)x
"""