#conda install pytorch torchvision -c pytorch
#pip3 install https://download.pytorch.org/whl/cu90/torch-1.0.0-cp37-cp37m-win_amd64.whl
#pip3 install torchvision
# This is code to download and install pytorch
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