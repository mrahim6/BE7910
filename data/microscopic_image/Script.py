# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 22:59:04 2019

@author: adil
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cv2
import os
import scipy
from scipy import ndimage
import matplotlib.pyplot as plt

path= 'd:/microscopic_image/'
names= os.listdir(path)

img= cv2.imread(path+names[1],0)
plt.imshow(img)

blobs, number_of_blobs = ndimage.label(img)
plt.imshow(blobs)


lists=[]
for i in range (1,number_of_blobs+1):
    sample=blobs[blobs==i]
    lists.append(len(sample))
    


df= pd.DataFrame(lists)
df.columns=['Values']
df1=pd.DataFrame()
df1['Particles'] = [i for i in range (1,number_of_blobs+1)]
df1['Pixel Count'] = df['Values']

df1.to_csv(path+'sample.csv',index=False)
np.savetxt(path+'sample_table.csv', blobs.astype(int), fmt='%i', delimiter=",")
cv2.imwrite(path+'blobs.png',blobs)


