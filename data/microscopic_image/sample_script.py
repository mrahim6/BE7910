import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cv2
import os
from scipy import ndimage


path= 'd:/microscopic_image/2day_segmented/sample_2/'
names= os.listdir(path)

finals= []
for a in names:
    img= cv2.imread(path+a,0)
    blobs, number_of_blobs = ndimage.label(img)
    finals.append(number_of_blobs)
    lists=[]
    for i in range (1,number_of_blobs+1):
        sample=blobs[blobs==i]
        lists.append(len(sample))
    df= pd.DataFrame(lists)
    df.columns=['Values']
    df1=pd.DataFrame()
    df1['Particles'] = [i for i in range (1,number_of_blobs+1)]
    df1['Pixel Count'] = df['Values']    
    df1['area']= np.square((99.9/719.0)) * df1['Pixel Count']
    df1['dia']= np.sqrt((4/np.pi) * df1['area'])
    df1.to_csv(path+a[:-4]+'.csv',index=False)



