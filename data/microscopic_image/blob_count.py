import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import cv2
#import os
from scipy import ndimage

#blob_numbers = []
def blob_count (path_to_image, image_name, per_pixel_length, path_to_save):
    """
    inputs
        path_to_image = path where the image is located
        image_name = name of the image with proper extensions
        per_pixel_length = length of reference line / pixels covered by reference line
        path_to_save = save directory for the dataframe
        
    outputs
        dataframe with area and diameter (saved in path_to_save)
    """   
    img= cv2.imread(path_to_image + "/" + image_name, 0)
    blobs, number_of_blobs = ndimage.label(img)
    #blob_numbers.append(number_of_blobs)
    lists=[]
    for i in range (1,number_of_blobs+1):
        sample=blobs[blobs==i]
        lists.append(len(sample))
    df= pd.DataFrame(lists)
    df.columns=['Values']
    df1=pd.DataFrame()
    df1['Particles'] = [i for i in range (1,number_of_blobs+1)]
    df1['Pixel Count'] = df['Values']    
    df1['area']= np.square(per_pixel_length) * df1['Pixel Count'] # Area calculation 
    #df1['area']= np.square((99.9/719.0)) * df1['Pixel Count']  , length of reference line = 99.9 um, pixels covered = 719
    df1['dia']= np.sqrt((4/np.pi) * df1['area'])  # Diameter calculation
    df1.to_csv(path_to_save+image_name[:-4]+'.csv',index=False)

#example
path = 'd:/microscopic_image/2day_segmented/sample_2/'
blob_count (path,'layer1.png', 99.9/719 , 'd:/')

