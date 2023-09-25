#!/usr/bin/env python
# coding: utf-8

# In[1]:
import struct
import numpy as np

#class Projecting:
def next_element_4ind(i,j,k,l,sizes):
    l+=1
    if not (l%sizes[3]):
        k+=1
        l=0
        if not (k%sizes[2]):
            j+=1
            k=0
            if not (j%sizes[1]):
                i+=1
                j=0
            
    return(i,j,k,l)

def next_element_3ind(i,j,k,sizes):
    k+=1
    if not (k%sizes[2]):
        j+=1
        k=0
        if not (j%sizes[1]):
            i+=1
            j=0

    return(i,j,k)

def next_element_2ind(i,j,sizes):
    j+=1
    if not (j%sizes[1]):
        i+=1
        j=0

    return(i,j)

def D4_new(file):
    
    count=0
    i=0
    j=0
    k=0
    l=0
    
    sizes=[[]]*4
    with open(file,"rb") as file:
        for i in range (0,4):
            sizes[i]=int(file.readline())
        precision=int(file.readline())
        print(sizes,precision)
        #print(file.readline())
        elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]
        density=np.zeros((sizes[0],sizes[1],sizes[2],sizes[3]))
        print(file.readline())
        #file.readline() 
        #file.readline()
        #file.readline()
        #file.readline()
        #file.readline()
        #file.readline()
        for i in range(0, sizes[0]):
            for j in range(0, sizes[1]):
                for k in range(0, sizes[2]):
                    for l in range(0,sizes[3]):
                        num=struct.unpack('d', file.read(precision))
                        density[i,j,k,l]=np.copy(num[0])
                        #if (np.abs(num[0])>1e-18):
                           # print(num[0],i,j,k,l)
                        count+=1
    if (count==elements):
        print(count)
        return density

    if (count!=elements):
        print("Four dimensional array wrongly readed")
        print(count)
        return density


def D4(file,sizes,precision):
    
    elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]

    density=np.zeros((sizes[0],sizes[1],sizes[2],sizes[3]))

    count=0
    i=0
    j=0
    k=0
    l=0
    

    with open(file,"rb") as file:
        file.read(12) #skip the lattice size and precision
        while (i<sizes[0]):
                count+=1
                num=struct.unpack('d', file.read(precision))
                density[i,j,k,l]=num[0]
                i,j,k,l = next_element_4ind(i,j,k,l,sizes)
                
    if (count==elements):
        print("correct reading")
        return density

    if (count!=elements):
        print("Four dimensional array wrongly readed")
        return 0
    
    return 0
        
def Zero_mode(file,sizes,precision):    
    
    elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]
     
    density=np.zeros((sizes[0],sizes[1],sizes[2],sizes[3]))
    
    count=0
    i=0
    j=0
    k=0
    l=0
    
    with open('text.txt', 'w') as f:
        with open(file,"rb") as file:
            file.read(12) #skip the lattice size and precision
            while (i<sizes[0]):
                    count+=1
                    num=struct.unpack('d', file.read(precision))
                    density[i,j,k,l]=num[0]
                    i,j,k,l = next_element_4ind(i,j,k,l,sizes)
                    print(num[0], file=f)
                
    if (count==elements):
        print("correct reading")
        return density

    if (count!=elements):
        print("Four dimensional array wrongly readed")
        return 0
    
    return 0


def D4_D3 (file, direction, sizes, precision):

    elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]

    sizes_left=sizes.copy()
    sizes_left.pop(direction)

    density_projected=np.zeros((sizes_left[0],sizes_left[1],sizes_left[2]))

    count=0
    i=0
    j=0
    k=0
    l=0



    with open(file,"rb") as file:
        file.read(12) #skip the lattice size and precision
        if (direction==0):
            while (i<sizes[0]):
                count+=1
                num=struct.unpack('d', file.read(precision))
                density_projected[j,k,l]+=num[0]
                i,j,k,l = next_element_4ind(i,j,k,l,sizes)
                
    
        if (direction==1):
            while (i<sizes[0]):
                count+=1
                num=struct.unpack('d', file.read(precision))
                density_projected[i,k,l]+=num[0]
                i,j,k,l = next_element_4ind(i,j,k,l,sizes)
        
        if (direction==2):
            while (i<sizes[0]):
                count+=1
                num=struct.unpack('d', file.read(precision))
                density_projected[i,j,l]+=num[0]
                i,j,k,l = next_element_4ind(i,j,k,l,sizes)
                
        if (direction==3):
            while (i<sizes[0]):
                count+=1
                num=struct.unpack('d', file.read(precision))
                density_projected[i,j,k]+=num[0]
                i,j,k,l = next_element_4ind(i,j,k,l,sizes)

    if (count==elements):
        return density_projected, sizes_left

    if (count!=elements):
        print("bad 4to 3 projection")
        return 0

def D3_D2 (file, direction, sizes, precision):

    elements=sizes[0]*sizes[1]*sizes[2]

    sizes_left=sizes.copy()
    sizes_left.pop(direction)

    density_projected=np.zeros((sizes_left[0],sizes_left[1]))

    count=0
    i=0
    j=0
    k=0


    with open(file,"rb") as file:
        file.read(12) #skip the lattice size and precision   
        if (direction==0):
            while (i<sizes[0]):
                count+=1
                num=struct.unpack('d', file.read(precision))
                density_projected[j,k]+=num[0]
                i,j,k = next_element_3ind(i,j,k,sizes) 

        if (direction==1):
            while (i<sizes[0]):
                count+=1
                num=struct.unpack('d', file.read(precision))
                density_projected[i,k]+=num[0]
                i,j,k = next_element_3ind(i,j,k,sizes) 
                
        if (direction==2):
            while (i<sizes[0]):
                count+=1
                num=struct.unpack('d', file.read(precision))
                density_projected[i,j]+=num[0]
                i,j,k = next_element_3ind(i,j,k,sizes)  

    if (count==elements):
        return density_projected, sizes_left

    if (count!=elements):
        print("bad 3 to 2 projection")
        return 0

def D4_D2 (file, direction1, direction2, sizes, precision):

    if (direction1<=direction2):
        print("direction 1 should be bigger than direction 2")
        return 0

    density3D, sizes=D4_D3(file, direction1, sizes, precision)

    sizes_left=sizes.copy()
    sizes_left.pop(direction2)
    density_projected=np.zeros((sizes_left[0],sizes_left[1]))
    
    elements=sizes_left[0]*sizes_left[1]
    count=0

    if (direction2==0):
        for i in range(0,sizes_left[0]):
            for j in range(0,sizes_left[1]):
                count+=1
                density_projected[i,j]=density3D[:,i,j].sum()

    if (direction2==1):
        for i in range(0,sizes_left[0]):
            for j in range(0,sizes_left[1]):
                count+=1
                density_projected[i,j]=density3D[i,:,j].sum()

    if (direction2==2):
        for i in range(0,sizes_left[0]):
            for j in range(0,sizes_left[1]):
                count+=1
                density_projected[i,j]=density3D[i,j,:].sum()
    
    if (count==elements):
        return density_projected, sizes_left

    if (count!=elements):
        print("bad 4 to 2 projection")
        return 0
    
        
