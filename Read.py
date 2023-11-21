import numpy as np
import struct
import array
import os 

def ascii_mode(directory):
    sizes=[[]]*4
    count=0
    i=0
    j=0
    k=0
    l=0
    
    with open(directory,"r") as file:
        file.readline()
        file.readline()
        file.readline()
        info=file.readline()
        info=info.split(" ")

        for m in range (0,4):
            sizes[m]=int(info[m+2])
        range_color=int(info[0])
        range_spin=int(info[1])
        elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]   
        zeromode=np.zeros((range_color, range_spin, sizes[0],sizes[1],sizes[2],sizes[3]),dtype=np.complex_)
        for l in range(0, sizes[3]):
            for k in range(0, sizes[2]):
                for j in range(0, sizes[1]):
                    for i in range(0,sizes[0]):
                        for spin in range(0,range_spin):
                            for color in range(0,range_color):
                                count+=1
                                num=file.readline()
                                num=num.replace("(","")
                                num=num.replace(")","")
                                num=num.replace("\n","")
                                split=num.split(",")
                                zeromode[color,spin,i,j,k,l]=complex(float(split[0]),float(split[1]))
    density=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])

    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for spin in range(0,range_spin):
                        for color in range(0,range_color):                    
                            density[i,j,k,l]+=np.copy(zeromode[color,spin,i,j,k,l].imag*zeromode[color,spin,i,j,k,l].imag + 
                                                      zeromode[color,spin,i,j,k,l].real*zeromode[color,spin,i,j,k,l].real)
    norm=density.sum()                       
    return(zeromode,density,sizes)

def mode_to_density(zeromode,range_color,range_spin,sizes):
    density=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for spin in range(0,range_spin):
                        for color in range(0,range_color):
                            density[i,j,k,l]+=np.copy(zeromode[color,spin, i,j,k,l].imag*zeromode[color,spin, i,j,k,l].imag + 
                                                      zeromode[color,spin, i,j,k,l].real*zeromode[color,spin, i,j,k,l].real)          
    #density_1d=density.sum(axis=(0,1,2))
    return(density)


def mode_real_density(zeromode,range_color,range_spin,sizes,chirality):
    density=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])
    if chirality==1:
        spin=0
    if chirality==0:
        spin=2        
    for l in range(0, sizes[3]):
        for k in range(0, sizes[2]):
            for j in range(0, sizes[1]):
                for i in range(0,sizes[0]):
                    for color in range(0,range_color):
                            zeromode[color,spin,i,j,k,l]=np.real(zeromode[color,spin,i,j,k,l])
    
    density=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])
    
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for spin in range(0,range_spin):
                        for color in range(0,range_color):
                            density[i,j,k,l]+=np.copy(zeromode[spin,color,i,j,k,l].imag*zeromode[spin,color,i,j,k,l].imag + 
                                                      zeromode[spin,color,i,j,k,l].real*zeromode[spin,color,i,j,k,l].real)          
    density_1d=density.sum(axis=(0,1,2))
    return(density_1d)

def bin_mode(file,sizes,range_color,range_spin):

    mode=np.memmap(file, dtype=np.double).byteswap()

    density=np.zeros([sizes[0],sizes[1],sizes[2],sizes[3]])
    index=0
    for l in range(0, sizes[3]):
        for k in range(0, sizes[2]):
            for j in range(0, sizes[1]):
                 for i in range(0,sizes[0]):
                    for spin in range(0,range_spin):
                        for color in range(0,range_color):
                            density[i,j,k,l]+=mode[index]*mode[index]+mode[index+1]*mode[index+1]
                            index+=2    
    #os.remove(temp_file)
    return density,sizes

def scalar_product(mode1, mode2,sizes,range_color,range_spin):
    scalar=0
    for i in range(0, sizes[0]):
        for j in range(0, sizes[1]):
            for k in range(0, sizes[2]):
                for l in range(0,sizes[3]):
                    for color in range(0,range_color):
                        scalar+=mode1[color,0, i,j,k,l].conjugate()*mode2[color,2, i,j,k,l] + mode1[color,1, i,j,k,l].conjugate()*mode2[color,3, i,j,k,l]
    return np.sqrt(scalar),sizes

def ascii_mode_1d(namefile):
    density,sizes=ascii_mode(namefile)
    density_1d=density.sum(axis=(0,1,2))
    return(density_1d,sizes)

def bin_mode_1d(namefile,sizes,colors,spin):
    density,sizes=bin_mode(namefile,sizes,colors,spin)
    density_1d=density.sum(axis=(0,1,2))
    return(density_1d,sizes)

def topology(file):
    sizes=[[]]*4
    precision=0

    count=0
    i=0
    j=0
    k=0
    l=0
    norm=0
    with open(file,"rb") as file:
        for i in range (0,4):
            sizes[i]=int(file.readline())
        precision=int(file.readline())
        colors=int(file.readline())
        elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]
        density=np.zeros((sizes[0],sizes[1],sizes[2],sizes[3]))
        for i in range(0, sizes[0]):
            for j in range(0, sizes[1]):
                for k in range(0, sizes[2]):
                    for l in range(0,sizes[3]):
                        count+=1
                        num=struct.unpack('d', file.read(precision))
                        density[i,j,k,l]=num[0]
                        norm+=num[0]
    if (count==elements):
        return (density,sizes)

    if (count!=elements):
        print("Four dimensional array wrongly readed")
        return (density,sizes)


def topology_1d(file):
    density,sizes=topology(file)
    density_1d=density.sum(axis=(0,1,2))
    return(density_1d,sizes)

