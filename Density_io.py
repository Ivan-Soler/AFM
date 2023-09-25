import struct
import sys

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

def D4(file,sizes,precision):
    elements=sizes[0]*sizes[1]*sizes[2]*sizes[3]

    density = [[[[0 for x in range(sizes[0])] for y in range(sizes[1])] for z in range(sizes[2])] for t in range(sizes[3])]
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
                density[i][j][k][l]=num[0]
                i,j,k,l = next_element_4ind(i,j,k,l,sizes)
                
    if (count==elements):
        return density

    if (count!=elements):
        print("Four dimensional array wrongly readed")
        return 0

def writting_tony(data,file_name,sizes,precision):
    file=open(file_name,"w")
    for t in range(0,sizes[3]):
        for x in range(0,sizes[0]):
            for y in range(0,sizes[1]):
                for z in range(0,sizes[2]):
                    file.write(str(t)+"\t"+str(x)+"\t"+str(y)+"\t"+str(z)+"\t"+str(data[x][y][z][t])+"\n")  
    file.close()
    

#MAIN program
if (len(sys.argv)>1):   
    file_name_1=str(sys.argv[1])
else:
    print("Error: provide name of file")
    sys.exit()
    
sizes=[0,0,0,0]

f = open(file_name_1, "rb")
for i in range (0,4):
    sizes[i]=int(f.readline())
precision=int(f.readline())
f.close()

output_file="density.dat"

density=D4(file_name_1,sizes,precision)
writting_tony(density,output_file,sizes,precision)