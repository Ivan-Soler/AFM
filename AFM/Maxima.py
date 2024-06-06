import numpy as np

def simple(density,sizes):
    maxima=np.zeros((4, ), dtype=int)
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            for k in range(0,sizes[2]):
                for l in range(0,sizes[3]):
                    #First the nearest neighbours
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],k,l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,k,(l-1)%sizes[3]])):
                        continue
                    maxima=np.column_stack((maxima, [i,j,k,l]))
    maxima=np.delete(maxima, 0, 1)
    return(maxima)

def simple_2d(density,sizes):
    maxima=np.zeros((2, ), dtype=int)
    for i in range(0,sizes[0]):
        for j in range(0,sizes[3]):
            #First the nearest neighbours
            if ((density[i,j] < density [(i+1)%sizes[0],j]) or (density[i,j] < density [(i-1)%sizes[0],j])):
                continue
            if ((density[i,j] < density [i,(j+1)%sizes[3]]) or (density[i,j] < density [i,(j-1)%sizes[3]])):
                continue
            maxima=np.column_stack((maxima, [i,j]))
    maxima=np.delete(maxima, 0, 1)
    return(maxima)


def simple_1d(density,sizes):
    maxima=np.zeros((1, ), dtype=int)
    for i in range(0,sizes):
            #First the nearest neighbours
            if ((density[i] < density [(i+1)%sizes]) or (density[i] < density [(i-1)%sizes])):
                continue
            maxima=np.column_stack((maxima, [i]))
    if (len(maxima)>0):
        maxima=np.delete(maxima, 0)
    return(maxima)

def improve_1d(density,sizes):
    maxima=np.zeros((1, ), dtype=int)
    for i in range(0,sizes):
            #First the nearest neighbours
            if ((density[i] < density [(i+1)%sizes]) or (density[i] < density [(i-1)%sizes]) or (density[i] < density [(i+2)%sizes]) or (density[i] < density [(i-2)%sizes])):
                continue
            maxima=np.column_stack((maxima, [i]))
    if (len(maxima)>0):
        maxima=np.delete(maxima, 0)
    return(maxima)

def improve_2d(density,sizes):
    maxima=np.zeros((2, ), dtype=int)
    for i in range(0,sizes[0]):
        for j in range(0,sizes[3]):
                    #First the nearest neighbours
                    if ((density[i,j] < density [(i+1)%sizes[0],j]) or (density[i,j] < density [(i-1)%sizes[0],j])):
                        continue
                    if ((density[i,j] < density [i,(j+1)%sizes[3]]) or (density[i,j] < density [i,(j-1)%sizes[3]])):
                        continue
                    #Now the diagonals of i j
                    if ((density[i,j] < density [(i+1)%sizes[0],(j+1)%sizes[1]]) or (density[i,j] < density [(i-1)%sizes[0],(j-1)%sizes[1]])) or (density[i,j] < density [(i+1)%sizes[0],(j-1)%sizes[1]]) or (density[i,j] < density [(i-1)%sizes[0],(j+1)%sizes[1]]):
                        continue
                    maxima=np.column_stack((maxima, [i,j]))

    maxima=np.delete(maxima, 0, 1)
    return(maxima)
                   
def improve(density,sizes):
    maxima=np.zeros((4, ), dtype=int)
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            for k in range(0,sizes[2]):
                for l in range(0,sizes[3]):
                    #First the nearest neighbours
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],k,l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,k,(l-1)%sizes[3]])):
                        continue
                    #Now the diagonals of i j
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],(j+1)%sizes[1],k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],(j-1)%sizes[1],k,l])) or (density[i,j,k,l] < density [(i+1)%sizes[0],(j-1)%sizes[1],k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],(j+1)%sizes[1],k,l]):
                        continue
                    #Now the diagonals of i k    
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,(k-1)%sizes[2],l])) or (density[i,j,k,l] < density [(i+1)%sizes[0],j,(k-1)%sizes[2],l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,(k+1)%sizes[2],l]):
                        continue
                    #Now the diagonals of i l           
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,(l-1)%sizes[3]])) or (density[i,j,k,l] < density [(i+1)%sizes[0],j,k,(l-1)%sizes[3]]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,(l+1)%sizes[3]]):
                        continue
                    #Now the diagonals of j k     
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],(k-1)%sizes[2],l])) or (density[i,j,k,l] < density [i,(j+1)%sizes[1],(k-1)%sizes[2],l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],(k+1)%sizes[2],l]):
                    #Now the diagonals j l
                        continue
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,(l-1)%sizes[3]])) or (density[i,j,k,l] < density [i,(j+1)%sizes[1],k,(l-1)%sizes[3]]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,(l+1)%sizes[3]]):
                        continue
                    #Now the diagonals k l
                    if ((density[i,j,k,l] < density [i,j,(k+1)%sizes[2],(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],(l-1)%sizes[3]])) or (density[i,j,k,l] < density [i,j,(k+1)%sizes[2],(l-1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],(l+1)%sizes[3]]):
                        continue
                    maxima=np.column_stack((maxima, [i,j,k,l]))

    maxima=np.delete(maxima, 0, 1)
    return(maxima)

def nnneighbours(density,sizes):
    maxima=np.zeros((4, ), dtype=int)
    for i in range(0,sizes[0]):
        for j in range(0,sizes[1]):
            for k in range(0,sizes[2]):
                for l in range(0,sizes[3]):
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,l]) or (density[i,j,k,l] < density [(i+2)%sizes[0],j,k,l]) or (density[i,j,k,l] < density [(i-2)%sizes[0],j,k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],k,l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,l]) or (density[i,j,k,l] < density [i,(j+2)%sizes[1],k,l]) or (density[i,j,k,l] < density [i,(j-2)%sizes[1],k,l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],l]) or (density[i,j,k,l] < density [i,j,(k+2)%sizes[2],l]) or (density[i,j,k,l] < density [i,j,(k-2)%sizes[2],l])):
                        continue
                    if ((density[i,j,k,l] < density [i,j,k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,k,(l-1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,k,(l+2)%sizes[3]]) or (density[i,j,k,l] < density [i,j,k,(l-2)%sizes[3]])):
                        continue
                    #Now the diagonals of i j
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],(j+1)%sizes[1],k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],(j-1)%sizes[1],k,l])) or (density[i,j,k,l] < density [(i+1)%sizes[0],(j-1)%sizes[1],k,l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],(j+1)%sizes[1],k,l]):
                        continue
                    #Now the diagonals of i k    
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,(k-1)%sizes[2],l])) or (density[i,j,k,l] < density [(i+1)%sizes[0],j,(k-1)%sizes[2],l]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,(k+1)%sizes[2],l]):
                        continue
                    #Now the diagonals of i l           
                    if ((density[i,j,k,l] < density [(i+1)%sizes[0],j,k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,(l-1)%sizes[3]])) or (density[i,j,k,l] < density [(i+1)%sizes[0],j,k,(l-1)%sizes[3]]) or (density[i,j,k,l] < density [(i-1)%sizes[0],j,k,(l+1)%sizes[3]]):
                        continue
                    #Now the diagonals of j k     
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],(k+1)%sizes[2],l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],(k-1)%sizes[2],l])) or (density[i,j,k,l] < density [i,(j+1)%sizes[1],(k-1)%sizes[2],l]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],(k+1)%sizes[2],l]):
                    #Now the diagonals j l
                        continue
                    if ((density[i,j,k,l] < density [i,(j+1)%sizes[1],k,(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,(l-1)%sizes[3]])) or (density[i,j,k,l] < density [i,(j+1)%sizes[1],k,(l-1)%sizes[3]]) or (density[i,j,k,l] < density [i,(j-1)%sizes[1],k,(l+1)%sizes[3]]):
                        continue
                    #Now the diagonals k l
                    if ((density[i,j,k,l] < density [i,j,(k+1)%sizes[2],(l+1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],(l-1)%sizes[3]])) or (density[i,j,k,l] < density [i,j,(k+1)%sizes[2],(l-1)%sizes[3]]) or (density[i,j,k,l] < density [i,j,(k-1)%sizes[2],(l+1)%sizes[3]]):
                        continue
                    maxima=np.column_stack((maxima, [i,j,k,l]))
    #maxima=np.delete(maxima, 0, 1)
    return(maxima)

def repeated(res,sizes):
    distance=np.zeros([len(res[0]),len(res[0])])
    for i in range(0,len(res[0])):         
            for j in range(0,len(res[0])):
                for k in range(0,4):
                       distance[i][j]+=np.copy((np.abs(res[k][i]-res[k][j])%sizes[k])*(np.abs(res[k][i]-res[k][j])%sizes[k]))
                distance[i][j]=np.copy(np.sqrt(distance[i][j]))
    repeated=0
    for i in range(0,len(res[0])):         
        for j in range(0,len(res[0])):
            if (distance[i][j]<=2 and i<j):
                print(i,j)
                repeated+=1    
    print("repeated maxima = ",repeated)
    return()