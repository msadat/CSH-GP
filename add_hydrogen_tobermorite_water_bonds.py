#this script adds hydrogen with central oxygen (Xx) atoms and writes final lammps data file
# this script adds hydrogen with oxygen atoms
#input data 'xyz'
#output is lammps format data with the bond togology and connectivity

import numpy as np
import glob


filename = 'Tobermorite_ortho_AMC.xyz' # xyz file 
flist = glob.glob(filename)

for f in flist:
    load1 = np.genfromtxt(f, skip_header=2, dtype=float, usecols=(1,2,3)) #dtype=("|S10", float, float, float),
    data1=np.array(load1)

for f in flist:
    load2 = np.genfromtxt(f, skip_header=2, dtype=str, usecols=(0)) #dtype=("|S10", float, float, float),
    data2=np.array(load2)

nXx = 0 #no of water molecules
idx = [] #id of the water molecules 

for j in range(len(data2)):

  if data2[j] == 'Xx':
            nXx +=1
            idx.append(j)

coord = np.zeros((len(idx), 3)) #coordinates for the central Xx atoms

for a in range(len(idx)):
    coord[a] = data1[idx[a]]

nangles = nXx
nbonds = 2*nXx    
    
#for each point a vector ([180-104]/2 = 38 deg wrt x-axis) will be chosen and one H atoms will be added along that direction at a distance of 1 A.
#after this another H atom will be added at an angle of 104 degree with the first vector


def addatom(d, center):
    """d = distance to add the atom from central atom,
    center = array of coordinates for the central Xx atom where hydrogen is to be added"""
    
    rot_x = np.pi*180/38 # rotation wrt x axis (38 degree)
    bx1 = d*np.cos(rot_x)
    by1 = d*np.sin(rot_x)

    H_array1 = np.zeros((len(center),3)) #array for the coordinates of Hydrogen1
    H_array2 = np.zeros((len(center),3)) #array for the coordinates of Hydrogen2 
    
    for i in range(len(center)):
        H_array1 [i,0] = center[i,0] + bx1 
    for j in range(len(center)):  
        H_array1 [j,1] = center[j,1] + by1
    for k in range(len(center)):  
        H_array1 [k,2] = center[k,2] + 0 #no country for old man 'z'

    for i in range(len(center)):
        H_array2 [i,0] = center[i,0] - bx1 
    for j in range(len(center)):  
        H_array2 [j,1] = center[j,1] + by1
    for k in range(len(center)):  
        H_array2 [k,2] = center[k,2] + 0 #no country for old man 'z'
        
    return np.concatenate((H_array1,H_array2))

H_data = addatom(1, coord) #H_atom array
    
#*******************write final lammps data with hydrogen added with Ow (Xx)************************

xmin = np.min(data1[:,0]) #+min(0.0,xy,xz,xy+xz)
xmax = np.max(data1[:,0]) #+max(0.0,xy,xz,xy+xz)

ymin = np.min(data1[:,1]) #+min(0.0,yz)
ymax = np.max(data1[:,1]) #+max(0.0,yz)

zmin = np.min(data1[:,2])
zmax = np.max(data1[:,2])

a = 11.26500  
b = 7.38600   
c = 44.97000

natoms = len(data1)+len(H_data)

alpha = 90.0 
beta  = 90.0  
gamma = 90.0


lx=a
xy = b*np.cos(gamma*np.pi/180)
ly=np.sqrt(b**2 - xy**2)
xz = c*np.cos(beta*np.pi/180)
yz = (b*c*np.cos(alpha*np.pi/180)-xy*xz)/ly
lz = np.sqrt(c**2 - xz**2 - yz**2)

outFile = open('data.Tobermorite_ortho_water_bonds', 'w')
outFile.write('LAMMPS data file written by Rafat Sadat using Python\n')
outFile.write('\n')
outFile.write('%i %s \n' %(natoms, 'atoms'))
outFile.write('6 atom types \n')
outFile.write('%i %s \n' %(nbonds, 'bonds'))
outFile.write('1 bond types \n')
outFile.write('%i %s \n' %(nangles, 'angles'))
outFile.write('1 angle types \n')
outFile.write('\n')
outFile.write('%f %f %s %s \n' %(-1.5, lx, 'xlo', 'xhi'))
outFile.write('%f %f %s %s \n' %(-1.5, ly, 'ylo', 'yhi'))
outFile.write('%f %f %s %s \n' %(-1.5, lz, 'zlo', 'zhi'))
outFile.write('\n')
outFile.write('%f %f %f %s %s %s \n' %(xy, xz, yz, 'xy', 'xz', 'yz'))
outFile.write('\n')
outFile.write('Atoms\n')
outFile.write('\n')


for j in range(len(data1)): #writing atoms without water
    if data2[j]=="Si":
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 1, 0, data1[j,0], data1[j,1], data1[j,2]))
    elif data2[j]=="O":
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 4, 0, data1[j,0], data1[j,1], data1[j,2]))
    elif data2[j]=="Ca":
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 3, 0, data1[j,0], data1[j,1], data1[j,2]))
    elif data2[j]=='Xx': #Ow
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 5, 0, data1[j,0], data1[j,1], data1[j,2]))     


for j in range(len(H_data)):
    outFile.write('%i %i %i %i %f %f %f \n' %(len(data1)+j+1, 0, 6, 0, H_data[j,0], H_data[j,1], H_data[j,2]))
    
              
outFile.write('\n')
outFile.write('Bonds\n')
outFile.write('\n')

count = 0
for j in range(0, len(idx), 1): #idx is the list of Xx atom serial no.
    count +=1
    outFile.write('%i %s %i %i \n' %(count,'1',idx[j]+1,len(data1)+j+1))
    count +=1
    outFile.write('%i %s %i %i  \n' %(count,'1',idx[j]+1,len(data1)+len(H_data)/2+j+1))

outFile.write('\n')
outFile.write('Angles\n')
outFile.write('\n')   

count = 0
for j in range(0, len(idx), 1):
    count +=1
    outFile.write('%i %s %i %i %i \n' %(count,'1', len(data1)+j+1, idx[j]+1, len(data1)+len(H_data)/2+j+1))    

outFile.close()
print "All done!"   

