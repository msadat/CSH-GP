#converting from lammps to xyz for Metakaolinite
#Author M. R. Sadat
#use this shell command:
# awk '{gsub("[0-9]*","",$1)}1' Lizardite_555.xyz  > Lizardite_corrected.xyz

import numpy as np
import glob


filename2 = 'data.CSH-nano-4nm' # lammps data file
flist2 = glob.glob(filename2)


natoms =  2944  

for f in flist2:
    #load2 = np.genfromtxt(f, dtype=float, skip_header =19, usecols=(4,5,6)) #dtype=("|S10", float, float, float),
    load2 = np.genfromtxt(f, dtype=float, skip_header =19, skip_footer = natoms+1, usecols=(4,5,6)) 
    dataovito1=np.array(load2)

for f in flist2:
    #load2 = np.genfromtxt(f, dtype=int, skip_header =19, usecols=(2)) #dtype=("|S10", float, float, float),
    load2 = np.genfromtxt(f, dtype=int, skip_header =19, skip_footer = natoms+1, usecols=(2))
    dataovito2=np.array(load2)

size2 = len(dataovito1) #total atoms

natoms = size2

xmin = np.min(dataovito1[:,0])
ymin = np.min(dataovito1[:,1])
zmin = np.min(dataovito1[:,2])

outFile = open('CSH-4nm.xyz', 'w')
outFile.write('%i \n' %(natoms))
outFile.write('\n')

#need to transfer the origin to zero, so deduct xmin, ymin and zmin with x, y and z

for j in range(size2): #writing atoms without water
    if dataovito2[j]==1:
        outFile.write('%s %f %f %f \n' %('Sic', dataovito1[j,0]-xmin, dataovito1[j,1]-ymin, dataovito1[j,2]-zmin))
    elif dataovito2[j]==2:
        outFile.write('%s %f %f %f \n' %('Ca', dataovito1[j,0]-xmin, dataovito1[j,1]-ymin, dataovito1[j,2]-zmin))
    elif dataovito2[j]==3:
        outFile.write('%s %f %f %f \n' %('Oc', dataovito1[j,0]-xmin, dataovito1[j,1]-ymin, dataovito1[j,2]-zmin))
    elif dataovito2[j]==4:
        outFile.write('%s %f %f %f \n' %('Hc', dataovito1[j,0]-xmin, dataovito1[j,1]-ymin, dataovito1[j,2]-zmin))
    #elif dataovito2[j]==4:
     #   outFile.write('%s %f %f %f \n' %('O', dataovito1[j,0], dataovito1[j,1], dataovito1[j,2]))

 
outFile.close()
print "All done!"        
