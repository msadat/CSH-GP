import numpy as np
import glob
import random

def delSiO2 (SiID, Otype):
    Si_list = [] #list of Si atoms (list of [lists of different Si atoms and their neighbors])
    Otype_full = []
    Obridging = [] #list of bridging O with two Si or Al neighbors (i.e. Si-O-Al)
    ONBO = []

    #***************Make a list of Si atoms and their neighbors***********************************
    
    with open(SiID,  'r') as f:
     for line in f:
                
        words_Si= line.split() # list of each lines
        
        for x in range(len(words_Si)):
            words_Si[x] = int(words_Si[x])
        Si_list.append(words_Si)
    #print len(Si_list)
    
    #***************Make a list of all oxygens**************************
    with open(Otype,  'r') as f:    
      for line in f:
                
        words_Otype = line.split() # list of each lines
      
        for x in range(len(words_Otype)):
            words_Otype[x] = int(words_Otype[x]) 
        Otype_full.append(words_Otype)
    #***************Make a list of bridging oxygens**************************

    for i in range(len(Otype_full)):
        bridgeO = 0
        
        for j in range(2,len(Otype_full[i])):
            if Otype_full[i][j] == 1: 
              bridgeO +=1
	    else:
	      bridgeO = bridgeO
        if bridgeO == 2:  
            Obridging.append(Otype_full[i][0])
            ONBO = ONBO
        elif bridgeO ==1:
            Obridging = Obridging  
            ONBO.append(Otype_full[i][0])
        else:
            ONBO = ONBO
            Obridging = Obridging
    print "%NBO= " , float(len(ONBO))/float(len(Otype_full))	
    #**********************Make a list of Si atoms with NBOs*****************

    Si_NBO = [] #list of Si with NBOs
    NBO_Si = [] #of list of corresponding NBOs
    SiNBO_combined_list = [] #combined list of [concateneted list of Si and their NBOs]
    
    for i in range(len(Si_list)):
         SiNBO_list = []
         for j in range(2,len(Si_list[i])):
	    #SiNBO_list = []
            if Si_list[i][j] in ONBO: #check if the O neighbor is in the ONBO list
                #Si_NBO.append(Si_list[i][0])
                #NBO_Si.append(Si_list[i][j])
                #SiNBO_combined_list.append(Si_NBO+NBO_Si) 
		SiNBO_list.append(Si_list[i][0])
 		SiNBO_list.append(Si_list[i][j])

            else:
                #Si_NBO = Si_NBO
                #NBO_Si = NBO_Si
    
                SiNBO_list = SiNBO_list
         SiNBO_combined_list.append(SiNBO_list)
    
    print "No. of Si with NBO =" , len(SiNBO_combined_list)
    #print SiNBO_combined_list
    #Si_random = np.random.random_integers(1,len(SiNBO_combined_list)) #randomly select a list of Si and its NBO neighbors
    #print SiNBO_combined_list[Si_random]
    return SiNBO_combined_list # [Si_random]
    
def delAtoms(x,t):
    """x = array of xyz data file with the atom types,
    t = list of atoms to be deleted """
    #np.delete won't work as the new list after deletion has changed row positions
    idx = []
    #print "TTT", len(x), len(t)
    for i in range(len(t)):
      for j in range(len(x)):
        if t[i] in x[j]:
          idx.append(j)
	  #x = np.delete(x, x[j],0)  #'0' is the axis notation for row, '1' is col
        else:
          idx = idx
    return np.delete(x,idx,0) 

if __name__ == "__main__":

    filename = 'data.Tobermorite-ortho_serial' #original lammps datafile

    flist = glob.glob(filename)
    natoms_GBP = 4773
    target_CS = 1.8
    for f in flist:
        #load = np.genfromtxt(f, dtype=float, skip_header=19, skip_footer=natoms_GBP+1, usecols=(0,2,4,5,6))
        load = np.genfromtxt(f, dtype=float, skip_header=18, usecols=(0,2,4,5,6))
	data=np.array(load)
    
    nSi = 0
    nCa = 0
    
    for j in range(len(data)):
        if data[j,1] == 1:
            nSi +=1
        elif data[j,1] == 2:
            nCa +=1  
    #print nSi, nCa     
    CS = float(nCa)/float(nSi)        
    print "Initial C/S = " , CS
    print "Initial No. of Si =", nSi
    print "Initial No. of Ca =", nCa

    dataSiNBO = delSiO2('Nlist-ID_Si_Al-Tobermorite', 'Nlist-type_O-Tobermorite') #list of Si and their NBOs
    #print len(dataSiNBO)
    #print len(dataSiNBO[11])
    #data = delAtoms(data, dataSiNBO[np.random.random_integers(0,len(dataSiNBO)-1)])
    #data = delAtoms(data, dataSiNBO[11]) #deletes ith row from the data

    print len(data)
    natoms = len(data)

    	
    while CS < target_CS:
      nSi = 0
      nCa = 0             
      #data = delAtoms(data, dataSiNBO[np.random.random_integers(0,len(dataSiNBO)-1)])
      data = delAtoms(data, dataSiNBO[random.sample(xrange(0,len(dataSiNBO)-1),1)[0]])

      #print len(data)
      for j in range(len(data)):
            if data[j,1] == 1:
                nSi +=1
            elif data[j,1] == 2:
                nCa +=1       
      CS = float(nCa)/float(nSi)  
      if CS==0:
	break
      print "Current C/S= ", CS
      print "Current No. of Si atoms = ", nSi
      print "Current No. of Ca atoms = ", nCa

      print "Getting a new C/S..."
      
      if CS >= target_CS:
         print "Desired C/S ratio reached"
         break
      
    print "Finally obtained C/S is = ", CS
    
    natoms = len(data)
    print "Finally remaining atoms = ", natoms

    

    #****************WRITE final data TO LAMMPS DATA FILE**********************************

    outfilename = 'Tobermorite_CS_'+str(target_CS)
     
    outFile = open('data.'+outfilename, 'w')
    outFile.write('LAMMPS data file written using Python script\n')
    outFile.write('\n')
    outFile.write('%i %s \n' %(natoms, 'atoms'))
    outFile.write('4 atom types \n')
    outFile.write('\n')
    outFile.write('%f %f %s %s \n' %(0.00, 42.0, 'xlo', 'xhi'))
    outFile.write('%f %f %s %s \n' %(0.00, 42.0, 'ylo', 'yhi'))
    outFile.write('%f %f %s %s \n' %(0.00, 42.0, 'zlo', 'zhi'))
    outFile.write('\n')
    outFile.write('%s \n' %('Masses'))
    outFile.write('\n')
    outFile.write('%i %f \n' %(1, 28.065)) #Si
    outFile.write('%i %f \n' %(2, 40.078)) #Ca
    outFile.write('%i %f \n' %(3, 1.0)) #H
    outFile.write('%i %f \n' %(4, 16.0)) #O
    outFile.write('\n')
    outFile.write('Atoms\n')
    outFile.write('\n')


    for j in range(len(data)):
        if data[j,1]==1:
            outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 1, 0, data[j,2], data[j,3], data[j,4]))
        elif data[j,1]==2:
            outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 2, 0, data[j,2], data[j,3], data[j,4]))
        elif data[j,1]==3:
            outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 3, 0, data[j,2], data[j,3], data[j,4]))
        elif data[j,1]==4:
            outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 4, 0, data[j,2], data[j,3], data[j,4]))

    outFile.close()
    print "All done!"         
