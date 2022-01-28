#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: CQML
"""
import re
import sys
import numpy as np
import os
#########################################################
#   This code make your poscar file having no shift edge! (default fraction=0.1)
#
#How to use it:
#
#   python reWritePOSCAR.py ${inputPOSCAR} ${savePOSCAR} ${defect1}
#
#   or
#  
#   python reWritePOSCAR.py ${inputPOSCAR} ${savePOSCAR} ${defect1} ${defect2} 
#
#########################################################
class AtomData:
    #Read the atom data and lattice parameter data from poscar file
    def __init__(self, input_POSCAR):
        if os.path.isfile(input_POSCAR):
            self.input=input_POSCAR
        else:
            print("There is no input POSCAR file!!")
            sys.exit(0)

        self.latt_constant = 1
        self.atom_data=[]       #name of atom species
        self.atom_data_num=[]   #number of atoms
        self.unit_cell_para=[]  #lattice parameter
        self.unit_cell_info=[]  #atom name & atom fracational coordination
        self.unit_cell_xyz=[]   #only atom fractional coordi
        self.unit_cell_xyz_real=[] #cartesian coordiante of atom position
        
        self.read_POSCAR()
    
    def read_POSCAR(self):
        print("!Input File name : ",self.input)
        with open(self.input, 'r') as f:
            temp_unit_cell_info=[]
        
            for i, line in enumerate(f):
                inputdata=line.strip().split()
                if i ==1:
                    self.latt_constant=list(map(float,inputdata))[0]
                elif i>=2 and i <=4:
                    self.unit_cell_para.append(list(map(float,inputdata)))
                elif i==5:
                    self.atom_data.append(inputdata)
                elif i==6:
                    self.atom_data_num.append(list(map(int,inputdata)))
                elif i>7:
                    unit_list=list(map(float,inputdata))
                    temp_unit_cell_info.append(unit_list)
        f.close()
        
        #consider the lattice constant and self.unit_cell_para 
        if self.latt_constant != 1:
            for i in range(len(self.unit_cell_para)):
                for j in range(3):
                    self.unit_cell_para[i][j]=self.unit_cell_para[i][j]*self.latt_constant

        #remove the 2D list of atom_data & atom_data_num
        self.atom_data=sum(self.atom_data,[])
        self.atom_data_num=sum(self.atom_data_num,[])

        #self.unit_cell_info
        for j in range(0,len(self.atom_data)):
            for i in range(0, self.atom_data_num[j]):
                temp_list=[self.atom_data[j],temp_unit_cell_info[0]]
                del temp_unit_cell_info[0]
                self.unit_cell_info.append(temp_list)

        #self.unit_cell_xyz & self.unit_cell_xyz_real
        for i in range(len(self.unit_cell_info)):
            self.unit_cell_xyz.append(self.unit_cell_info[i][1][:3])
        self.unit_cell_xyz=np.array(self.unit_cell_xyz,dtype=np.float64)
        self.unit_cell_xyz_real=(np.dot(self.unit_cell_para,self.unit_cell_xyz.T)).T
#########################################################


#########################################################
class WritePOSCAR(AtomData):
    #using the AtomData class information, 
    #consider the defect atoms and nearest atoms
    #write the new poscar file
    def __init__(self,inputPOSCAR,savefile):
        self.inputfile=inputPOSCAR
        self.savefile=savefile

        #inheritance of AtomData information
        super().__init__(inputPOSCAR)

        #consider the Defect position using 
        #self.Make_Defect(DefectList)
        #self.Find_Nearest_atom()

        #write the poscar file considering the defect position
        #self.Write_POSCAR_wiDefect()

    def Make_Defect(self,DefectList):
        self.DefectList=DefectList
        ######################
        #consider first defect 
        defect1=self.DefectList[0]
        ref_posi_frac=[0.5,0.5,0.5]
        temp_position=(np.dot(self.unit_cell_para,np.array(ref_posi_frac).T)).T

        #find the defect index
        Index_defect1=self.Find_Nearest_atom(temp_position)[0]
        Posi_defect1=self.unit_cell_xyz_real[Index_defect1]

        #find the nearest index
        Index_Near_defect1=self.Find_Nearest_atom(Posi_defect1)

        self.Posi_index_defect=[Index_defect1]
        Posi_index_Near_defect=Index_Near_defect1

        ######################
        #condisder second defect
        if len(DefectList)==2:
            #find the defect index
            Index_defect2=self.Find_Nearest_atom(Posi_defect1)[0]
            Posi_defect2=self.unit_cell_xyz_real[Index_defect2]

            #find the nearest index
            Index_Near_defect2=self.Find_Nearest_atom(Posi_defect2)

            self.Posi_index_defect.append(Index_defect2)
            Posi_index_Near_defect+=Index_Near_defect2
        
        print("Defect index in POSCAR :",self.Posi_index_defect)
        #final defect and nearest atom index with no overlapping
        self.Posi_index_Near_defect=list(set(Posi_index_Near_defect).difference(set(self.Posi_index_defect)))
        print("Nearest index in POSCAR :", self.Posi_index_Near_defect)

    def ReName_defect_atom(self):
        ######################
        #re-name of defect and nearest atom 
        #self.atom_data & self.atom_data_num & self.unit_cell_info
        #
        #re-nameing of nearest atom (self.Posi_index_Near_defect)
        Nearest_unit_cell_info=[]
        for temp_index in self.Posi_index_Near_defect:
            target_index=self.Find_index_unit_cell_info(temp_index)
            temp_cell_info=self.unit_cell_info[target_index]
            temp_cell_info[0]=temp_cell_info[0].lower()*2

            self.unit_cell_info.pop(target_index)
            Nearest_unit_cell_info.append(temp_cell_info)
            Nearest_unit_cell_info.sort()
        print("Nearest atom infomation :")
        print(Nearest_unit_cell_info)
 

        #re-nameing of defect atom (self.Posi_index_defect)
        Defect_unit_cell_info=[]
        count=0
        for temp_index in self.Posi_index_defect:
            target_index=self.Find_index_unit_cell_info(temp_index)
            temp_cell_info=self.unit_cell_info[target_index]
            
            temp_cell_info[0]=self.DefectList[count] 
            self.unit_cell_info.pop(target_index)
            Defect_unit_cell_info.append(temp_cell_info)
            Defect_unit_cell_info.sort()
            count+=1
        print("Defect atom infomation :")
        print(Defect_unit_cell_info)

        #summary unit_cell_info
        self.unit_cell_info+=Nearest_unit_cell_info
        self.unit_cell_info+=Defect_unit_cell_info

        #update the atom data
        counts=dict()
        for i in range(len(self.unit_cell_info)):
            name = self.unit_cell_info[i][0]
            if name not in counts:
                counts[name] =1
            else:
                counts[name] +=1
        #print(counts)
 
        temp_atom_data=list(counts.keys())
        self.atom_data=temp_atom_data

        temp_atom_data_num=list(counts.values())
        self.atom_data_num=temp_atom_data_num
        #print(temp_atom_data)
        #print(temp_atom_data_num)


    def Find_index_unit_cell_info(self,temp_index):
        for i in range(len(self.unit_cell_info)):
            if (np.array(self.unit_cell_info[i][1][:3]) == self.unit_cell_xyz[temp_index][:3]).all():
                temp=i
        return temp

    
    def Find_Nearest_atom(self,defect):
        #using defect position and self.unit_cell_xyz_real
        #we will find the nearest atoms
        Distances=[]
        for i in range(len(self.unit_cell_xyz_real)):
            temp_dist=self.Check_short_length(defect,self.unit_cell_xyz_real[i])

            Distances.append([i,temp_dist])

        Distances.sort(key=lambda x: (x[1],x[0]))
    
        #no consider the dist==0 atom position
        if Distances[0][1]!=0:
            min_dist=Distances[0][1]
            min_dist_index=Distances[0][0]
        else:
            min_dist=Distances[1][1]
            min_dist_index=Distances[1][0]

        index=[min_dist_index]

        err=0.6 #angstrom
        for i in range(1,len(Distances)):
            if (Distances[i][1]-min_dist)**2 < err:
                index.append(Distances[i][0])
            else:
                break
        #remove the overlap index in index list
        index=list(set(index))

        return index
    
    def Check_short_length(self,X_defect, X_atom ):
        #in cell
        temp_dist=self.Cal_dist(X_defect,X_atom)
        #print(temp_dist)

        #a-axis
        ext_a=(np.dot(self.unit_cell_para,np.array([1,0,0]))).T
        #b-aixs
        ext_b=(np.dot(self.unit_cell_para,np.array([0,1,0]))).T
        #c-aixs
        ext_c=(np.dot(self.unit_cell_para,np.array([0,0,1]))).T

        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    ext_X_atom=X_atom+i*ext_a+j*ext_b+k*ext_c
                    ext_temp_dist=self.Cal_dist(X_defect,ext_X_atom)
                    if temp_dist > ext_temp_dist:
                        temp_dist = ext_temp_dist
        return temp_dist

    def Cal_dist(self,list1,list2):
        dist=((list1[0]-list2[0])**2+(list1[1]-list2[1])**2+(list1[2]-list2[2])**2)**(0.5)
        return dist
     
    def Write_POSCAR_wiDefect(self):
        print('!Save File name is : ',self.savefile)
        with open(self.savefile, 'w') as f:
            f.write("POSCAR\n")
            f.write("1.00000\n")
            for i in range(len(self.unit_cell_para)):
                line="{:>10.7f}    {:>10.7f}    {:>10.7f}\n".format(
                    float(self.unit_cell_para[i][0]),
                    float(self.unit_cell_para[i][1]),
                    float(self.unit_cell_para[i][2]))
                f.write(line)
            
            for i in range(len(self.atom_data)):
                f.write("{} ".format(self.atom_data[i]))
            f.write("\n")

            for i in range(len(self.atom_data_num)):
                f.write("{} ".format(self.atom_data_num[i]))
            f.write("\n")
            
            f.write("Direct\n")
            for i in range(len(self.unit_cell_info)):
                for j in range(len(self.unit_cell_info[i][1])):
                    if j > 2:
                        line=" {:>10.2f}\t".format(self.unit_cell_info[i][1][j])
                    else:
                        line=" {:>10.7f}\t".format(self.unit_cell_info[i][1][j])
                    f.write(line)
                f.write("\n")
        f.close()        

###########################################################

if __name__ == "__main__":

    if len(sys.argv)!=4 and len(sys.argv)!=5:
        print("This code make the no shift zone in POSCAR file!!")
        print("Parameter : [1] input-POSCAR, [2] save-POSCAR name")
        print("            [3] Defect candidate 1, [4] Defect candidate 2")
        sys.exit()
    
    #parameter
    inputPOSCAR=sys.argv[1]
    savePOSCAR=sys.argv[2]
    
    List_defect=[]
    if len(sys.argv)>=4:
        defect1=sys.argv[3] #name of defect candidate 1
        List_defect.append(defect1)
    if len(sys.argv)>=5:
        defect2=sys.argv[4] #name of defect cancdidate 2
        List_defect.append(defect2)
    
    
    #read Atom data
    #process1=AtomData(inputPOSCAR)
    
    #print(process1.unit_cell_info)
    #print(process1.unit_cell_xyz_real)
    
    
    #find defect atom position
    
    
    #sys.exit()
    
    #save new POSCAR file
    process2=WritePOSCAR(inputPOSCAR,savePOSCAR)

    #make defect in new POSCAR file
    process2.Make_Defect(List_defect)
    process2.ReName_defect_atom()


    #finally save the poscar file wi defect position
    process2.Write_POSCAR_wiDefect()
    print()

