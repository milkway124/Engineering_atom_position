#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wes Jan 19 02:196:00 2022

@author: CQML
"""
import re
import sys
import numpy as np
import os
from collections import Counter
#########################################################
#   This code make stacking new poscar file using input poscar file!
#
#How to use it:
#
#   python reWritePOSCAR.py  ${savePOSCAR} ${inputPOSCAR} ${fractional_posi}
#
#   or
#  
#   python reWritePOSCAR.py  ${savePOSCAR} ${inputPOSCAR2} ${fractional_posi1} ${inputPOSCAR2} ${fractional_posi2} 
#
#########################################################
class AtomData:
    def __init__(self, input_POSCAR,fractional_posi):
        if os.path.isfile(input_POSCAR):
            self.input=input_POSCAR
        else:
            print("There is no input POSCAR file!!")
            sys.exit(0)

        self.latt_constant = 1
        self.unit_cell_para=[]
        self.atom_data=[]
        self.atom_data_num=[]
        self.atom_line=[]
 
        self.read_POSCAR()
    
        self.atom_dict=dict(zip(self.atom_data,self.atom_data_num))
    
    def read_POSCAR(self):
        print("!Input File name : ",self.input)
        with open(self.input, 'r') as f:
            temp_line=[]
        
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
                    temp_line.append(list(map(float,line.split())))
        f.close()
        
        if self.latt_constant != 1:
            for i in range(len(self.unit_cell_para)):
                for j in range(3):
                    self.unit_cell_para[i][j]=self.unit_cell_para[i][j]*self.latt_constant

        self.atom_data=sum(self.atom_data,[])
        self.atom_data_num=sum(self.atom_data_num,[])

        tmp=0
        for j in range(0,len(self.atom_data)):
            for i in range(0, self.atom_data_num[j]):
                self.atom_line.append([self.atom_data[j],temp_line[tmp]])
                tmp+=1
    
#########################################################

#########################################################
class TotalAtomData:
    def __init__(self):
        self.unit_cell_para=[]
        self.atom_line=[]
        self.atom_dict={}
        self.atom_count=[]

    def add_atom_dict(self,atom_dict):
        self.atom_dict=dict(Counter(self.atom_dict)+Counter(atom_dict))

    def add_atom_line(self,atom_line):
        self.atom_line.extend(atom_line)
        self.atom_count.append(len(atom_line))

    def add_unit_cell_para(self,unit_cell_para):
        self.Check_same_unit_cell_para(unit_cell_para)

    def Check_same_unit_cell_para(self,unit_cell_para):
        if len(self.unit_cell_para) == 0:
            self.unit_cell_para.extend(unit_cell_para)
        else:
            for j in range(3):
                for i in range(3):
                    if (self.unit_cell_para[i][j]) != (unit_cell_para[i][j]):
                        print("The unit cell parameter[%d][%d] is not Same!!!"%(i,j))
                        print("unit cell parameter[%d][%d] : %lf & %lf !!"%(i,j,self.unit_cell_para[i][j],unit_cell_para[i][j]))
                        sys.exit(0)
    
    def rePositions_Z(self,new_Z,frational_list):
        if len(list(set(frational_list)))!=len(fractional_list):
            print("The fractional position is wrong!! Maybe overlap!!")
            print(" ---> ",fractional_list)
            sys.exit()

        #re-write atomic positioin
        line=0
        new_atom_Z=[]
        for i in range(len(self.atom_count)):
            min_atomic_posi=1
            fraction_line=[line,line+self.atom_count[i]]
            #find minimum atomic posi
            for j in range(fraction_line[0],fraction_line[1]):
                temp=self.atom_line[j][1][2]
                if min_atomic_posi > temp:
                    min_atomic_posi=temp
                line+=1
            #shift atomic position and re-size Z constant
            min_atom_z=self.unit_cell_para[2][2]*min_atomic_posi
            for j in range(fraction_line[0],fraction_line[1]):
                new_atom_Z.append(float(fractional_list[i])+ (self.unit_cell_para[2][2]*self.atom_line[j][1][2]-min_atom_z)/(new_Z))

        for j in range(len(new_atom_Z)):
            self.atom_line[j][1][2]=new_atom_Z[j]    

        #sort atomic position along atom_dict's key
        self.unit_cell_para[2][2]=new_Z

        #re-write unit cell parameter         
        new_atom_line=[]
        list_atom_dict=list(self.atom_dict.keys())
        for i in range(len(list_atom_dict)):
            for j in range(len(self.atom_line)):
                if self.atom_line[j][0] == list_atom_dict[i]:
                    new_atom_line.append(self.atom_line[j])
        self.atom_line=new_atom_line
            

 
#########################################################

#########################################################
class WritePOSCAR:
    def __init__(self,savefile,unit_cell_para,atom_dict,atom_line):
        self.savefile=savefile
        
        self.Write_POSCAR(unit_cell_para,atom_dict,atom_line)

    def Check_in_unitCell(self,posi):
        if posi >=0 and posi < 1:
            return posi

        if posi >= 1:
            posi=posi-1
        elif posi <0:
            posi=posi+1
        posi=self.Check_in_unitCell(posi)
        return posi
        
    def Write_POSCAR(self,unit_cell_para,atom_dict,atom_line):
        print('!Save File name is : ',self.savefile)
        with open(self.savefile, 'w') as f:
            f.write("POSCAR\n")
            f.write("1.00000\n")
            for i in range(len(unit_cell_para)):
                line="{:>10.9f}    {:>10.9f}    {:>10.9f}\n".format(
                    float(unit_cell_para[i][0]),
                    float(unit_cell_para[i][1]),
                    float(unit_cell_para[i][2]))
                f.write(line)
            
            for i in range(len(atom_dict)):
                f.write("{}\t".format(list(atom_dict.keys())[i]))
            f.write("\n")

            for i in range(len(atom_dict)):
                f.write("{}\t".format(list(atom_dict.values())[i]))
            f.write("\n")
            
            f.write("Direct\n")
            for i in range(len(atom_line)):
                for j in range(len(atom_line[i][1])):
                    if j < 3:
                        line="{:>10.9f} \t ".format(self.Check_in_unitCell(atom_line[i][1][j])) 
                        #line="{:>10.9f} \t ".format(atom_line[i][1][j]) 
                        f.write(line)
                    else:
                        line="{} ".format(atom_line[i][1][j])
                        f.write(line)
                line="\n"
                f.write(line)

        f.close()        

###########################################################

if len(sys.argv)==1:
    print("This code make stacking new poscar file using input poscar file!!")
    print("Parameter : [1] save-POSCAR name, [2] input-POSCAR, [3] fractional position along minimal atomic position")
    print("            if you want to using more input-poscar...")
    print("            [4] input-POSCAR2 [5] fractional position 2 ")
    sys.exit()

if len(sys.argv)%2 != 0:
    print("This code have to input the below parameters!!")
    print("Parameter : [1] save-POSCAR name, [2] input-POSCAR, [3] fractional position")
    sys.exit()

#TotalAtomdata
Total=TotalAtomData()
new_Z=6.65
print("\nNew Z lattice :",new_Z)


#parameter
savePOSCAR=sys.argv[1]
fractional_list=[]
for i in range(int((len(sys.argv)-2)/2)):
    inputPOSCAR=sys.argv[(i+1)*2]

    fractional_posi=sys.argv[(i+1)*2+1]
    fractional_list.append(fractional_posi)

    #insert Atom data
    process1=AtomData(inputPOSCAR,fractional_posi)

    Total.add_atom_dict(process1.atom_dict)
    Total.add_atom_line(process1.atom_line)
    Total.add_unit_cell_para(process1.unit_cell_para)


#re position using frational-position
Total.rePositions_Z(new_Z,fractional_list)

print("\nTotal Class Check :")
print(Total.unit_cell_para)
print(Total.atom_dict)
#print((Total.atom_line))
print()

#save new stacking POSCAR file
process2=WritePOSCAR(savePOSCAR,Total.unit_cell_para,Total.atom_dict,Total.atom_line)



