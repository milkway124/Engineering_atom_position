#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 22:06:00 2021

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
#   python reWritePOSCAR.py ${inputPOSCAR} ${savePOSCAR}
#
#   or
#  
#   python reWritePOSCAR.py ${inputPOSCAR} ${savePOSCAR}  False True
#
#   or
#
#   python reWritePOSCAR.py ${inputPOSCAR} ${savePOSCAR}  True True
#########################################################
class AtomData:
    def __init__(self, input_POSCAR):
        if os.path.isfile(input_POSCAR):
            self.input=input_POSCAR
        else:
            print("There is no input POSCAR file!!")
            sys.exit(0)

        self.latt_constant = 1
        self.atom_data=[]
        self.atom_data_num=[]
        self.unit_cell_para=[]
        self.unit_cell_info=[]
        
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
        
        if self.latt_constant != 1:
            for i in range(len(self.unit_cell_para)):
                for j in range(3):
                    self.unit_cell_para[i][j]=self.unit_cell_para[i][j]*self.latt_constant

        self.atom_data=sum(self.atom_data,[])
        self.atom_data_num=sum(self.atom_data_num,[])

        for j in range(0,len(self.atom_data)):
            for i in range(0, self.atom_data_num[j]):
                temp_list=[self.atom_data[j],temp_unit_cell_info[0][0],
                                            temp_unit_cell_info[0][1],
                                            temp_unit_cell_info[0][2]]
                del temp_unit_cell_info[0]
                self.unit_cell_info.append(temp_list)

        self.unit_cell_info=np.array(self.unit_cell_info)
        
        self.unit_cell_xyz=self.unit_cell_info[:,1:].astype(float)
        self.unit_cell_xyz_real=(np.dot(self.unit_cell_para,self.unit_cell_xyz.T)).T
#########################################################
    
#########################################################
class WritePOSCAR:
    def __init__(self,savefile,ClassAtomData,fraction=0.1,Along_Z=False,test=False):
        self.savefile=savefile
        self.Read_previous_data(ClassAtomData)

        if not Along_Z:
            self.EdgeFraction=fraction
            print("Edgefraction :",self.EdgeFraction)
            self.Write_FullEdge_POSCAR()
            if test:
                self.Write_FullEdge_POSCAR_NoShift() #for check, No shift zone
        else:
            Height_standard=self.Maxium_Height()
            err=0.2 #angstrom
            print("Standard Height : ",Height_standard,"+-",err, "[Angstrom]")
            self.Write_Along_Z_POSCAR(Height_standard,err)
            if test:
                self.Write_Along_Z_POSCAR_OnlyShift(Height_standard,err) #for check, only shift zone

    def Read_previous_data(self,ClassAtomData):
        self.unit_cell_para=ClassAtomData.unit_cell_para
        self.atom_data=ClassAtomData.atom_data
        self.atom_data_num=ClassAtomData.atom_data_num
        self.unit_cell_info=ClassAtomData.unit_cell_info
        self.unit_cell_xyz=ClassAtomData.unit_cell_xyz
        self.unit_cell_xyz_real=ClassAtomData.unit_cell_xyz_real

    def EdgeSide(self, unit_cell_xyz):
        LowLimit=self.EdgeFraction
        HighLimit=1-self.EdgeFraction
        for i in range(len(unit_cell_xyz)):
            if( unit_cell_xyz[i] <= LowLimit and unit_cell_xyz[i] >= -LowLimit) or (unit_cell_xyz[i] >= HighLimit and unit_cell_xyz[i] <= 2-HighLimit):
                return True

        return False

    def Check_Standard(self,standard,unit_cell_xyz_real,err):
        if (unit_cell_xyz_real[2] <= standard+err) and (unit_cell_xyz_real[2] >= standard-err):
            return True
        return False

    def Maxium_Height(self):
        temp_Z=0
        for i in range(len(self.unit_cell_xyz_real)):
            if (self.unit_cell_xyz_real[i][2] > temp_Z):
                temp_Z=self.unit_cell_xyz_real[i][2]
        return temp_Z

        
    def Write_FullEdge_POSCAR(self):
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
            for i in range(len(self.unit_cell_xyz)):
                if self.EdgeSide(self.unit_cell_xyz[i]):
                    line="{:>10.9f}    {:>10.9f}    {:>10.9f}  0  0  0\n".format(
                        self.unit_cell_xyz[i][0],
                        self.unit_cell_xyz[i][1],
                        self.unit_cell_xyz[i][2])
                else:
                    line="{:>10.9f}    {:>10.9f}    {:>10.9f}\n".format(
                        self.unit_cell_xyz[i][0],
                        self.unit_cell_xyz[i][1],
                        self.unit_cell_xyz[i][2])
                f.write(line)
            
        f.close()        

    def Write_FullEdge_POSCAR_NoShift(self):
        print('!Save File name(NoShift) is : ',self.savefile+'_test')
        ListEdge=[]
        for i in range(len(self.unit_cell_xyz)):
            if self.EdgeSide(self.unit_cell_xyz[i]):
                ListEdge.append(i)

        new_atom_data=[]
        new_atom_data_num=[]
        for i in ListEdge:
            new_atom_data.append(self.unit_cell_info[i][0])
        count={}
        for i in new_atom_data:
            try: count[i] += 1
            except: count[i] =1
        new_atom_data=list(count.keys())
        new_atom_data_num=list(count.values())

        with open(self.savefile+'_test', 'w') as f:
            f.write("POSCAR\n")
            f.write("1.00000\n")
            for i in range(len(self.unit_cell_para)):
                line="{:>10.7f}    {:>10.7f}    {:>10.7f}\n".format(
                    float(self.unit_cell_para[i][0]),
                    float(self.unit_cell_para[i][1]),
                    float(self.unit_cell_para[i][2]))
                f.write(line)
            
            for i in range(len(new_atom_data)):
                f.write("{} ".format(new_atom_data[i]))
            f.write("\n")

            for i in range(len(new_atom_data_num)):
                f.write("{} ".format(new_atom_data_num[i]))
            f.write("\n")
            
            f.write("Direct\n")
            for i in ListEdge:
                line="{:>10.9f}    {:>10.9f}    {:>10.9f}  0  0  0\n".format(
                    self.unit_cell_xyz[i][0],
                    self.unit_cell_xyz[i][1],
                    self.unit_cell_xyz[i][2])
                f.write(line)
            
        f.close()        

    def Write_Along_Z_POSCAR(self,Height_standard,err):
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
            for i in range(len(self.unit_cell_xyz)):
                if self.Check_Standard(Height_standard,self.unit_cell_xyz_real[i],err):
                    line="{:>10.9f}    {:>10.9f}    {:>10.9f}\n".format(
                        self.unit_cell_xyz[i][0],
                        self.unit_cell_xyz[i][1],
                        self.unit_cell_xyz[i][2])
                else:
                    line="{:>10.9f}    {:>10.9f}    {:>10.9f}  0  0  0 \n".format(
                        self.unit_cell_xyz[i][0],
                        self.unit_cell_xyz[i][1],
                        self.unit_cell_xyz[i][2])
                f.write(line)
            
        f.close()        

    def Write_Along_Z_POSCAR_OnlyShift(self,Height_standard,err):
        print('!Save File name(OnlyShift) is : ',self.savefile+'_test')
        ListShift=[]
        for i in range(len(self.unit_cell_xyz)):
            if self.Check_Standard(Height_standard,self.unit_cell_xyz_real[i],err):
                ListShift.append(i)
        new_atom_data=[]
        new_atom_data_num=[]
        for i in ListShift:
            new_atom_data.append(self.unit_cell_info[i][0])
        count={}
        for i in new_atom_data:
            try: count[i] += 1
            except: count[i] = 1
        new_atom_data=list(count.keys())
        new_atom_data_num=list(count.values())


        with open(self.savefile+'_test', 'w') as f:
            f.write("POSCAR\n")
            f.write("1.00000\n")
            for i in range(len(self.unit_cell_para)):
                line="{:>10.7f}    {:>10.7f}    {:>10.7f}\n".format(
                    float(self.unit_cell_para[i][0]),
                    float(self.unit_cell_para[i][1]),
                    float(self.unit_cell_para[i][2]))
                f.write(line)
            
            for i in range(len(new_atom_data)):
                f.write("{} ".format(new_atom_data[i]))
            f.write("\n")

            for i in range(len(new_atom_data_num)):
                f.write("{} ".format(new_atom_data_num[i]))
            f.write("\n")
            
            f.write("Direct\n")
            for i in ListShift:
                line="{:>10.9f}    {:>10.9f}    {:>10.9f}\n".format(
                    self.unit_cell_xyz[i][0],
                    self.unit_cell_xyz[i][1],
                    self.unit_cell_xyz[i][2])
                f.write(line)
        f.close()        
###########################################################

if len(sys.argv)==1:
    print("This code make the no shift zone in POSCAR file!!")
    print("Parameter : [1] input-POSCAR, [2] save-POSCAR name, [3] Along Z opt (True or False)")
    print("            [4] option for print no-shift or only shift atom POSCAR (True or False)")
    print("            Default opt [3], [4] : False")
    sys.exit()

#parameter
inputPOSCAR=sys.argv[1]
savePOSCAR=sys.argv[2]
Along_Z=False
test=False

if len(sys.argv)==4:
    Along_Z=sys.argv[3] #consider z-aixs plane
    test=False

if len(sys.argv)==5:
    Along_Z=sys.argv[3] #consider z-aixs plane
    test=sys.argv[4]    #see the no shift zone or only shift zone
#print(sys.argv)

#insert Atom data
process1=AtomData(inputPOSCAR)

#save new POSCAR file
process2=WritePOSCAR(savePOSCAR,process1,Along_Z=Along_Z,test=test)



