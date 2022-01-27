# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 15:39:31 2021

@author: CQML
"""

import numpy as np
from numpy import linalg as LA
import sys
import os
import NamePoint2Bar as NPB

########################################################################
#this code make you a curvature structure (gaussian shape only : x-plane & y-plane & xy-plane)
#How to use it:
#   import BubbleMaker_module as BMm
#
##  read the input data
#   A=BMm.AtomData(inputfile)
#
##  modify the bubble shape
#   B=BMm.MakeBubble(A,hightest_Z,FWHM)
#   B.MakeBubbleStructure()
#   B.SearchNoShiftArea(limit_err)
#
##   write the bubble shape
#   C=BMm.WritePOSCAR(B,output)
#   C.Write_NoShift_POSCAR(outputnoShift)
#   C.Write_OnlyShift_POSCAR(outputOnlyShift)
#
########################################################################
#   First. Class AtomData
#   1.input      : only POSCAR file  // output format is POSCAR
########################################################################
#   Second, Class MakeBubble
#   1.self.hightest_move    : The Bubble hight , unit : Angstrom
#   2.self.FWHM             : FWHM of Gaussian function , unit : Angstrom
#   3.self.origin           : the center(x,y) of Gaussian function, atomic position-->[0~1],[0~1]
#   4.self.limit_err        : How long atoms move in selected area to 
#                             not move atom by the curvature , unit : angstrom //default : 0.01 angstrom
########################################################################
#   Third, Class WritePOSCAR
#   1.output    : POSCAR file name // format is POSCARfile
########################################################################

# Read atom data in POSCAR file
class AtomData:
    
    def __init__(self,input_POSCAR):
        if os.path.isfile(input_POSCAR):
            self.input=input_POSCAR
        else:
            print("There is no input POSCAR file!!")
            sys.exit(0)
        
        
        #set condition
        self.latt_constant = 1
        
        self.atom_data=[]           #atom specise
        self.atom_data_num=[]       #atom num
        self.unit_cell_para=[]      #unit-cell a,b,c
        self.unit_cell_info=[]      #[atom Name, x, y, z] ; x & y & z is fraction position

        self.read_POSCAR()
        
        # self.unit_cell_xyz        #[x,y,z] ; x & y & z is fraction position
        # self.unit_cell_xyz_real   #[x_Real,y_Real,z_Real] ; x_Real & y_Real & z_Real is real position
    
    
    def read_POSCAR(self):
        print('!Input File name is : ',self.input)
        with open(self.input, 'r') as f:
            temp_unit_cell_info=[]
            
            for i, line in enumerate(f):
                inputdata=line.strip().split()
                if i ==1:
                    self.latt_constant=list(map(float,inputdata))[0]
                elif i>=2 and i<=4:
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
                temp_list=[self.atom_data[j], temp_unit_cell_info[0][0],
                            temp_unit_cell_info[0][1],
                            temp_unit_cell_info[0][2]]
                del temp_unit_cell_info[0]
                self.unit_cell_info.append(temp_list)
                
        self.unit_cell_info=np.array(self.unit_cell_info)            
        
        self.unit_cell_xyz=self.unit_cell_info[:,1:].astype(np.float64)
        self.unit_cell_xyz_real=(np.dot(self.unit_cell_para,self.unit_cell_xyz.T)).T
        
########################################################################


# main code to make bubble structure
class MakeBubble:
    ####################################################################
    def __init__(self, AtomData,hightest_move,FWHM ):
        #Set previous AtomData         
        self.Read_previou_data(AtomData)

        #Conditions#
        #1. atom is not move along z-axis : X-range [0~(move_limit_range)%] & [(1-move_limit_range)%~100%]
        
        
        #2. Bubble shape condition
        self.Bubble_hight=float(hightest_move) #unit : angstrom
        self.Bubble_FWHM=float(FWHM) #unit : angstrom
        
        #3. set the zero point (origin in gaussian function)
        self.origin=None
        
        print('!atom_data',self.atom_data)
        print('!atom_data_num',self.atom_data_num)
        print('!Read the unit_cell_para',self.unit_cell_para)
        self.min_length_before=self.Find_minlength(self.unit_cell_info[:,1:],self.unit_cell_para)
        print()

        self.limit_range_atom=[]
        self.Distored_unit_cell_info=[]

    ####################################################################        
    def Read_previou_data(self,AtomData):
        self.atom_data=AtomData.atom_data
        self.atom_data_num=AtomData.atom_data_num
        self.unit_cell_para=AtomData.unit_cell_para
        self.unit_cell_info=AtomData.unit_cell_info
        
    def SetOrigin(self,origin=None):
        if origin is not None:
            if len(origin)==2:
                self.origin=[[origin[0],origin[1],0]]
            elif len(origin)==3:
                self.origin=[[origin[0],origin[1],origin[2]]]
            else:
                print("The length of origin should be than 1")
                sys.exit(1)
        else:
            self.FindOrigin()
            
    def FindOrigin(self):
        temp_origin=[0.5,0.5]
        distance=100
        for i in range(len(self.unit_cell_info)):
            tempX=float(self.unit_cell_info[i][1])-temp_origin[0]
            tempY=float(self.unit_cell_info[i][2])-temp_origin[1]
            temp_distance=(tempX**2+tempY**2)**(0.5)
            if temp_distance <= distance:
                distance = temp_distance
                j=i
            
        self.origin=[[float(self.unit_cell_info[j][1]),float(self.unit_cell_info[j][2]),0]]
       
    ####################################################################
    def ShrinkCellPara(self,shrink_streng):
        #strained structure
        print("!Before lattice parameter : ",self.unit_cell_para)
            
        print("!shrink the X & Y-axis lattice parameter ({}%)".format(shrink_streng))
        for i in range(3):
            self.unit_cell_para[0][i]=self.unit_cell_para[0][i]*(1-shrink_streng/100)
            self.unit_cell_para[1][i]=self.unit_cell_para[1][i]*(1-shrink_streng/100)
        
        print("!After lattice parameter : ",self.unit_cell_para)
        
    ####################################################################
    def MakeBubbleStructure(self):
        # set the origin to cal gaussian
        if self.origin is None:
            self.SetOrigin()
        self.RealOrigin=self.AtomPosi2RealPosi(self.unit_cell_para,self.origin)
        self.RealOrigin=np.array(self.RealOrigin)
        self.RealOrigin=np.squeeze(self.RealOrigin,axis=0)
        
        #transform the data type [1]~[3] from str to float
        self.Distored_unit_cell_info=self.unit_cell_info[:,1:].astype(np.float) 
        #atom position -> real position
        Realatomposi=self.AtomPosi2RealPosi(self.unit_cell_para,self.Distored_unit_cell_info)
        Realatomposi=np.array(Realatomposi)
        
        #bubble distored(in view of real-atom position)
        self.sigma=self.FindSigma(self.Bubble_FWHM)
        
        high_z=self.gaussian_2D(Realatomposi[:,0],Realatomposi[:,1],self.Bubble_hight,self.sigma,self.RealOrigin)
        
        Realatomposi[:,2]=Realatomposi[:,2]+high_z
        
        #real position -> atom position
        self.Distored_unit_cell_info=self.RealPosi2AtomPosi(self.unit_cell_para,Realatomposi)
        self.Distored_unit_cell_info=np.array(self.Distored_unit_cell_info)


    def MakeConeStructure(self):
        # set the origin to cal gaussian
        if self.origin is None:
            self.SetOrigin()
        self.RealOrigin=self.AtomPosi2RealPosi(self.unit_cell_para,self.origin)
        self.RealOrigin=np.array(self.RealOrigin)
        self.RealOrigin=np.squeeze(self.RealOrigin,axis=0)
        
        #transform the data type [1]~[3] from str to float
        self.Distored_unit_cell_info=self.unit_cell_info[:,1:].astype(np.float) 
        #atom position -> real position
        Realatomposi=self.AtomPosi2RealPosi(self.unit_cell_para,self.Distored_unit_cell_info)
        Realatomposi=np.array(Realatomposi)
        
        #bubble distored(in view of real-atom position)
        # self.sigma=self.FindSigma(self.Bubble_FWHM)
        
        high_z=self.Corn_2D(Realatomposi[:,0],Realatomposi[:,1],self.Bubble_hight,self.Bubble_FWHM,self.RealOrigin)
        Realatomposi[:,2]=Realatomposi[:,2]+high_z
        
        #real position -> atom position
        self.Distored_unit_cell_info=self.RealPosi2AtomPosi(self.unit_cell_para,Realatomposi)
        self.Distored_unit_cell_info=np.array(self.Distored_unit_cell_info)

    def MakePillarStructure(self,slop=1):
        # set the origin to cal gaussian
        if self.origin is None:
            self.SetOrigin()
        self.RealOrigin=self.AtomPosi2RealPosi(self.unit_cell_para,self.origin)
        self.RealOrigin=np.array(self.RealOrigin)
        self.RealOrigin=np.squeeze(self.RealOrigin,axis=0)
        
        #transform the data type [1]~[3] from str to float
        self.Distored_unit_cell_info=self.unit_cell_info[:,1:].astype(np.float64) 
        #atom position -> real position
        Realatomposi=self.AtomPosi2RealPosi(self.unit_cell_para,self.Distored_unit_cell_info)
        Realatomposi=np.array(Realatomposi)
        
        #bubble distored(in view of real-atom position)
        # self.sigma=self.FindSigma(self.Bubble_FWHM)
        
        high_z=self.Pillar_2D(Realatomposi[:,0],Realatomposi[:,1],self.Bubble_hight,self.Bubble_FWHM,self.RealOrigin,slop)
        Realatomposi[:,2]=Realatomposi[:,2]+high_z
        
        #real position -> atom position
        self.Distored_unit_cell_info=self.RealPosi2AtomPosi(self.unit_cell_para,Realatomposi)
        self.Distored_unit_cell_info=np.array(self.Distored_unit_cell_info)


    def SearchNoShiftArea(self,limit_err=None):
        if len(self.Distored_unit_cell_info)==0:
            print("Must Run MakeBubbleStructure function!!")
            sys.exit(1)
        
        #we want to the special atom which is no shift by gaussian function
        #So, we use the limit_er, this value means the atom is 
        if limit_err is None:
            limit_err=0.01 #unit : angstrom

        
        Before_Distored_unit_cell_info=self.unit_cell_info[:,1:].astype(np.float64) 
        Before_Distored_unit_cell_info=np.array(Before_Distored_unit_cell_info)
        
        After_Distored_unit_cell_info=self.Distored_unit_cell_info
        
        if len(Before_Distored_unit_cell_info) != len(After_Distored_unit_cell_info):
            print("There are some prooblems in MakeBubbleStrcuture function!! Need to Check code!!")
            sys.exit(1)

        #check the no shift atom num and Save num in self.limit_range_atom        
        self.limit_range_atom=[]
        for i in range(len(After_Distored_unit_cell_info)):
            Before=self.AtomPosi2RealPosi(self.unit_cell_para,[Before_Distored_unit_cell_info[i]])
            After=self.AtomPosi2RealPosi(self.unit_cell_para,[After_Distored_unit_cell_info[i]])
            
            Before=np.array(Before)            
            Before=np.squeeze(Before,axis=0)
            After=np.array(After)            
            After=np.squeeze(After,axis=0)

            temp_shift=After[2]-Before[2]
            if temp_shift <= limit_err:
                self.limit_range_atom.append(i)

        # print(self.limit_range_atom)
        

    ####################################################################         
    def AtomPosi2RealPosi(self,unit_cell_para,AtomPosi):
        numpy_unit_cell_para=np.array(unit_cell_para)
        RealPosi=[]
        for i in range(len(AtomPosi)):
            temp=np.dot(np.array(AtomPosi[i]),numpy_unit_cell_para)
            RealPosi.append(temp)
        return RealPosi

    def RealPosi2AtomPosi(self,unit_cell_para,RealPosi):
        numpy_unit_cell_para=np.array(unit_cell_para)
        numpy_unit_inv=np.linalg.inv(numpy_unit_cell_para)
        AtomPosi=[]
        for i in range(len(RealPosi)):
            temp=np.dot(np.array(RealPosi[i]),numpy_unit_inv)
            AtomPosi.append(temp)
        return AtomPosi

    ####################################################################
    def gaussian(self,A,B,x):
        return np.exp(-A*((x-0.5)**2))/B
    
    def gaussian_2D(self,x,y,A,sigma,origin):
        X_term=((x-origin[0])**2)/(2*sigma**2)
        Y_term=((y-origin[1])**2)/(2*sigma**2)
        return A* np.exp( -(X_term + Y_term) )
    
    def gaussian_2D_X_Y(self,x,y,A,sigma_x,sigma_y,origin):
        X_term=((x-origin[0])**2)/(2*sigma_x**2)
        Y_term=((y-origin[1])**2)/(2*sigma_y**2)
        return A* np.exp( -(X_term + Y_term) )
    
    def FindSigma(self,FWHM):
        return FWHM/(2*(2*np.log(2))**(0.5))
    
    def Corn_2D(self,x,y,A,FWHM,origin):
        #line shape
        X_term=(x-origin[0])
        Y_term=(y-origin[1])
        Shift_Z=-A/(2*FWHM/2)*(X_term**2+Y_term**2)**(0.5)+A
        
        #if mius --> no shift
        if len(Shift_Z) == 1:
            if Shift_Z <= 0:
                Shift_Z=0
        else :
            for i in range(len(Shift_Z)):
                if Shift_Z[i] <= 0:
                    Shift_Z[i]=0
                    
        return Shift_Z
    
    def Pillar_2D(self,x,y,hight,Cutoff,origin,slop):
        #using fermi-dirac distribution function
        #exponential distribution
        X_term=(x-origin[0])
        Y_term=(y-origin[1])
        Radius=(X_term**2+Y_term**2)**(0.5)
        
        
        Shift_Z=hight/(1+np.exp((Radius-Cutoff)/slop))
        
        return Shift_Z
        
        
    
    ####################################################################
    def CheckMinlength(self):
        self.min_length=self.Find_minlength(self.Distored_unit_cell_info,self.unit_cell_para)
        print('\n!Minium length between atoms(Before) :',self.min_length_before,'Angstrom')
        print('!Minium length between atoms(After)  :',self.min_length,'Angstrom')
    
    def Find_minlength(self,Distored_unit_cell_info,unit_cell_para):

        unit_information=np.array(Distored_unit_cell_info,dtype='f')
        unit_cell_para=np.array(unit_cell_para)

        New_curved_unit_cell_information=[]
        real_atom=[]
        for i in range(len(Distored_unit_cell_info)):
            for j in range(3):
                real_atom=unit_information[i]*unit_cell_para
            New_curved_unit_cell_information.append(real_atom)
   
        length=[]
        for i in range(len(New_curved_unit_cell_information)-1):
            for row in range(i+1,len(New_curved_unit_cell_information)):
                temp=LA.norm(New_curved_unit_cell_information[i]-New_curved_unit_cell_information[row])
                length.append(temp)

        return min(length)

########################################################################

# write new POSCAR file
class WritePOSCAR:
    def __init__(self,MakeBubble,output,ExpandZ=None,NoOverlimit=None):
        self.output=output
        
        self.Read_previous_data(MakeBubble)
        
        if NoOverlimit is None:
            self.NoOverlimit=True
        else:
            self.NoOverlimit=NoOverlimit
            
        
        if len(self.Distored_unit_cell_info)==0:
            print("They are no unit cell information to save output file")
            sys.exit(2)
        
        if ExpandZ is None or ExpandZ is True:
            self.Expand_Zaxis()
        
        self.Write_new_POSCAR()
    
    ####################################################################
    def Read_previous_data(self,ClassMakeBubble):
        self.Distored_unit_cell_info    = ClassMakeBubble.Distored_unit_cell_info
        self.unit_cell_para             = ClassMakeBubble.unit_cell_para
        self.atom_data                  = ClassMakeBubble.atom_data
        self.atom_data_num              = ClassMakeBubble.atom_data_num
        self.limit_range_atom           = ClassMakeBubble.limit_range_atom
        self.unit_cell_info             = ClassMakeBubble.unit_cell_info
    
    
    def Expand_Zaxis(self):
        limit_z=0.8
        if(max(self.Distored_unit_cell_info[:,2])>limit_z):
            print('\nThe maxium value is over',limit_z, 'So, we have to change the Z value')
            print('Before unit_cell_para :\n',self.unit_cell_para)            
            print('Before max :',max(self.Distored_unit_cell_info[:,2]))
            print('Before min :',min(self.Distored_unit_cell_info[:,2]))
            
            
            diff_cell_z=(max(self.Distored_unit_cell_info[:,2])-1)+(max(self.Distored_unit_cell_info[:,2]-0.5)*(1-limit_z)/(limit_z-0.5))
            plus_z=np.round_(diff_cell_z*np.array(self.unit_cell_para[2]),2)
            new_unit_cell_parameter_Z=self.unit_cell_para[2]+plus_z*2

            for i in range(len(self.Distored_unit_cell_info)):
                origin_z=(self.Distored_unit_cell_info[i][2]*np.array(self.unit_cell_para[2]))
                new_z=origin_z+plus_z
                self.Distored_unit_cell_info[i][2]=LA.norm(new_z)/LA.norm(new_unit_cell_parameter_Z)

            self.unit_cell_para[2]=new_unit_cell_parameter_Z.tolist()

            print('\n!After unit_cell_para :\n',self.unit_cell_para)            
            print('!After max :',max(self.Distored_unit_cell_info[:,2]))
            print('!After min :',min(self.Distored_unit_cell_info[:,2]))

        else:
            print('Maxium or Minium Z of cell information : ')
            print('max :',max(self.Distored_unit_cell_info[:,2]))
            print('min :',min(self.Distored_unit_cell_info[:,2]))
    

    def Write_new_POSCAR(self):
        print('\n!Output File name is : ',self.output)
        with open(self.output, 'w') as f:
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
            
            for i in range(len(self.Distored_unit_cell_info)):
                if i in self.limit_range_atom:
                    line="{:>10.9f}    {:>10.9f}    {:>10.9f}  0  0  0\n".format(
                        self.Inside_limitRange(self.Distored_unit_cell_info[i][0]),
                        self.Inside_limitRange(self.Distored_unit_cell_info[i][1]),
                        self.Inside_limitRange(self.Distored_unit_cell_info[i][2]))
                else:
                    line="{:>10.9f}    {:>10.9f}    {:>10.9f}\n".format(
                        self.Inside_limitRange(self.Distored_unit_cell_info[i][0]),
                        self.Inside_limitRange(self.Distored_unit_cell_info[i][1]),
                        self.Inside_limitRange(self.Distored_unit_cell_info[i][2]))
                f.write(line)
        f.close()        
        
        
    def Write_NoShift_POSCAR(self,outputfile=None):
        if outputfile is None:
            outputfile = self.ouput
        
        #Get atom name in No shift atom list
        list_atom=[]
        for i in self.limit_range_atom:
            list_atom.append(self.unit_cell_info[i][0])
        #Count the overlap atom
        new_atom_list=self.Get_overlap_element(list_atom)
        #sort the atom order along new_atom_list
        temp_limit_range=[]
        for key in new_atom_list.keys():
            for i in range(len(self.limit_range_atom)):
                temp=self.limit_range_atom[i]
                if self.unit_cell_info[temp][0] == key:
                    temp_limit_range.append(temp)
        if len(self.limit_range_atom) != len(temp_limit_range):
            print("The limit_range atom : ",len(self.limit_range_atom))
            print("The temp_limit_range : ",len(temp_limit_range))
            print("Plz, Check limit_range_atom code line!!!")
            sys.exit(2)
        
        
        print('\n!Output File(NoShift POSCAR) name is : ',outputfile)
        with open(outputfile, 'w') as f:
            f.write("POSCAR\n")
            f.write("1.00000\n")
            
            for i in range(len(self.unit_cell_para)):
                line="{:>10.7f}    {:>10.7f}    {:>10.7f}\n".format(
                    float(self.unit_cell_para[i][0]),
                    float(self.unit_cell_para[i][1]),
                    float(self.unit_cell_para[i][2]))
                f.write(line)
            
            for key in new_atom_list.keys():
                f.write("{} ".format(key))
            f.write("\n")

            for num in new_atom_list.values():
                f.write("{} ".format(num))
            f.write("\n")
            
            f.write("Direct\n")
            
            for i in temp_limit_range:
                line="{:>10.9f}    {:>10.9f}    {:>10.9f}  0  0  0\n".format(
                    self.Inside_limitRange(self.Distored_unit_cell_info[i][0]),
                    self.Inside_limitRange(self.Distored_unit_cell_info[i][1]),
                    self.Inside_limitRange(self.Distored_unit_cell_info[i][2]))
                f.write(line)
        f.close()        
    
    def Write_OnlyShift_POSCAR(self,outputfile=None):
        if outputfile is None:
            outputfile = self.ouput
        
        #Get atom name in Only shift atom list
        list_atom=[]
        for i in range(len(self.Distored_unit_cell_info)):
            if i not in self.limit_range_atom:
                list_atom.append(self.unit_cell_info[i][0])
        #Count the overlap atom
        new_atom_list=self.Get_overlap_element(list_atom)
        #sort the atom order along new_atom_list
        temp_limit_range=[]
        for key in new_atom_list.keys():
            for i in range(len(self.Distored_unit_cell_info)):
                if i not in self.limit_range_atom:
                    temp=i
                    if self.unit_cell_info[temp][0] == key:
                        temp_limit_range.append(temp)
        if (len(self.limit_range_atom)+len(temp_limit_range)) != len(self.Distored_unit_cell_info):
            print("The limit_range atom + The temp_limit_range : ",len(self.limit_range_atom)," + ",len(temp_limit_range))
            print("The Distored_unit_cell_info : ",len(self.Distored_unit_cell_info))
            print("Plz, Check limit_range_atom code line!!!")
            sys.exit(2)
  
        print('\n!Output File(NoShift POSCAR) name is : ',outputfile)
        with open(outputfile, 'w') as f:
            f.write("POSCAR\n")
            f.write("1.00000\n")
            
            for i in range(len(self.unit_cell_para)):
                line="{:>10.7f}    {:>10.7f}    {:>10.7f}\n".format(
                    float(self.unit_cell_para[i][0]),
                    float(self.unit_cell_para[i][1]),
                    float(self.unit_cell_para[i][2]))
                f.write(line)
            
            for key in new_atom_list.keys():
                f.write("{} ".format(key))
            f.write("\n")

            for num in new_atom_list.values():
                f.write("{} ".format(num))
            f.write("\n")
            
            f.write("Direct\n")
            
            for i in temp_limit_range:
                line="{:>10.9f}    {:>10.9f}    {:>10.9f}\n".format(
                    self.Inside_limitRange(self.Distored_unit_cell_info[i][0]),
                    self.Inside_limitRange(self.Distored_unit_cell_info[i][1]),
                    self.Inside_limitRange(self.Distored_unit_cell_info[i][2]))
                f.write(line)
        f.close()        
    
    def Get_overlap_element(self,list_atom):
        new_list={}
        for i in list_atom:
            try: new_list[i] +=1
            except: new_list[i] =1
        return new_list

    def Inside_limitRange(self,fraction):
        if self.NoOverlimit is True:
            if float(fraction) >= 1:
                return float(fraction) - 1
            elif float(fraction) < 0:
                return float(fraction) + 1
            else:
                return float(fraction)
        

        
    ####################################################################

########################################################################
# # test module code
inputfile='POSCAR_UnitCell_h_BN_pbesol_30'

hightest_Z=3.5    #angstrom
FWHM=10         #angstrom, r_cut, upper_flat
value_slop=0.5

output='POSCAR_pillarshape_Z_'+NPB.str_point2bar(hightest_Z)+'_FWHM_'+str(FWHM)+'_slop_'+NPB.str_point2bar(value_slop)
outputnoShift='POSCAR_noshift_pillarshape_Z_'+NPB.str_point2bar(hightest_Z)+'_FWHM_'+str(FWHM)
outputOnlyShift='POSCAR_Onlyshift_pillarshape_Z_'+NPB.str_point2bar(hightest_Z)+'_FWHM_'+str(FWHM)

#read the input data
A=AtomData(inputfile)

#modify the bubble shape
B=MakeBubble(A,hightest_Z,FWHM)
# B.MakeBubbleStructure()
# B.MakeConeStructure()
B.MakePillarStructure(slop=value_slop)
B.SearchNoShiftArea()

#write the bubble shape
C=WritePOSCAR(B,output)
C.Write_NoShift_POSCAR(outputnoShift)
C.Write_OnlyShift_POSCAR(outputOnlyShift)
########################################################################