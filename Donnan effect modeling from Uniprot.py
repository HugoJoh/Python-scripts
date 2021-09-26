# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 17:23:04 2021

@author: HugoJoh
"""

#Open the data file
filename="P02768.fasta.txt"
a=open(filename)
rawdata=a.read()

#Data preparation
rawdata=rawdata.replace('\n','')
indx=rawdata.find("SV=")#position of 'SV='
rawdata=rawdata[indx:]
rawdata=rawdata[4:]

#Protein molecular weight calculation, use molecular formula, take in account that the entire aa is not present as there are bindings
Acount=rawdata.count('A')
Rcount=rawdata.count('R')
Ncount=rawdata.count('N')
Dcount=rawdata.count('D')
Ccount=rawdata.count('C')
Ecount=rawdata.count('E')
Qcount=rawdata.count('Q')
Gcount=rawdata.count('G')
Hcount=rawdata.count('H')
Icount=rawdata.count('I')
Lcount=rawdata.count('L')
Kcount=rawdata.count('K')
Mcount=rawdata.count('M')
Fcount=rawdata.count('F')
Pcount=rawdata.count('P')
Scount=rawdata.count('S')
Tcount=rawdata.count('T')
Wcount=rawdata.count('W')
Ycount=rawdata.count('Y')
Vcount=rawdata.count('V')
import numpy as np
list1 = [1,Rcount,Hcount,Kcount,Dcount,Ecount,Ycount,Ccount,1]
counts=np.array(list1)
#All aa without a H2O molecule ; add this H2O molecule at the end of the calculation
ProteinMM=18.01528+Acount*(89.1-18.01528)+Rcount*(174.2-18.01528)+Ncount*(132.1-18.01528)+Dcount*(133.1-18.01528)+Ccount*(121.2-18.01528)+Ecount*(147.1-18.01528)+Qcount*(146.2-18.01528)+Gcount*(75.1-18.01528)+Hcount*(155.2-18.01528)+Icount*(131.2-18.01528)+Lcount*(131.2-18.01528)+Kcount*(146.2-18.01528)+Mcount*(149.2-18.01528)+Fcount*(165.2-18.01528)+Pcount*(115.1-18.01528)+Scount*(105.1-18.01528)+Tcount*(119.1-18.01528)+Wcount*(204.2-18.01528)+Ycount*(181.2-18.01528)+Vcount*(117.1-18.01528)
print('Protein molecular weight:',round(ProteinMM,2),'g/mol')

#Create pHs
pHs=np.array(np.arange(3,11,0.01))

#Compute protein charges for 298K
pkaslist=[8.00,12.10,6.30,10.00,3.85,4.40,10.10,8.50,3.10]
pkas=np.array(pkaslist)
def chargefunc(x,y):
    if x > 4:
        if pkas[x] < pHs[y]:
            return -1 * counts[x]
        else:
            return 0
    elif x <= 4:
        if pkas[x] > pHs[y]:
            return 1 * counts[x]
        else:
            return 0
vchargefunc=np.vectorize(chargefunc)
list2 = []
Charges=np.array(list2)
for i in range(0,800):
    for j in range(0,9):
        Charges=np.append(Charges,vchargefunc(j,i))
rows, = Charges.shape
rows=int(rows/9)
Charges=Charges.reshape(rows,9)
Charges=np.sum(Charges,axis=1)

# Plot protein charges versus pHs
import matplotlib.pyplot as plt
Chargefit=np.polyfit(pHs,Charges,3)
a,b,c,d=Chargefit
import statistics as stat
list9=[]
Yhat=np.array(list9)
Yhat=a*pow(pHs,3)+b*pow(pHs,2)+c*pHs+d
R2=1-((sum(pow((Charges-Yhat),2)))/(sum(pow((Charges-stat.mean(Charges)),2))))
print('R2 of the "Protein charges = f(pH)" graph:',round(R2,2))
figure,axes= plt.subplots(4)
figure.set_figwidth(7)
figure.set_figheight(14)
figure.subplots_adjust(hspace=0.4)
axes[0].plot(pHs,Charges)
axes[0].plot(pHs,a*pow(pHs,3)+b*pow(pHs,2)+c*pHs+d)
axes[0].set_title('Protein charges = f(pH)')
axes[0].set(xlabel="pH",ylabel="Charges")

#Values to enter by the user
ProteinConc_during_DF_gbykg=70.0
print('[Protein] during DF:',round(ProteinConc_during_DF_gbykg,5),'g/L')
ProteinConc_during_DF_molbykg=ProteinConc_during_DF_gbykg/ProteinMM
print('[Protein] during DF:',round(ProteinConc_during_DF_molbykg,5),'mol/L')
DF_buffer_HistConc=10.0
print('DF buffer [Histidine]:',DF_buffer_HistConc,'mM')
DF_buffer_pH=5.80
print('DF buffer pH:',DF_buffer_pH)

#Donnan effect model per se
#Assumptions and notes:
# 1L solution = 1 kg solution for [Protein] in g/kg.
# Retentate [Cl-] is underestimated as Cl- ions from HCl are not taken in account if DF buffer pH differs from 6.
import math
ClConc_Start=2.4*1.02*0.0037/36.46 #As pH is lowered with a 3.7% HCl solution, this [Cl-] results from HCl solution volume added*HCl solution density*[HCl solution]/HCL MW
ClConc_Permeate=ClConc_Start
HisPlusConc_Start=(DF_buffer_HistConc/1000)/(1+(pow(10,(DF_buffer_pH-pkaslist[2]))))
HisConc_Start=DF_buffer_HistConc/1000-HisPlusConc_Start
HisConc_Retentate=HisConc_Start
HisConc_Permeate=HisConc_Retentate
HisPlusConc_Permeate=HisPlusConc_Start
list3=[28]
ProteinCharge=np.array(list3)
list4=[]
IRatio=np.array(list4)
list5=[]
DonnanRatio=np.array(list5)
list6=[]
ClConc_Retentate=np.array(list6)
list7=[]
HisPlusConc_Retentate=np.array(list7)
list8=[]
pH_Retentate=np.array(list8)
for i in range(0,1000):
    IRatio=np.append(IRatio,(ProteinCharge[i]*ProteinConc_during_DF_molbykg/HisPlusConc_Permeate))
    DonnanRatio=np.append(DonnanRatio,((-1)*IRatio[i]/2+math.sqrt(pow((IRatio[i]/2),2)+1)))
    ClConc_Retentate=np.append(ClConc_Retentate,(pow(DonnanRatio[i]*(pow(ClConc_Start,-1)),-1)))
    HisPlusConc_Retentate=np.append(HisPlusConc_Retentate,DonnanRatio[i]*(pow(HisPlusConc_Permeate,1)))
    pH_Retentate=np.append(pH_Retentate,pkaslist[2]+math.log10(HisConc_Retentate/HisPlusConc_Retentate[i]))
    ProteinCharge=np.append(ProteinCharge,a*pow(pH_Retentate[i],3)+b*pow(pH_Retentate[i],2)+c*pH_Retentate[i]+d)

#Outputting values
axes[1].plot(range(0,1000),HisPlusConc_Retentate)
axes[1].set_title('[His+] in retentate = f(iteration)')
axes[1].set(xlabel="iteration",ylabel="[His+] in mol/kg")
axes[2].plot(range(0,1000),ClConc_Retentate)
axes[2].set_title('[Cl-] in retentate = f(iteration)')
axes[2].set(xlabel="iteration",ylabel="[Cl-] in mol/kg")
axes[3].plot(range(0,1000),pH_Retentate)
axes[3].set_title('Retentate pH = f(iteration)')
axes[3].set(xlabel="iteration",ylabel="pH")
print('[His] in retentate:',round((HisConc_Retentate),5),'mol/kg')
print('[His] in permeate:',round((HisConc_Permeate),5),'mol/kg')
print('[His+] in retentate:',round((stat.mean(HisPlusConc_Retentate[997:1000])),5),'mol/kg')
print('[His+] in permeate:',round((HisPlusConc_Permeate),5),'mol/kg')
print('[Cl-] in permeate:',round((ClConc_Permeate),5),'mol/kg')
print('[Cl-] in retentate:',round((stat.mean(ClConc_Retentate[997:1000])),5),'mol/kg')
print('Retentate pH:',round((stat.mean(pH_Retentate[997:1000])),2))
print('[Histidine] in retentate:',round((((stat.mean(HisPlusConc_Retentate[997:1000]))+(HisConc_Retentate))*1000),5),'mM')