#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 18:51:37 2017

@author: Romain
"""

import os
from shutil import copyfile
import unicodedata
from database.create_atom import atom


core = int(input("number of core : "))
charge = int(input("total charge : "))
label = []
for i in range(core):
    val = input("core label "+ str(i+1) +" : ")
    label.append(val)
netat = int(input("number of states : "))
algo = input("use FullCI (t/f) : ")
if(algo!='t') and (algo!='f'):
    raise Exception("invalid algorithm")
methodf = input("use f-orbitals (y/n) : ")
if(methodf!='y') and (methodf!='n'):
    raise Exception("invalid method")
elif(methodf=='y'):
    orbf=True
elif(methodf=='n'):
    orbf=False

homo=False
if(core==2)and(label[0]==label[1]):
    homo=True
    label=[label[0]]
print(label)
output=""
atoms=""
name=""
basis=""
basisf=""
psdinp=""
psdinpf=""
rcut=""
calfa=""
active=1
if(algo=='t'):active=3
electron=0
n=0
for val in label:
    n+=1
    output += val
    atoms += "'"+val+"',"
    name += unicodedata.normalize("NFKD", val.casefold())
    elem=atom(val)
    electron+=elem.ele_val
    number = elem.number
    (coef,rcutval,pola,pseudo)=elem.get_basis('cipsi')
    s_shell=len(coef[0])
    p_shell=len(coef[1])
    d_shell=len(coef[2])
    f_shell=len(coef[3])
    if(pseudo[0]=="'pssl'"): psdinp = pseudo[1]
    basis += "&CENT\n"
    basis += " ATOM='"+str(val)+"',IATN="+str(number)+",ZNUC="+str(elem.ele_val)+",\n"
    basis += " APSD="+pseudo[0]+",\n"

    if(f_shell>0):
        basis += " NTYPG="+str(s_shell)+"*1,"+str(p_shell)+"*2,"+str(d_shell)+"*3,"+str(f_shell)+"*4,\n"
        basis += " NCONG="+str(s_shell)+"*1,"+str(p_shell)+"*1,"+str(d_shell)+"*1,"+str(f_shell)+"*1,\n"
        basis += " EXG="
        for i in range(4):
            for expo in coef[i]:
                basis += str(expo)+","
            basis +="\n"
        basis += " CSG="+str(s_shell)+"*1,\n"
        basis += " CPG="+str(p_shell)+"*1,\n"
        basis += " CDG="+str(d_shell)+"*1,\n"
        basis += " CFG="+str(f_shell)+"*1,\n"
    else:
        basis += " NTYPG="+str(s_shell)+"*1,"+str(p_shell)+"*2,"+str(d_shell)+"*3,\n"
        basis += " NCONG="+str(s_shell)+"*1,"+str(p_shell)+"*1,"+str(d_shell)+"*1,\n"
        basis += " EXG="
        for i in range(3):
            for expo in coef[i]:
                basis += str(expo)+","
            basis +="\n"
        basis += " CSG="+str(s_shell)+"*1,\n"
        basis += " CPG="+str(p_shell)+"*1,\n"
        basis += " CDG="+str(d_shell)+"*1,\n"
    basis += "&END\n\n"

    calfa += str(pola)+","
    rcut += " RCUT1("+str(n)+",0)="+str(rcutval[0])
    rcut += ",RCUT1("+str(n)+",1)="+str(rcutval[1])
    rcut += ",RCUT1("+str(n)+",2)="+str(rcutval[2])
    rcut += ",RCUT1("+str(n)+",3)="+str(rcutval[3])+",\n"

    (coef,temp,temp1,pseudo)=elem.get_basis('cipsi-nof')
    s_shell=len(coef[0])
    p_shell=len(coef[1])
    d_shell=len(coef[2])
    if(pseudo[0]=="'pssl'"): psdinpf = pseudo[1]
    basisf += "&CENT\n"
    basisf += " ATOM='"+str(val)+"',IATN="+str(number)+",ZNUC="+str(elem.ele_val)+",\n"
    basisf += " APSD="+pseudo[0]+",\n"
    basisf += " NTYPG="+str(s_shell)+"*1,"+str(p_shell)+"*2,"+str(d_shell)+"*3,\n"
    basisf += " NCONG="+str(s_shell)+"*1,"+str(p_shell)+"*1,"+str(d_shell)+"*1,\n"
    basisf += " EXG="
    for i in range(3):
        for expo in coef[i]:
            basisf += str(expo)+","
        basisf +="\n"
    basisf += " CSG="+str(s_shell)+"*1,\n"
    basisf += " CPG="+str(p_shell)+"*1,\n"
    basisf += " CDG="+str(d_shell)+"*1,\n"
    basisf += "&END\n\n"

if(homo):
    electron *= 2
    name += '2'
    output += '2'
    atoms += "'"+label[0]+"',"
    calfa += str(pola)+","
    rcut += " RCUT1(2,0)="+str(rcutval[0])
    rcut += ",RCUT1(2,1)="+str(rcutval[1])
    rcut += ",RCUT1(2,2)="+str(rcutval[2])
    rcut += ",RCUT1(2,3)="+str(rcutval[3])+","
electron-=charge

namesym = name+"sym"
namedip = name+"dip"
if(orbf):nameorb=name
else:nameorb=name+'-nof'
atoms=atoms[:-1]
if(charge==1):output+="p"

if(core == 3):
    group = 'CS'
    nsym = 2
    sym = "'a1','a2'"
    ndip = 3
    dipole = "'a1','a1' \n 'a1','a2' \n 'a2','a2'"
elif(core == 2):
    if(homo):
      group = 'DINFH'
      nsym = 6
      sym = "'spg','spu','pxg','pxu','smg','smu'"
      ndip = 7
      dipole = """'spg','spu' \n'spg','pxu'\n'pxg','pxu'\n'pxg','smu'\n'smg','smu'
'spu','pxg'\n'pxu','smg'"""
    else:
      group = 'CINFV'
      nsym = 3
      sym = "'sp','px','sm'"
      ndip = 5
      dipole = "'sp','sp' \n'sp','px' \n'px','px'\n'px','sm'\n'sm','sm'"
else:
    group = 'C1'
    nsym = 1
    sym = "'a'"
    ndip = 1
    dipole = "'a','a'"

if(electron==1):
    spin = "0,1,0,0,0,0"
elif(electron==2):
    spin = "1,0,1,0,0,0"
elif(electron==3):
    spin = "0,1,0,1,0,0"
elif(electron==4):
    spin = "1,0,1,0,1,0"

dir_path = os.path.dirname(os.path.realpath(__file__))
template_file = dir_path+"../template/auto.in"
template_name = dir_path+"../template/name.dat"
template_dip = dir_path+"../template/dip.dat"
template_sym = dir_path+"../template/sym.dat"
out_file = 'input/auto.in'
out_name = 'input/'+name+'.dat'
out_namef = 'input/'+name+'-nof.dat'
out_sym = 'input/'+namesym+'.dat'
out_dip = 'input/'+namedip+'.dat'

os.makedirs(os.path.dirname(out_file), exist_ok=True)
os.makedirs(os.path.dirname(out_name), exist_ok=True)

with open(template_file,"r") as file_in:
    data = file_in.read().format(
        output = output
        ,core = core
        ,electron = electron
        ,atoms = atoms
        ,name = nameorb
        ,namesym = namesym
        ,namedip = namedip
        ,group = group
        ,nsym = nsym
        ,sym = sym
        ,ndip = ndip
        ,dipole = dipole
        ,spin = spin
        ,netat = netat
        ,algo = algo
        ,active = active
        )
with open(out_file,"w") as file_out:
    file_out.write(data)

with open(template_name,"r") as file_in:
    data = file_in.read().format(
        name = name
        ,basis = basis
        ,psdinp = psdinp
        ,rcut = rcut
        ,calfa = calfa
        )
with open(out_name,"w") as file_out:
    file_out.write(data)

with open(template_name,"r") as file_in:
    data = file_in.read().format(
        name = name
        ,basis = basisf
        ,psdinp = psdinpf
        ,rcut = rcut
        ,calfa = calfa
        )
with open(out_namef,"w") as file_out:
    file_out.write(data)

copyfile(template_sym,out_sym)
copyfile(template_dip,out_dip)
