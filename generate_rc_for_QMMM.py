#!/bin/python
import argparse
import os.path
import math
def get_pos(fxyz,ser,pos):
    f=open(fxyz,'r')
    coor=0
    out=0
    if (pos=="x"):
        coor=1
        
    if (pos=="y"):
        coor=2
        
    if (pos=="z"):
        coor=3
    
    if (pos=="name"):
        coor=0
    	       
    atm=0
    count=0
    for line in f:
        count += 1
        if (count==1):
            natom=line.split()[0]
# second line in xyz is a title    
        if (count>2):
            atm += 1
        
        if (atm==ser):
            out=line.split()[coor]
    
    f.close()
    if (ser==0) or (pos=="natom"):
        return natom
    else:
        return out 
            
            

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_name", dest='input', default='init.xyz', help='input xyz file')
parser.add_argument("-o", "--output_prefix", dest='outpre', default='moved', help='output prefix name')
parser.add_argument("-w", "--window_size", dest='ws', default=0.1, help='window size for moving atoms default 0.1 ang')
parser.add_argument("-di", "--initial_dist", dest='di', help='initial distance')
parser.add_argument("-df", "--final_dist", dest='df', help='final distance')
parser.add_argument("-mov", "--move_atoms", dest='mov', help='serial number of atoms to move together. hint. enter the interacting atom at first. the atoms should be separated by space in quotes')
parser.add_argument("-A", "--A_atom", dest='Aatom', help="reference A atom, serial number")
parser.add_argument("-B", "--B_atom", dest='Batom', help="reference B atom serial number. hint. B atom is close to the atom you are moving")

args = parser.parse_args()

ws=float(args.ws)
mov=args.mov.split()

inp=args.input
out=args.outpre
A=int(args.Aatom)
B=int(args.Batom)
di=float(args.di)
df=float(args.df)

if os.path.isfile(inp):
    print("file {} exists. good!".format(inp))
    print("input parameters")
    print("A atom: {}".format(str(A)))
    print("B atom: {}".format(str(B)))
    print("d0    : {}".format(str(di)))
    print("df    : {}".format(str(df)))
    print("move  : {}".format(mov))
else:
    print("file {} does not exist. bad!".format(inp))
    exit()

ax=float(get_pos(inp,A,'x'))
ay=float(get_pos(inp,A,'y'))
az=float(get_pos(inp,A,'z'))
aname=get_pos(inp,A,'name')


bx=float(get_pos(inp,B,'x'))
by=float(get_pos(inp,B,'y'))
bz=float(get_pos(inp,B,'z'))
aname=get_pos(inp,B,'name')

natom=int(get_pos(inp,0,'natom'))

orig={}
counti=0
for i in mov:
    i=int(i)
    orig[i]={}
    orig[i]['x']=float(get_pos(inp,i,'x'))
    orig[i]['y']=float(get_pos(inp,i,'y'))
    orig[i]['z']=float(get_pos(inp,i,'z'))
    orig[i]['name']=get_pos(inp,i,'name')
    if (counti==0):
        orig0=i
    
    counti += 1
    
    
abx=bx-ax
aby=by-ay
abz=bz-az

sqrt=math.sqrt
abdist=sqrt((abx**2)+(aby**2)+(abz**2))

vecx=abx/abdist
vecy=aby/abdist
vecz=abz/abdist

i=di
counti=0
while (i<=df):
    i=di+(counti*ws)
    fout=open(out+"-d"+str(counti)+".xyz", 'w')
    fout.write('{}\n'.format(str(natom)))
    fout.write('\t X-O distance {}\n'.format(str(i)))
    
    new={}
    countj=0
    for j in mov:
        j=int(j)
        new[j]={}
        if (countj==0):
            new[j]['x']=bx+(vecx*i)
            new[j]['y']=by+(vecy*i)
            new[j]['z']=bz+(vecz*i)
            atm0=j
        else:
            new[j]['x']=new[atm0]['x']-orig[orig0]['x']+orig[j]['x']
            new[j]['y']=new[atm0]['y']-orig[orig0]['y']+orig[j]['y']
            new[j]['z']=new[atm0]['z']-orig[orig0]['z']+orig[j]['z']
        
        new[j]['name']=get_pos(str(inp),int(j),'name')
        print("{}\t{}\t{}\t{}\t{}\t".format(str(j),new[j]['name'],str(new[j]['x']),str(new[j]['y']),str(new[j]['z'])))
        countj += 1
        #print("{}\t{}\t{}\t{}\t{}\t".format('67',get_pos('init.xyz',67,"name"),'x','y','z'))
        
    k=1
    while (k<=natom):
        if k in new:
            print('found atom {} to move\n'.format(str(k)))	
            fout.write('{}\t{}\t{}\t{}\n'.format(new[k]['name'],new[k]['x'],new[k]['y'],new[k]['z']))
        else:
            name=get_pos(inp,k,'name')
            x=get_pos(inp,k,'x')
            y=get_pos(inp,k,'y')
            z=get_pos(inp,k,'z')
            fout.write('{}\t{}\t{}\t{}\n'.format(name,x,y,z))
            
        k += 1
                
    counti += 1
    fout.write("\n")
    fout.close()
    
    
