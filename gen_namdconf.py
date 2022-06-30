#!/bin/python
import argparse
from random import randint

def get_pbc_info( pdbfile,pbcinfo ):
  pdb=open(pdbfile,'r')
  x=[]
  y=[]
  z=[]
  cen=[]
  for line in pdb:
    if (line[0:4]=="ATOM"):
      xyz=line[28:56].split()
      x.append(float(xyz[0]))
      y.append(float(xyz[1]))
      z.append(float(xyz[2]))

#print x
  n=len(x)
  o="{:.3f}".format(0.0)
  avgx=sum(x)/n
  cbx=max(x)-min(x)
  x1="{:.3f}".format(cbx)
  c1="{:.3f}".format(avgx)

  avgy=sum(y)/n
  cby=max(y)-min(y)
  y1="{:.3f}".format(cby)
  c2="{:.3f}".format(avgy)

  avgz=sum(z)/n
  cbz=max(z)-min(z)
  z1="{:.3f}".format(cbz)
  c3="{:.3f}".format(avgz)

  if (pbcinfo=="x"):
    return "{:<7} {:<7} {:<7}".format(x1,o,o)
  if (pbcinfo=="y"):
    return "{:<7} {:<7} {:<7}".format(o,y1,o)
  if (pbcinfo=="z"):
    return "{:<7} {:<7} {:<7}".format(o,o,z1)
  if (pbcinfo=="center"):
    return "{:<7} {:<7} {:<7}".format(c1,c2,c3)



parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_name", dest='input', help='input basename for psf and pdb')
parser.add_argument("-o", "--output_conf", dest='conf', help='output conf file name')
parser.add_argument("-en", "--ensemble", dest='ensemble', default='npt', help='ensemble to use: npt / nvt')
parser.add_argument("-T", "--temperature", dest='temperature', default='300', help='simulation temperature in K')
parser.add_argument("-re", "--restart", dest='restart', default='None', help='restart simulation? default: None')
parser.add_argument("-p", "--prmfile", action='append', dest='prmfile', help='parameter file (mandatory)')
parser.add_argument("-dt", "--dt", dest='dt', default='2', help='time step in fs')
parser.add_argument("-pbc", "--pbc", dest='pbc', default='yes', help= 'use pbc? (yes/no)')
parser.add_argument("-memb" "--memb", dest='memb', default='no', help='membrane simulation? (yes/no)')
parser.add_argument("-em", "--min", dest='min', default='1000', help='minimization steps')
parser.add_argument("-md", "--md", dest='md', default='5', help='simulation time in ns')
parser.add_argument("-eq", "--equilibration", dest='equil', default='no', help='equilibration run? (yes/no)')
parser.add_argument("-wf", "--writefreq", dest='writefreq', default=2500, type=int, help='write frequency for trajectory files')
parser.add_argument("-pme", "--pme", dest='pme', default='yes', help='use pme? (yes/no)')
parser.add_argument("-on", "--output_name", dest='outputname', default='md1', help='name for the output trajectory files')
parser.add_argument("-eb", "--ebonds", dest='ebonds', default='None', help='extrabonds file. Default: none')
parser.add_argument("-rvel", "--random_vel", dest='randvel', default='None', help='use random velocities? (yes/no)')

args = parser.parse_args()
fn = open(args.conf, 'w')
fn.write('{:<30} {}\n'.format('structure', args.input+'.psf'))
fn.write('{:<30} {}\n'.format('coordinates', args.input+'.pdb'))

if (args.restart == 'None'):
  fn.write('{:<30} {}\n'.format('temperature', args.temperature))
else:
  fn.write('{:<30} {}\n'.format('bincoordinates', args.restart+'.coor'))
  fn.write('{:<30} {}\n'.format('binvelocities', args.restart+'.vel'))
  fn.write('{:<30} {}\n'.format('extendedsystem', args.restart+'.xsc'))


fn.write('{:<30} {}\n'.format('firsttimestep', '0'))
fn.write('{:<30} {}\n'.format('paraTypeCharmm', 'on'))

for i in args.prmfile:
  fn.write('{:<30} {}\n'.format('parameters', i))

fn.write('\n')
fn.write('{:<30} {}\n'.format('exclude', 'scaled1-4'))
fn.write('{:<30} {}\n'.format('1-4scaling','1.0'))

if (args.pme == "yes"):
  cutoff=10.0
else:
  cutoff=14.0

swdist=cutoff-1.0
pdist=cutoff+4.0
  
fn.write('{:<30} {}\n'.format('cutoff', str(cutoff)))
fn.write('{:<30} {}\n'.format('switching','on'))
fn.write('{:<30} {}\n'.format('vdwForceSwitching','yes'))
fn.write('{:<30} {}\n'.format('switchdist', str(swdist)))
fn.write('{:<30} {}\n'.format('pairlistdist', str(pdist)))
fn.write('\n')
fn.write('{:<30} {}\n'.format('timestep',args.dt))
fn.write('{:<30} {}\n'.format('rigidbonds','all'))
fn.write('{:<30} {}\n'.format('nonbondedFreq','1'))
fn.write('{:<30} {}\n'.format('fullElectFrequency','2'))
fn.write('{:<30} {}\n'.format('stepspercycle','20'))
fn.write('\n')
fn.write('{:<30} {}\n'.format('langevin','on'))
if (args.equil=="yes"):
  ldamp="5"
else:
  ldamp="1"
  

fn.write('{:<30} {}\n'.format('langevinDamping', ldamp))
fn.write('{:<30} {}\n'.format('langevinTemp',args.temperature))
fn.write('{:<30} {}\n'.format('langevinHydrogen','off'))

if (args.ensemble=="nvt"):
  lpist='off'
else:
  lpist='on'
fn.write('\n')  
fn.write('{:<30} {}\n'.format('langevinPiston',lpist))
fn.write('{:<30} {}\n'.format('langevinPistonTarget','1.01325'))
fn.write('{:<30} {}\n'.format('langevinPistonPeriod','100.0'))
fn.write('{:<30} {}\n'.format('langevinPistonDecay','50.0'))
fn.write('{:<30} {}\n'.format('langevinPistonTemp',args.temperature))
fn.write('\n')
fn.write('{:<30} {}\n'.format('PME',args.pme))
fn.write('{:<30} {}\n'.format('PMEGridSpacing','1.0'))

# setting up periodic boundary condition values
if (args.pbc=="yes") and (args.restart=="None"):
  pdb=args.input+'.pdb'
  x=get_pbc_info(pdb,"x")
  y=get_pbc_info(pdb,"y")
  z=get_pbc_info(pdb,"z")
  cen=get_pbc_info(pdb,"center")
  
  fn.write('{:<30} {}\n'.format('CellbasisVector1',x))
  fn.write('{:<30} {}\n'.format('CellbasisVector2',y))
  fn.write('{:<30} {}\n'.format('CellbasisVector3',z))
  fn.write('{:<30} {}\n'.format('CellOrigin',cen))

if (args.pbc=="yes") or (args.restart!="None"):  
  fn.write('{:<30} {}\n'.format('wrapWater','on'))
  fn.write('{:<30} {}\n'.format('wrapAll','on'))

fn.write('{:<30} {}\n'.format('useGroupPressure','yes'))
if (args.memb=="yes"):
  flexcell="yes"
else:
  flexcell="no"
  
fn.write('{:<30} {}\n'.format('useFlexibleCell',flexcell))
fn.write('{:<30} {}\n'.format('useConstantRatio','no'))
fn.write('\n')  
fn.write('{:<30} {}\n'.format('dcdfreq',str(args.writefreq)))
fn.write('{:<30} {}\n'.format('xstfreq',str(args.writefreq)))
fn.write('{:<30} {}\n'.format('outputEnergies','500'))
fn.write('{:<30} {}\n'.format('outputPressure','500'))
fn.write('{:<30} {}\n'.format('restartfreq',str(2*args.writefreq)))
fn.write('\n')
if (args.ebonds!="None"):
  fn.write('{:<30} {}\n'.format('extrabonds','on'))
  fn.write('{:<30} {}\n'.format('extrabondsfile',args.ebonds))
  fn.write('\n')
  
fn.write('{:<30} {} \n'.format('outputName',args.outputname))

if (args.restart=="None"):
  fn.write('{:<30} {}\n'.format('minimize',args.min))

if (args.randvel=="yes"):
  fn.write('{:<30} {}\n'.format('reinitvels', str(args.temperature)))

totmd=int(float(args.md)*(10**6)/float(args.dt))
fn.write('{:<30} {}\n'.format('run',str(totmd)))
  


fn.close()


