import argparse
import os
import numpy as np
import pickle 
import units

parameters = units.get_parameters()

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-phi',type=float,default=1.0,help='packing fraction')
parser.add_argument('-v',type=float,default=0.1,help='propulsion speed')
parser.add_argument('-f',type=float,default=300.0,help='constant force outside channel')
parser.add_argument('-al',type=float,default=1.0,help='slowdown factor')
parser.add_argument('-cycle',type=float,default=0.1,help='cycle time')
parser.add_argument('-run',type=int,default=2,help='total run time')
parser.add_argument('-dt',type=float,default=2e-6,help='time step')
parser.add_argument('-outdir',type=str,default='./out',help='output directory')
args = parser.parse_args()

outdir = args.outdir 
phi = args.phi
v = args.v
force = args.f
dt = args.dt
sig = 1.0
rcut = 2**(1./6)*sig
eps = 100.0
frcut = 1.5*sig
al = args.al
D0 = parameters['Dt']
taur = 1./parameters['Dr']
rcon = parameters['rcon']
width = parameters['width']
run = args.run
cycle = args.cycle

filename = 'v{}phi{}f{}al{}.par'.format(args.v,phi,force,al)
f = open(filename,'w')
f.write("outdir\t{}\n".format(outdir))
f.write("phi\t{}\n".format(phi))
f.write("sig\t{}\n".format(sig))
f.write("rcut\t{}\n".format(rcut))
f.write("eps\t{}\n".format(eps))
f.write("frcut\t{}\n".format(frcut))
f.write("al\t{}\n".format(al))
f.write("v0\t{}\n".format(v))
f.write("D0\t{}\n".format(D0))
f.write("taur\t{}\n".format(taur))
f.write("rcon\t{}\n".format(rcon))
f.write("width\t{}\n".format(width))
f.write("f\t{}\n".format(force))
f.write("run\t{}\n".format(run))
f.write("cycle\t{}\n".format(cycle))
f.write("dt\t{}\n".format(dt))
f.write("t0\t{}\n".format(parameters['t0']))
f.write("l0\t{}\n".format(parameters['l0']))
f.close()

os.system('./sim '+filename)

r_all = []
phi_all = []
for i in range(int(run/cycle)):
    x,y,Phi = np.loadtxt('{}/v{}phi{}f{}al{}_{}.dat'.format(outdir,v,phi,force,al,i),unpack=True)
    r_all.append(np.column_stack((x,y)))
    phi_all.append(Phi)
r_all = np.array(r_all)
phi_all = np.array(phi_all)

traj = '{}/v{}phi{}f{}al{}.pickle'.format(outdir,v,phi,force,al)
data = {}
data['r'] = r_all
data['phi'] = phi_all
with open(traj,'wb') as handle:
    pickle.dump(data,handle)

os.system('rm {}/v{}phi{}f{}al{}_*.dat'.format(outdir,v,phi,force,al))
