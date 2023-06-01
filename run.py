import argparse
import os
import numpy as np
import pickle 
import units

parameters = units.get_parameters()

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-phi',type=float,default=1.2,help='packing fraction')
parser.add_argument('-v',type=float,default=0.1,help='propulsion speed')
parser.add_argument('-f',type=float,default=500.0,help='constant force outside channel')
parser.add_argument('-al',type=float,default=0.0,help='slowdown factor (translational diffusion)')
parser.add_argument('-alr',type=float,default=0.0,help='slowdown factor (rotational diffusion)')
parser.add_argument('-kap',type=float,default=1.0,help='alignment strength')
parser.add_argument('-eps',type=float,default=100.0,help='energy scale of WCA')
parser.add_argument('-weps',type=float,default=0.0,help='energy scale of WHDF')
parser.add_argument('-cycle',type=float,default=0.1,help='cycle time')
parser.add_argument('-run',type=int,default=50,help='total run time')
parser.add_argument('-dt',type=float,default=2e-5,help='time step')
parser.add_argument('-outdir',type=str,default='./out',help='output directory')
parser.add_argument('-potential',type=str,default='wca',help='interaction potential')
args = parser.parse_args()

outdir = args.outdir + "/" + args.potential 
potential = args.potential 
phi = args.phi
v = args.v
force = args.f
dt = args.dt
rcut = 1.0
whdf_rcut = 1.05*rcut
sig = rcut/2**(1./6)
eps = args.eps
whdf_eps = args.weps
frcut = 1.5*rcut
al = args.al
al_r = args.alr
kap = args.kap
D0 = parameters['Dt']
taur = 1./parameters['Dr']
rcon = parameters['rcon']
width = parameters['width']
run = args.run
cycle = args.cycle

filename = 'v{}phi{}f{}al{}alr{}kap{}'.format(args.v,phi,force,al,al_r,kap)
if potential=="att":
    filename = 'v{}phi{}f{}al{}alr{}weps{}rc{}kap{}'.format(args.v,phi,force,al,al_r,whdf_eps,whdf_rcut,kap)
f = open(filename+'.par','w')
f.write("outdir\t{}\n".format(outdir))
f.write("potential\t{}\n".format(potential))
f.write("phi\t{}\n".format(phi))
f.write("sig\t{}\n".format(sig))
f.write("eps\t{}\n".format(eps))
f.write("rcut\t{}\n".format(rcut))
f.write("whdf_rcut\t{}\n".format(whdf_rcut))
f.write("whdf_eps\t{}\n".format(whdf_eps))
f.write("frcut\t{}\n".format(frcut))
f.write("al\t{}\n".format(al))
f.write("al_r\t{}\n".format(al_r))
f.write("kap\t{}\n".format(kap))
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

os.system('./sim '+filename+'.par')

r_all = []
phi_all = []
for i in range(int(run/cycle)):
    x,y,Phi = np.loadtxt('{}/{}_{}.dat'.format(outdir,filename,i),unpack=True)
    r_all.append(np.column_stack((x,y)))
    phi_all.append(Phi)
r_all = np.array(r_all)
phi_all = np.array(phi_all)

traj = '{}/{}.pickle'.format(outdir,filename)
data = {}
data['r'] = r_all
data['phi'] = phi_all
with open(traj,'wb') as handle:
    pickle.dump(data,handle)

os.system('rm {}/{}_*.dat'.format(outdir,filename))
