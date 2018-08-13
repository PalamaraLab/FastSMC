import numpy as np
import argparse, sys, os
from smcpp import _smcpp, util, estimation_tools
from smcpp.model import *

parser = argparse.ArgumentParser(description='getCSFS')
parser.add_argument("-D", "--demographicFile", help="File with demographic model", action="store")
parser.add_argument("-n", "--samples", help="Number of distinguished+undistinguished lineages", action="store", type=int)
parser.add_argument("-o", "--out", help="Output file", action="store", default="")
parser.add_argument("-d", "--discretizationFile", help="File with time intervals where to compute the CSFS", action="store")
args=parser.parse_args()

# parse time and size
demo=np.loadtxt(args.demographicFile)
arrayTime=demo[:,0]
arraySize=demo[:,1]
arrayDisc=np.loadtxt(args.discretizationFile)

# add dummy last time to get np.diff
arrayTimeAppend = np.append(arrayTime, arrayTime[-1]+100)

# set n and N0
n = args.samples # number of total haploids, distinguished+undistinguished 
mu = 1.65E-8
N0 = arraySize[0]
# Population scaled mutation rate
theta = mu * 2. * N0
om = OldStyleModel(arraySize / (2. * N0), arraySize / (2. * N0), np.diff(arrayTimeAppend / (2. * N0)), N0)
arrayDiscOriginal = arrayDisc
arrayDisc = arrayDisc / (2. * N0)

f = open(args.out + ".csfs",'w')
arrayDisc = np.append(arrayDisc, np.inf)
arrayDiscOriginal = np.append(arrayDiscOriginal, np.inf)
for i in range(len(arrayDisc)-1):
        t0 = arrayDiscOriginal[i]
        t1 = arrayDiscOriginal[i+1]
        res = _smcpp.raw_sfs(om, n-2, t0 / (2. * N0), t1 / (2. * N0)) * theta
        res[0,0] = 1 - np.sum(res)
        if args.out != "":
                 f.write("Time:\t" + ' '.join(map(str, arrayTime)) + "\n")
                 f.write("Size:\t" + ' '.join(map(str, arraySize)) + "\n")
                 f.write("Mu:\t" + str(mu) + "\n")
                 f.write("Samples:\t" + str(args.samples) + "\n")
                 t0s = str(t0)
                 if (t1 != np.inf):
                     t1s = str(t1)
                 else:
                     t1s = "Infinity"
                 f.write("Interval:\t" + t0s + "\t" + t1s + "\n")
                 f.write('\n'.join(' '.join(str(cell) for cell in row) for row in res) + "\n")
f.close()
