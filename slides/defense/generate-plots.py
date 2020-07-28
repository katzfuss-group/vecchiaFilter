from colour import Color
import matplotlib as mpl
import numpy as np
import sys, os
sys.path.append('/home/marcin/pyMRA/pyMRA/')
sys.path.append('/home/marcin/MRF/python')
import MRATools as mt
import pdb
from MRA.MRD import MRD


Nx = 15
Ny = 15
locs = mt.genLocations2d(Nx=Nx, Ny=Ny)



colors = [(112.0/255, 115.0/255, 115.0/255), (1, 1, 1), (80.0/255, 0, 0)] # maroon to white
COLOR_MAP = mpl.colors.LinearSegmentedColormap.from_list("TAMU_map", colors, N=1000)



# create maroon color map for visualizing matrices
white = Color('white')
clist = list(white.range_to(Color('#800000'),256))
ccodes = [color.hex_l for color in clist]
cm = mpl.colors.LinearSegmentedColormap.from_list('maroon',ccodes)



cov = mt.ExpCovFun(locs, l=0.2)
mt.dispMat(cov, cmap=cm, colorbar=False, fName="images/trueCov.png")

M = 3; J = [4]; r = [3]
mrd = MRD(locs, M, J, r, cov)
B = mrd.getBasisFunctionsMatrix(timesKC=True, order="leaves")
color_end = max(np.max(B), np.abs(np.min(B)))

mt.dispMat(B*B.T, cmap=cm, colorbar=False, fName="images/BBT.png")

mt.dispMat(B, cmap=COLOR_MAP, vmax=color_end, vmin=-color_end, colorbar=False, fName="images/B.png")

mt.dispMat(B.T, cmap=COLOR_MAP, vmax=color_end, vmin=-color_end, colorbar=False, fName="images/BT.png")

mt.dispMat(mt.filterNNZ(B.T*B), colorbar=False, cmap=cm, fName="images/BTB.png")

BTBc = np.linalg.cholesky(B.T * B + np.eye(Nx*Ny))
mt.dispMat(mt.filterNNZ(BTBc), colorbar=False, cmap=cm, fName="images/BTBc.png")

BTBcInv = np.linalg.solve(BTBc, np.eye(Nx*Ny))
mt.dispMat(mt.filterNNZ(BTBcInv, tol=0.00), colorbar=False, cmap=cm, fName="images/BTBcInv.png")
