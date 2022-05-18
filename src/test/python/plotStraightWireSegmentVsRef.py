#!/usr/bin/env python

# Plot a comparison between demo outputs of a given ABSCAB implementation
# and the reference data provided in this repository.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator, MaxNLocator

# machine precision (ca. 2.22e-16 for 64-bit double)
eps = np.finfo(np.float64).eps

testKnotsRp = np.loadtxt("../resources/testKnotsRpStraightWireSegment.dat")
testKnotsZp = np.loadtxt("../resources/testKnotsZpStraightWireSegment.dat")

numR = len(testKnotsRp)
numZ = len(testKnotsZp)

idxRp = np.loadtxt("../resources/idxRpStraightWireSegment.dat", dtype=int)
idxZp = np.loadtxt("../resources/idxZpStraightWireSegment.dat", dtype=int)

numCases = len(idxRp)

# A_z
##ref1d = np.loadtxt("../resources/refDataStraightWireSegment.dat")[:,0]
##act1d = np.loadtxt("../../../data/straightWireSegmentAz.dat")

# B_phi
ref1d = np.loadtxt("../resources/refDataStraightWireSegment.dat")[:,1]
act1d = np.loadtxt("../../../data/straightWireSegmentBphi.dat")

ref = np.zeros([numZ, numR])
act = np.zeros([numZ, numR])

for i in range(numCases):
    ref[idxZp[i], idxRp[i]] = ref1d[i]
    act[idxZp[i], idxRp[i]] = act1d[i]

# compute rel. error between ref and act
good = -16

data = np.zeros([numZ, numR])
for i,rp in enumerate(testKnotsRp):
    for j,zp in enumerate(testKnotsZp):
        if rp == 0.0 and 0.0 <= zp and zp <= 1.0:
            # replace uninitialized points at location of wire segment with np.nan
            data[j,i] = np.nan
        elif abs(ref[j,i]) > 0:
            if act[j,i] != ref[j,i]:
                data[j,i] = np.log10(min(1, abs((act[j,i] - ref[j,i])/ref[j,i])))
            else:
                data[j,i] = good
        else:
            data[j,i] = 1.0 if abs(act[j,i]) > 0.0 else good
            
nClusters = 18
cmap = plt.get_cmap("viridis", nClusters)

fig=plt.figure(figsize=(4.5, 6.5))
ax = plt.gca()

im = plt.imshow(data, origin="lower", interpolation=None, cmap=cmap)

cbar = plt.colorbar(im, drawedges = True, fraction=0.0755, pad=-0.02,
                    anchor=(0,0.55), extend="min")
cbar.set_label(r"$\log_{10}(\mathrm{relative~deviation~from~reference})$")

plt.axis("off")

xTicksAndLabels = [
    [0.0,   r"$0$"],
    [1e-30, r"$10^{-30}$"],
    [1e-15, r"$10^{-15}$"],
    [1.0,   r"$1$"],
    [1e+15, r"$10^{15}$"],
    [1e+30, r"$10^{30}$"]
]

yTicksAndLabels = [
    [-1e+30,    r"$-10^{30}$"],
    [-1e+15,    r"$-10^{15}$"],
    [-1.0,      r"$-1$"],
    [-1e-15,    r"$-10^{-15}$"],
    [-1e-30,    r"$-10^{-30}$"], 
    [0.0,       r"$0$"],
    [1e-30,     r"$10^{-30}$"],
    [1e-15,     r"$10^{-15}$"],
    [1e-1,      r"$10^{-1}$"],
    [0.5,       r"$1/2$"],
    [1-(1e-1),  r"$1-10^{-1}$"],
    [1-(1e-15), r"$1-10^{-15}$"],
    [1-eps/2,   r"$1-\epsilon/2$"],
    [1.0,       r"$1$"],
    [1+eps,     r"$1+\epsilon$"],
    [1+(1e-15), r"$1+10^{-15}$"],
    [1+(1e-1),  r"$1+10^{-1}$"],
    [2.0,       r"$2$"],
    [1e+1,      r"$10^{1}$"],
    [1e+15,     r"$10^{15}$"],
    [1e+30,     r"$10^{30}$"]
]

# location and dimensions of imshow result
e = im.get_extent()
x0 = e[0]
y0 = e[2]
dx = e[1]-e[0]
dy = e[3]-e[2]

arrowHead = 3.0 # surplus length of arrow past extent of imshow result
headWidth = 2
headLength = 3
fc = "k" # black-filled arrow heads

# x axis arrow and label
ax.arrow(x0, y0, dx + arrowHead, 0.0, fc=fc,
         head_width=headWidth, head_length = headLength)
ax.text(dx+arrowHead+2, -5, r"$\rho'$")

# y axis arrow
ax.arrow(x0, y0, 0.0, dy + arrowHead, fc=fc,
         head_width=headWidth, head_length = headLength)
ax.text(-5, dy + arrowHead, r"$z'$")

majorTickLength = 3.5
minorTickLength = 2.0

# offset from end of tick to start of label text
labelOffsetX = 2.0
labelOffsetY = 2.0

# squeezed labels additional displacement
xAddUD = 3 # displacement in up/down direction
xAddLR = 5 # displacement in left/right direction
xAddTL = 3 # length of tick at end of displaced tick

yAddUD = 3.2 # displacement in up/down direction
yAddLR = 3 # displacement in left/right direction
yAddTL = 3 # length of tick at end of displaced tick
fac2 = 0.8 # shrinking factor for yAddLR at 2 off

# major ticks where wanted and minor ticks everywhere else
for i,k in enumerate(testKnotsRp):
    found = False
    for tl in xTicksAndLabels:
        if k == tl[0]:
            found = True

            if k == 0.0:
                # move "0" tick to bottom left

                # draw major x tick
                ax.plot([i, i, i-xAddLR, i-xAddLR],
                        [y0, -majorTickLength, -(majorTickLength+xAddUD), -(majorTickLength+xAddUD+xAddTL)],
                        "k-", lw=0.8)

                # print ticklabel
                # -1.5: need additional y space to account for dist from centerline to top
                ax.text(i-xAddLR, -majorTickLength - 1.5 - labelOffsetY - (xAddUD+xAddTL),
                        tl[1], ha='center', va='center_baseline')

            elif k == 1e-30:
                # move "1e-30" tick to bottom right
                
                # draw major x tick
                ax.plot([i, i, i+xAddLR, i+xAddLR],
                        [y0, -majorTickLength, -(majorTickLength+xAddUD), -(majorTickLength+xAddUD+xAddTL)],
                        "k-", lw=0.8)

                # print ticklabel
                # -1.5: need additional y space to account for dist from centerline to top
                ax.text(i+xAddLR, -majorTickLength - 1.5 - labelOffsetY - (xAddUD+xAddTL),
                        tl[1], ha='center', va='center_baseline')
            else:                
                # draw major x tick
                ax.plot([i, i], [y0, -majorTickLength], "k-", lw=0.8)

                # print ticklabel
                # -1.5: need additional y space to account for dist from centerline to top
                ax.text(i, -majorTickLength - 1.5 - labelOffsetY, tl[1], ha='center', va='center_baseline')

            # no need to search further...
            break

    if not found:
        # draw minor x tick
        ax.plot([i, i], [y0, -minorTickLength], "k-", lw=0.5)

for i,k in enumerate(testKnotsZp):
    found = False
    for tl in yTicksAndLabels:
        if k == tl[0]:
            found = True

            if k in [1-1e-15]:
                # move down-left by 2


                # draw major x tick
                ax.plot([x0, -fac2*majorTickLength, -(majorTickLength+yAddLR+(1-fac2)*yAddLR), -(majorTickLength+yAddLR+yAddLR)],
                        [i, i, i-2*yAddUD, i-2*yAddUD],
                        "k-", lw=0.8)

                # print ticklabel
                ax.text(-majorTickLength - labelOffsetX - (yAddLR+yAddLR),
                        i-2*yAddUD,
                        tl[1], ha='right', va='center_baseline')

            elif k in [-1e-30, 1e-1, 1-eps/2, 1+1e-1]:
                # move down-left by 1
                
                # draw major x tick
                ax.plot([x0, -majorTickLength, -(majorTickLength+yAddLR), -(majorTickLength+yAddLR+yAddLR)],
                        [i, i, i-yAddUD, i-yAddUD],
                        "k-", lw=0.8)

                # print ticklabel
                ax.text(-majorTickLength - labelOffsetX - (yAddLR+yAddLR),
                        i-yAddUD,
                        tl[1], ha='right', va='center_baseline')
                
            elif k in [0.0, 0.5, 1, 2]:
                # move left
                
                # draw major x tick
                ax.plot([x0, -(majorTickLength+yAddLR+yAddLR)],
                        [i, i],
                        "k-", lw=0.8)

                # print ticklabel
                ax.text(-majorTickLength - labelOffsetX - (yAddLR+yAddLR),
                        i,
                        tl[1], ha='right', va='center_baseline')
                
            elif k in [1e-30, 1-1e-1, 1+eps, 1e1]:
                # move up-left by 1

                # draw major x tick
                ax.plot([x0, -majorTickLength, -(majorTickLength+yAddLR), -(majorTickLength+yAddLR+yAddLR)],
                        [i, i, i+yAddUD, i+yAddUD],
                        "k-", lw=0.8)

                # print ticklabel
                ax.text(-majorTickLength - labelOffsetX - (yAddLR+yAddLR),
                        i+yAddUD,
                        tl[1], ha='right', va='center_baseline')

            elif k in [1+1e-15]:
                # move up-left by 2

                # draw major x tick
                ax.plot([x0, -fac2*majorTickLength, -(majorTickLength+yAddLR+(1-fac2)*yAddLR), -(majorTickLength+yAddLR+yAddLR)],
                        [i, i, i+2*yAddUD, i+2*yAddUD],
                        "k-", lw=0.8)

                # print ticklabel
                ax.text(-majorTickLength - labelOffsetX - (yAddLR+yAddLR),
                        i+2*yAddUD,
                        tl[1], ha='right', va='center_baseline')

            else:   
                # draw major y tick
                ax.plot([x0, -majorTickLength], [i, i], "k-", lw=0.8)

                # print ticklabel
                ax.text(-majorTickLength - labelOffsetX, i, tl[1], ha='right', va='center_baseline')

            # no need to search further...
            break
        
    if not found:
        # draw minor y tick
        ax.plot([x0, -minorTickLength], [i, i], "k-", lw=0.5)

plt.tight_layout()
plt.subplots_adjust(left=0.1,
                    right=0.8,
                    bottom=0.02,
                    top=0.98)

# EPS has no transparency
#plt.savefig(os.path.join(runFolder, "numCorrectDigits_newA.eps"), dpi=300)

##plt.savefig(os.path.join(runFolder, "numCorrectDigits_newA.pdf"), dpi=300)

plt.show()
