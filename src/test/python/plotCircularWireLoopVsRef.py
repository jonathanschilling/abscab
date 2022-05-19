#!/usr/bin/env python

# Plot a comparison between demo outputs of a given ABSCAB implementation
# and the reference data provided in this repository.

if __name__ == "__main__":
    import sys

    numArgs = len(sys.argv)
    if numArgs < 3:
        print("usage: " + sys.argv[0] + " <ref> <act> [savefig filename]")
        print(" where <ref> is the reference data (probably '../resources/CircularWireLoop_*_*_ref.dat')\n" +
              " and   <act> is the output to test (probably '../../../data/CircularWireLoop_*_*_*.dat')\n" +
              " Optionally, the generated plot is saved to the given filename (3rd parameter).");
        sys.exit(1)

    # filename of reference data
    refFilename = sys.argv[1]
    print("loading reference data from '%s'"%(refFilename,))
    
     # filename of actual data file to test
    actFilename = sys.argv[2]
    print("loading data to test from '%s'"%(actFilename,))
    
    # optional filename into which to save the plot
    savefigFilename = None
    if numArgs > 3:
        savefigFilename = sys.argv[3]
        print("will save plot to '%s'"%(savefigFilename,))

    # here we go...
    
    # lazy loading
    import numpy as np
    
    # machine precision (ca. 2.22e-16 for 64-bit double)
    eps = np.finfo(np.float64).eps
    
    testKnotsRp = np.loadtxt("../resources/testKnotsRpCircularWireLoop.dat")
    testKnotsZp = np.loadtxt("../resources/testKnotsZpCircularWireLoop.dat")
    
    numR = len(testKnotsRp)
    numZ = len(testKnotsZp)
    
    idxRp = np.loadtxt("../resources/idxRpCircularWireLoop.dat", dtype=int)
    idxZp = np.loadtxt("../resources/idxZpCircularWireLoop.dat", dtype=int)
    
    numCases = len(idxRp)

    # A_phi
    # ref1d = np.loadtxt("../resources/CircularWireLoop_A_phi_ref.dat")
    # act1d = np.loadtxt("../../../data/CircularWireLoop_A_phi_Java.dat")
    
    # B_rho
    # ref1d = np.loadtxt("../resources/CircularWireLoop_B_rho_ref.dat")
    # act1d = np.loadtxt("../../../data/CircularWireLoop_B_rho_Java.dat")
    
    # B_z
    # ref1d = np.loadtxt("../resources/CircularWireLoop_B_z_ref.dat")
    # act1d = np.loadtxt("../../../data/CircularWireLoop_B_z_Java.dat")
    
    ref1d = np.loadtxt(refFilename)
    act1d = np.loadtxt(actFilename)

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
            if rp == 1.0 and zp == 0.0: # circular wire loop
                # replace uninitialized points at location of wire segment with np.nan
                data[j,i] = np.nan
            elif abs(ref[j,i]) > 0:
                if act[j,i] != ref[j,i]:
                    data[j,i] = np.log10(min(1, abs((act[j,i] - ref[j,i])/ref[j,i])))
                else:
                    data[j,i] = good
            else:
                data[j,i] = 1.0 if abs(act[j,i]) > 0.0 else good
        
    # lazy loading
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import MultipleLocator, MaxNLocator
    
    nClusters = 18
    cmap = plt.get_cmap("viridis", nClusters)
    
    fig=plt.figure(figsize=(6.5, 4.5))
    ax = plt.gca()
    
    im = plt.imshow(data, origin="lower", interpolation=None, cmap=cmap)
    
    cbar = plt.colorbar(im, drawedges = True, fraction=0.025, pad=-0.02,
                        anchor=(0, 0.53), extend="min")
    cbar.set_label(r"$\log_{10}(\mathrm{relative~deviation~from~reference})$")
    
    plt.axis("off")
    
    xTicksAndLabels = [
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
    
    yTicksAndLabels = [
        [0.0,   r"$0$"],
        [1e-30, r"$10^{-30}$"],
        [1e-15, r"$10^{-15}$"],
        [1.0,   r"$1$"],
        [1e+15, r"$10^{15}$"],
        [1e+30, r"$10^{30}$"]
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
    ax.arrow(x0, y0, dx + arrowHead, 0.0, fc=fc, head_width=headWidth, head_length = headLength)
    ax.text(dx+arrowHead, -5, r"$\rho'$")
    
    # y axis arrow
    ax.arrow(x0, y0, 0.0, dy + arrowHead, fc=fc, head_width=headWidth, head_length = headLength)
    ax.text(-5, dy + arrowHead, r"$z'$")
    
    
    majorTickLength = 3.5
    minorTickLength = 2.0
    
    # offset from end of tick to start of label text
    labelOffsetX = 2.0
    labelOffsetY = 2.0
    
    # squeezed labels additional displacement
    xAddUD = 3 # displacement in up/down direction
    xAddLR = 3 # displacement in left/right direction
    xAddTL = 2 # length of tick at end of displaced tick
    fac2 = 0.8 # shrinking factor for yAddLR at 2 off
    
    yAddUD = 3 # displacement in up/down direction
    yAddLR = 2 # displacement in left/right direction
    yAddTL = 2 # length of tick at end of displaced tick
    
    # major ticks where wanted and minor ticks everywhere else
    for i,k in enumerate(testKnotsRp):
        found = False
        for tl in xTicksAndLabels:
            if k == tl[0]:
                found = True
    
                if k in [1-1e-15]:
                    # move down-left by 2
    
                    # draw major x tick
                    ax.plot([i, i, i-2*xAddLR, i-2*xAddLR],
                            [y0, -fac2*majorTickLength, -(majorTickLength+xAddUD+(1-fac2)*xAddTL), -(majorTickLength+xAddUD+xAddTL)],
                            "k-", lw=0.8)
    
                    # print ticklabel
                    ax.text(i-2*xAddLR,
                            -majorTickLength - labelOffsetY - (xAddUD+xAddTL) - 1.5,
                            tl[1], ha='center', va='center_baseline', rotation=-90)
    
                elif k in [1e-1, 1-eps/2, 1+1e-1]:
                    # move to bottom left
    
                    # draw major x tick
                    ax.plot([i, i, i-xAddLR, i-xAddLR],
                            [y0, -majorTickLength, -(majorTickLength+xAddUD), -(majorTickLength+xAddUD+xAddTL)],
                            "k-", lw=0.8)
    
                    # print ticklabel
                    # -1.5: need additional y space to account for dist from centerline to top
                    ax.text(i-xAddLR, -majorTickLength - 1.5 - labelOffsetY - (xAddUD+xAddTL),
                            tl[1], ha='center', va='center_baseline', rotation=-90)
    
                elif k in [0.0, 0.5, 1, 2]:
                    # move to bottom
    
                    # draw major x tick
                    ax.plot([i, i],
                            [y0, -(majorTickLength+xAddUD+xAddTL)],
                            "k-", lw=0.8)
    
                    # print ticklabel
                    # -1.5: need additional y space to account for dist from centerline to top
                    ax.text(i, -majorTickLength - 1.5 - labelOffsetY - (xAddUD+xAddTL),
                            tl[1], ha='center', va='center_baseline', rotation=-90)
                    
    
                elif k in [1e-30, 1-1e-1, 1+eps, 1e1]:
                    # move to bottom right
                    
                    # draw major x tick
                    ax.plot([i, i, i+xAddLR, i+xAddLR],
                            [y0, -majorTickLength, -(majorTickLength+xAddUD), -(majorTickLength+xAddUD+xAddTL)],
                            "k-", lw=0.8)
    
                    # print ticklabel
                    # -1.5: need additional y space to account for dist from centerline to top
                    ax.text(i+xAddLR, -majorTickLength - 1.5 - labelOffsetY - (xAddUD+xAddTL),
                            tl[1], ha='center', va='center_baseline', rotation=-90)
    
                elif k in [1+1e-15]:
                    # move up-left by 2
    
                    # draw major x tick
                    ax.plot([i, i, i+2*xAddLR, i+2*xAddLR],
                            [y0, -fac2*majorTickLength, -(majorTickLength+xAddUD+(1-fac2)*xAddTL), -(majorTickLength+xAddUD+xAddTL)],
                            "k-", lw=0.8)
    
                    # print ticklabel
                    ax.text(i+2*xAddLR,
                            -majorTickLength - labelOffsetY - (xAddUD+xAddTL) - 1.5,
                            tl[1], ha='center', va='center_baseline', rotation=-90)
                    
                else:                
                    # draw major x tick
                    ax.plot([i, i], [y0, -majorTickLength], "k-", lw=0.8)
    
                    # print ticklabel
                    # -1.5: need additional y space to account for dist from centerline to top
                    ax.text(i, -majorTickLength - 1.5 - labelOffsetY,
                            tl[1], ha='center', va='center_baseline', rotation=-90)
    
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
    
                if k in [0.0]:
                    # move left
                    
                    # draw major x tick
                    ax.plot([x0, -(majorTickLength+yAddLR+yAddTL)],
                            [i, i],
                            "k-", lw=0.8)
    
                    # print ticklabel
                    ax.text(-majorTickLength - labelOffsetX - (yAddLR+yAddTL),
                            i,
                            tl[1], ha='right', va='center_baseline')
                    
                elif k in [1e-30]:
                    # move up-left by 1
    
                    # draw major x tick
                    ax.plot([x0, -majorTickLength, -(majorTickLength+yAddLR), -(majorTickLength+yAddLR+yAddTL)],
                            [i, i, i+yAddUD, i+yAddUD],
                            "k-", lw=0.8)
    
                    # print ticklabel
                    ax.text(-majorTickLength - labelOffsetX - (yAddLR+yAddTL),
                            i+yAddUD,
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
    plt.subplots_adjust(left=0.05,
                        right=0.86,
                        bottom=0.17,
                        top=0.98)
    
    if savefigFilename is not None:
        plt.savefig(savefigFilename, dpi=300)
    
    plt.show()
