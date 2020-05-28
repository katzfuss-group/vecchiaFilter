from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pdb
import os



def readData( family ):

    it = 1
    
    RMSPEfile = os.path.join(family, "RMSPE." + str(it))

    RMSPE = {"MRA" : 0, "LR" : 0}

    count = 0

    while( os.path.exists( RMSPEfile ) ):

        count += 1
        with open(RMSPEfile, "r") as ipFile:
            scores = [l.strip().split(',')[1:] for l in ipFile.readlines()][2:]
            MRA = np.array([float(score[1]) for score in scores])
            LR = np.array([float(score[2]) for score in scores])


        RMSPE["MRA"] = RMSPE["MRA"]*(count-1)/count + MRA/count
        RMSPE["LR"] = RMSPE["LR"]*(count-1)/count + LR/count

        it += 1
        RMSPEfile = os.path.join(family, "RMSPE." + str(it))


    return RMSPE



def plotScore(scoresDict, name):
    
    m = 1e6; M = -1e6
    for family in families:
        score = scoresDict[family]
        m = min(min(score['MRA']), min(score['LR']), m)
        M = max(max(score['MRA']), max(score['LR']), M)


    #print("min = %f, max = %f" % (m, M))

    fig = plt.figure(figsize=(9, 3))
    fig.suptitle("time", x=0.5, y=0.09, fontsize=10)
    for idx, family in enumerate(families):

        if family=="gauss":
            familyName = "Gaussian"
        elif family=="poisson":
            familyName = "Poisson"
        else:
            familyName = family
            
        score = scoresDict[family]
        Tmax = score['MRA'].shape[0]
        time = np.arange(3,Tmax)

        ax = fig.add_subplot(1, 4, idx+1)
        l1, = ax.plot(time, score['MRA'][3:], color="red", linestyle="solid", label="HV")
        l2, = ax.plot(time, score['LR'][3:], color="blue", linestyle=":", label="LR")
        ax.set_title(familyName)
        if(name=="dLS"):
            ax.set_ylim(-0.1*M, 1.05*M)
            l3 = ax.axhline(y=0, color="black", linestyle="dashed")
        elif(name=="RRMSPE"):
            ax.set_ylim(0.95*m, 1.05*M)
            l3 = ax.axhline(y=1.0, color="black", linestyle="dashed")
        elif(name=="RMSPE"):
            ax.set_ylim(0, 1.05*M)

            
        if(idx==0):
            ax.set_ylabel(name)
        else:
            ax.get_yaxis().set_visible(False)


    fig.legend([l1, l2], ["HV", "LR"], "upper center", bbox_to_anchor=[0, 0.02, 1, 1], ncol=2)
    plt.tight_layout(pad=2)
    
    plt.savefig('linear-' + name + '.pdf')  
    
    plt.show()


HOME = Path.home()
os.chdir(os.path.join(HOME, "HVLF/simulations-big/paper_results"))

families = ["gauss", "logistic", "poisson", "gamma"]
RMSPE = {}

for family in families:

    RMSPE[family] = readData(family)

plotScore(RMSPE, "RMSPE")


