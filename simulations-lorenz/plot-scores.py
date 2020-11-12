import matplotlib.pyplot as plt
import numpy as np
import pdb
import os



def readData( family ):

    it = 1
    
    RRMSPEfile = os.path.join(family, "RRMSPE." + str(it))
    LogScfile  = os.path.join(family, "LogSc."  + str(it))

    RRMSPE = {"MRA" : 0, "LR" : 0}
    LogSc = {"MRA" : 0, "LR" : 0}

    count = 0
    while( os.path.exists( RRMSPEfile ) and it<40):

        count += 1
            
        with open(RRMSPEfile, "r") as ipFile:
            scores = [l.strip().split(',')[1:] for l in ipFile.readlines()][2:]
            MRA = np.array([float(score[1]) for score in scores])
            LR = np.array([float(score[2]) for score in scores])


        RRMSPE["MRA"] = RRMSPE["MRA"]*(count-1)/count + MRA/count
        RRMSPE["LR"] = RRMSPE["LR"]*(count-1)/count + LR/count
        
        with open(LogScfile, "r") as ipFile:
            scores = [l.strip().split(',')[1:] for l in ipFile.readlines()][2:]
            MRA = np.array([float(score[1]) for score in scores])
            LR = np.array([float(score[2]) for score in scores])
            
        LogSc["MRA"] = LogSc["MRA"]*(count-1)/count + MRA/count
        LogSc["LR"] = LogSc["LR"]*(count-1)/count + LR/count
        
        it += 1
        RRMSPEfile = os.path.join(family, "RRMSPE." + str(it))
        LogScfile  = os.path.join(family, "LogSc."  + str(it))


    return RRMSPE, LogSc




def plotScore(scoresDict, name):


    m = 1e6; M = -1e6
    for family in families:
        score = scoresDict[family]
        m = min(min(score['MRA']), min(score['LR']), m)
        M = max(max(score['MRA']), max(score['LR']), M)


    

    fig = plt.figure(figsize=(9, 3))
    if max(abs(scoresDict['gauss']['MRA']))>10:
        fig.suptitle("time", x=0.5, y=0.09, fontsize=10)
    #fig.suptitle("N, size of the conditioning set", x=0.5, y=0.09, fontsize=10)
    for idx, family in enumerate(families):

        if family=="gauss":
            familyName = "Gaussian"
        elif family=="poisson":
            familyName = "Poisson"
        else:
            familyName = family
            
        score = scoresDict[family]
        Tmax = score['MRA'].shape[0]
        time = np.arange(1,Tmax+1)
        ax = fig.add_subplot(1, 4, idx+1)
        l1, = ax.plot(time, score['MRA'], color="red", linestyle="solid", label="HV")
        l2, = ax.plot(time, score['LR'], color="blue", linestyle=":", label="low-rank")
        ax.set_title(familyName)
        if(name=="dLS"):
            ax.set_ylim(-0.1*M, 1.05*M)
            l3 = ax.axhline(y=0, color="black", linestyle="dashed")
        elif(name=="RRMSPE"):
            ax.set_ylim(1-0.05*M, 1.05*M)
            l3 = ax.axhline(y=1.0, color="black", linestyle="dashed")

        if(idx==0):
            ax.set_ylabel(name)
        else:
            ax.get_yaxis().set_visible(False)
        print(max(abs(scoresDict['gauss']['MRA'])))
        if max(abs(scoresDict['gauss']['MRA']))<5:
            ax.get_xaxis().set_visible(False)
            ax.set_title(familyName)

            
    if max(abs(scoresDict['gauss']['MRA']))<10:
        fig.legend([l1, l2, l3], ["HV", "LR", "DL"], "upper center", bbox_to_anchor = [0, 0.02, 1, 1], ncol=3)
    plt.tight_layout(pad=2)

    #fig.legend([l1, l2], labels=["HV", "low-rank", "Laplace"], ncol=3, bbox_to_anchor=(-0.3, -0.022, 1, 1))
    #plt.tight_layout(pad=2)

    plt.savefig('lorenz-' + name + '.pdf')  
    
    plt.show()



os.chdir("/home/marcin/vecchiaFilter/simulations-lorenz")

families = ["gauss", "logistic", "poisson", "gamma"]
RRMSPE = {}
LogSc = {}

for family in families:

    RRMSPE[family], LogSc[family] = readData(family)

plotScore(RRMSPE, "RRMSPE")
plotScore(LogSc, "dLS")

