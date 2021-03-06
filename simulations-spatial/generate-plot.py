import matplotlib.pyplot as plt
import numpy as np
import pdb
import pandas as pd


data = pd.read_csv("results").filter(items=["Mod", "Neighbors", "Seed_off", "MSE_Laplace", "MSE_VL", "MSE_LowRank", "LS_Laplace", "LS_VL", "LS_LowRank"])


MSE_L = np.zeros(len(data))
dLS_L = np.zeros(len(data))
MSE = 0
dLS = 0

for row_id in range(len(data)):
    if data.MSE_Laplace[row_id] != 0:
        MSE = data.MSE_Laplace[row_id]
        dLS = data.LS_Laplace[row_id]

    MSE_L[row_id] = MSE
    dLS_L[row_id] = dLS

    
data["MSE_Laplace"] = MSE_L
data["LS_Laplace"] = dLS_L

families = ["Gaussian", "logistic", "Poisson", "gamma"]

data["RRMSPE_VL"] =  data["MSE_VL"].div(data.MSE_Laplace, axis=0)
data["RRMSPE_LR"] =  data["MSE_LowRank"].div(data.MSE_Laplace, axis=0)
data["dLS_VL"]    = -data["LS_VL"].sub(data.LS_Laplace, axis=0)
data["dLS_LR"]    = -data["LS_LowRank"].sub(data.LS_Laplace, axis=0)

RRMSPEs = data.filter(["RRMSPE_VL", "RRMSPE_LR", "Neighbors", "Mod"]).groupby(["Mod", "Neighbors"]).mean()
dLSs = data.filter(["dLS_VL", "dLS_LR", "Neighbors", "Mod"]).groupby(["Mod", "Neighbors"]).mean()


fig = plt.figure(figsize=(9, 3))
#fig.suptitle("N, size of the conditioning set", x=0.5, y=0.09, fontsize=10)
#plt.text(0.5, 0.13, 'matplotlib', horizontalalignment='center',verticalalignment='center', transform=fig.transFigure.transform())
#plt.text(0.5, 0.13, 'matplotlib', transform=fig.transFigure)
for idx, family in enumerate(families):

    scores = RRMSPEs[RRMSPEs.index.get_level_values('Mod') == idx+1]
    Neighbors = np.array(scores.index.get_level_values('Neighbors'))

    
    
    ax = fig.add_subplot(1, 4, idx+1)
    l1, = ax.plot(Neighbors, scores['RRMSPE_VL'], color="red", linestyle="solid", label="HV")
    l2, = ax.plot(Neighbors, scores['RRMSPE_LR'], color="blue", linestyle=":", label="low-rank")
    ax.set_title(family)
    l3 = ax.axhline(y=1.0, color="black", linestyle="dashed")
    ax.set_ylim(1-0.1*max(RRMSPEs.max()), 1.05*max(RRMSPEs.max()))
    ax.get_xaxis().set_visible(False)
    
    if(idx==0):
        ax.set_ylabel("RRMSPE")
    else:
        ax.get_yaxis().set_visible(False)



#fig.legend([l1, l2, l3], labels=["HV", "LR", "DL"], ncol=3, bbox_to_anchor=(-0.3, 0.022, 1, 1))
fig.legend([l1, l2, l3], ["HV", "LR", "DL"], "upper center", ncol=3)
plt.tight_layout(pad=2.5)

plt.savefig('spatial-RRMSPE.pdf')  

plt.show()



fig = plt.figure(figsize=(9, 3))
fig.suptitle("N, size of the conditioning set", x=0.5, y=0.09, fontsize=10)
for idx, family in enumerate(families):


    scores = dLSs[dLSs.index.get_level_values('Mod') == idx+1]
    Neighbors = np.array(scores.index.get_level_values('Neighbors'))

    ax = fig.add_subplot(1, 4, idx+1)
    l1, = ax.plot(Neighbors, scores['dLS_VL'], color="red", linestyle="solid", label="HV")
    l2, = ax.plot(Neighbors, scores['dLS_LR'], color="blue", linestyle=":", label="low-rank")
    #ax.set_title(family)
    l3 = ax.axhline(y=0, color="black", linestyle="dashed")

    
    m = min(dLSs.min())
    M = max(dLSs.max())
    ax.set_ylim(-0.1*M, 1.05*M)
    
    if(idx==0):
        ax.set_ylabel("dLS")
    else:
        ax.get_yaxis().set_visible(False)

#fig.legend([l1, l2, l3], labels=["HV", "LR", "DL"], ncol=3, bbox_to_anchor=(-0.3, 0.022, 1, 1), loc="upper center")
#fig.legend([l1, l2, l3], ["HV", "LR", "DL"], loc="upper center", ncol=3)
plt.tight_layout(pad=2.5)

plt.savefig('spatial-dLS.pdf')  


        
plt.show()


