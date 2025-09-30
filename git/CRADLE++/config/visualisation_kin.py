import numpy as np
import matplotlib.pyplot as plt
import sys

def histParam_kin(name):
    var_names = ["px","py","pz"]  
    data_mom = data[[i for i in range(len(names)) if names[i] == name],5:8]
    for j in range(3):
        plt.figure()
        plt.hist(data_mom[:,j], bins=100)
        plt.xlabel(f'{var_names[j]} (keV)')
        plt.ylabel('Counts')
        #plt.yscale('log')
        plt.title(name)
        plt.savefig(name + f'_{var_names[j]}.pdf')
    print(name + " Total # events: " +  str(data_mom.shape[0]))


plt.close('all')
plt.ion()

fileName = sys.argv[1]
try:
    ang_config = sys.argv[2]
    ang_data = False
    if ang_config == "true":
        ang_data = True
except:
    ang_data = False

print(ang_data)

data = np.genfromtxt(fileName, dtype=None)
names = np.array([str(i[2])[2:-1] for i in data])
data = np.array([list(i)[:2] + list(i)[3:] for i in data])
events = data[:,0]
data[np.isnan(data)]=0.

for name in (set(names)&{'e+','e-','enubar','enu'}): #plot only lepton momentum distributions
    histParam_kin(name)

if ang_data: #angular distribution between electron and neutrino
    chargedl_data = data[np.logical_or(names == 'e+',names == 'e-'),5:8]
    neutrino_data = data[np.logical_or(names == 'enu',names == 'enubar'),5:8]
    chargedl_data = chargedl_data/np.tile(np.linalg.norm(chargedl_data,axis=1),(3,1)).T
    neutrino_data = neutrino_data/np.tile(np.linalg.norm(neutrino_data,axis=1),(3,1)).T
    angle = np.sum(chargedl_data*neutrino_data,axis=1)
    plt.figure()
    plt.hist(angle, bins=100)
    plt.xlabel('cos(theta_{e,enu})')
    plt.ylabel('Counts')
    #plt.yscale('log')
    plt.title(name)
    plt.savefig('angular_correlation_enu_e.pdf')

plt.show(block=True)
