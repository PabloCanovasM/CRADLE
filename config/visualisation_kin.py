import numpy as np
import matplotlib.pyplot as plt
import sys

def histParam_kin(p_names):
    for name in p_names:
        var_names = ["px","py","pz"]
        data_mom = data[names == name,5:]
        for j, var_name in enumerate(var_names):
            plt.figure()
            plt.hist(data_mom[:,j], bins=100)
            plt.xlabel(f'{var_name} (keV)')
            plt.ylabel('Counts')
            #plt.yscale('log')
            plt.title(name)
            plt.savefig(f'{name}_{var_name}.pdf')
        print(f"{name} Total # events: {data_mom.shape[0]}")

def histParam_kin2(p_names):
    for name in p_names:
        var_names = ["px","py","pz"]
        data_mom = data[names == name,5:8]
        data2_mom = data2[names2 == name,5:8]
        for j, var_name in enumerate(var_names):
            plt.figure()
            plt.hist(data_mom[:,j], bins=100,histtype="step",label=fileName[:-4])
            plt.hist(data2_mom[:,j], bins=100,histtype="step",label=fileName2[:-4])
            plt.xlabel(f'{var_name} (keV)')
            plt.ylabel('Counts')
            #plt.yscale('log')
            plt.title(name)
            plt.legend()
            plt.savefig(f'{name}_{var_name}.pdf')
        print(f"{name} Total # events 1st file: {data_mom.shape[0]}")
        print(f"{name} Total # events 2nd file: {data2_mom.shape[0]}")        
        data_E = data[names == name,4]
        data_v = data_mom/data_E[:,np.newaxis]        
        data2_E = data2[names2 == name,4]
        data2_v = data2_mom/data2_E[:,np.newaxis]        
        var_names2 = ["vx","vy","vz"]
        for j, var_name in enumerate(var_names2):
            plt.figure()
            plt.hist(data_v[:,j], bins=100,histtype="step",label=fileName[:-4])
            plt.hist(data2_v[:,j], bins=100,histtype="step",label=fileName2[:-4])
            plt.xlabel(f'{var_name} (keV)')
            plt.ylabel('Counts')
            #plt.yscale('log')
            plt.title(name)
            plt.legend()
            plt.savefig(f'{name}_{var_name}.pdf')
        print(f"{name} Total # events 1st file: {data_mom.shape[0]}")
        print(f"{name} Total # events 2nd file: {data2_mom.shape[0]}")



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

try:
    fileName2 = sys.argv[3]
except:
    fileName2 = None

data = np.genfromtxt(fileName, dtype=str)
names = data[:,2]
data = data[:, np.arange(data.shape[1]) != 2].astype(dtype=np.float32)
events = data[:,0]
data[np.isnan(data)]=0.

if fileName2 != None:
    data2 = np.genfromtxt(fileName2, dtype=str)
    names2 = data2[:,2]
    data2 = data2[:, np.arange(data2.shape[1]) != 2].astype(dtype=np.float32)
    events2 = data2[:,0]
    data2[np.isnan(data2)]=0.

if fileName2 == None:
    p_names = (set(names)&{'e+','e-','enubar','enu'})
    histParam_kin(p_names)
else:
    if set(names) == set(names2):
        p_names = (set(names)&{'e+','e-','enubar','enu'})
        histParam_kin2(p_names)
    else:
        print("Comparing decays with different particles. Aborting")

if ang_data: #angular distribution between electron and neutrino
    chargedl_data = data[np.logical_or(names == 'e+',names == 'e-'),5:8]
    neutrino_data = data[np.logical_or(names == 'enu',names == 'enubar'),5:8]
    chargedl_data = chargedl_data/np.tile(np.linalg.norm(chargedl_data,axis=1),(3,1)).T
    neutrino_data = neutrino_data/np.tile(np.linalg.norm(neutrino_data,axis=1),(3,1)).T
    angle = np.sum(chargedl_data*neutrino_data,axis=1)
    if fileName2 != None:
        chargedl_data2 = data2[np.logical_or(names2 == 'e+',names2 == 'e-'),5:8]
        neutrino_data2 = data2[np.logical_or(names2 == 'enu',names2 == 'enubar'),5:8]
        chargedl_data2 = chargedl_data2/np.tile(np.linalg.norm(chargedl_data2,axis=1),(3,1)).T
        neutrino_data2 = neutrino_data2/np.tile(np.linalg.norm(neutrino_data2,axis=1),(3,1)).T
        angle2 = np.sum(chargedl_data2*neutrino_data2,axis=1)
        plt.figure()
        plt.hist(angle, bins=100, histtype="step",label=fileName[:-4])
        plt.hist(angle2, bins=100, histtype="step",label=fileName2[:-4])
        plt.xlabel('cos(theta_{e,enu})')
        plt.ylabel('Counts')
        #plt.yscale('log')
        plt.legend()
        plt.savefig('angular_correlation_enu_e.pdf')
    else:
        plt.figure()
        plt.hist(angle, bins=100)
        plt.xlabel('cos(theta_{e,enu})')
        plt.ylabel('Counts')
        #plt.yscale('log')
        plt.savefig('angular_correlation_enu_e.pdf')

plt.show(block=True)
