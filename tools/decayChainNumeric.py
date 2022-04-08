import numpy as np
from EventGenerator import generator
import json

"""
module responsible for carrying out a chain of decays, returning a list of gammabeta, 
returns the average gametes with their uncertainties and can create a file with all the gammabetas in the decay chain.

author: Felipe Zanzin.

"""



def NDecay(V, Ms):
    """
    Module with function to return the 4-moment vector of a decay directly
    ....
    :param V: speeds List (3-vector)
    :param Ms: masses List
    """

    Vmom = V
    moments = []

    for i, Mmom in enumerate(Ms[:-1]):
        Mdaughter1 = Ms[i + 1]
        Mdaughter2 = 0.0
        events = generator(mx=Mmom, m1=Mdaughter1, m2=Mdaughter2, v=[Vmom], n=1).ev
        event = events[0]
        daughter1 = event.p1  # 4-moment daugther1
        moments.append(daughter1)
        Vmom = list(daughter1.mspeeds)
    return moments


def calcAvgDsvNumgb(n, V, Ms, outputFile=None):
    """
    Module with function to directly return the average beta gamma and its standard deviation of a decay series
    ....
    :param V: speeds List (3-vector)
    :param Ms: masses List
    """

    GammaB = []

    for j in range(0, n, 1):
        momentos = NDecay(V, Ms)
        gammaBeta = [round(np.sqrt(np.dot(V, V)) / (np.sqrt(1 - (np.dot(V, V)))),4)]
        #Calculating the beta gamma for a Event
        for p in momentos:
            gb = round(p.beta * (1 / (np.sqrt(1 - (p.beta ** 2)))),4)
            gammaBeta.append(gb)
        GammaB.append(gammaBeta)


    GammaBa=np.array(GammaB)
    GammaBaverage=[(np.around(np.sqrt(np.dot(V, V)) / (np.sqrt(1 - (np.dot(V, V)))),4),0)]
    for i in range (1,len(GammaBa[0]),1):
        gbDecay = np.around(np.mean(GammaBa[:,i][:]),4)
        DSVDecay = np.around(np.std(GammaBa[:,i][:]),4)
        GammaBaverage.append((gbDecay,DSVDecay)) #Calculating the average beta gamma and the gamma beta standard deviation of each decay

    if outputFile != None :
    #file with essential decay information and all gammabetas generated in the decay chain.
        file = open(outputFile, 'w')
        file.write("#gammabeta (%i particles, %i decays), first line masses List of the entire decay chain, second line the velocity of the mother particle, third line number of events \n" %(len(Ms),len(Ms)-1))
        file.write("#The gammabetas are distributed as follows (mother's gammabeta, Decay gammabetas)" + '\n')
        file.write(json.dumps(Ms)+'\n'+json.dumps(V)+'\n'+str(n)+'\n')

        for gamma in GammaB:
            file.write(json.dumps(gamma)+'\n')

        file.close()

    return GammaBa


def fromFile(arquivoEntrada):
#function to read files originating from calcAvgDsvNumgb.

    evList = []
    arquivo = open(arquivoEntrada, 'r')
    data = arquivo.read()
    decayBlocks = data.split("[")
    decayBlocks = [d.split("]") for d in decayBlocks if d.strip()] #GammaBetas for decay chain.
    for i, b in enumerate(decayBlocks):
        evList.append(Evento(blockStr=b))
    return evList


def Evento(blockStr):
    lines = blockStr.strip().split("\n")
    gammabl = []
    for l in lines[:]:
        gammab = [eval(x) for x in l.split(",")]
        gammabl.append(gammab)
    return gammabl
