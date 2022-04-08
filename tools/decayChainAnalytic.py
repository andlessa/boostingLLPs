import numpy as np
import sys
from scipy.interpolate import interp1d

def getAvgGammaBeta(Ms,beta0Process,computeStd=bool,V=[]):
    """
    Method to return the average beta gamma of a series of decay using an analytical approximation.
    ....
    :param V: velocity vector for first parent (e.g. [vx,vy,vz]).
    :param Ms: list with masses for all particles in the decay chain, with each decay being a tuple with the first entry being the particle that will decay (M)
               and the second its companion that will not decay ([(M0,m0),(M1,m1),...,(MN,mn]).The masses must be sorted in descending order.
    :param ComputerStd: boolean variable to define whether or not to calculate the standard deviation.

    :return: List with average gamma*beta values with their respective standard deviations for each step of the cascade decay
             ([(gb0,std0),(gb1,std1),...(gbN,stdN)]).
    """
    if len(V) == 0:
        beta = velocity(Ms[0][0],beta0Process)
    else:
        beta = np.sqrt(np.dot(V,V))

    gamma = 1/np.sqrt(1-beta**2)
    GammaBeta = [(beta*gamma,0)] #tuple with the gammabeta and the uncertainty of the mother particle
    Vmom = [0,0,beta]

    for i, Mmom in enumerate(Ms[:-1]):
        Mdaughter1 = Ms[i + 1][0]
        Mdaughter2 = Ms[i + 1][1]
        Mmom = Mmom[0]
        gb = getDecayGammaBeta(MX=Mmom, M1=Mdaughter1, M2=Mdaughter2, VX=Vmom)
        if computeStd:
            Dsv = getDecayStdDev(MX=Mmom, M1=Mdaughter1, M2=Mdaughter2,
                                 VX=Vmom, gbAvg=gb)
        else:
            Dsv = 0
        beta = gb/np.sqrt(1+gb**2) #Compute beta from gamma*beta
        Vmom = [0, 0, beta] #Define parent velocity for the next decay
        GammaBeta.append((gb,Dsv))

    return np.array(GammaBeta)



def bintres(x, a):
    # Result of integral in theta
    return (x*np.sqrt(abs(x ** 2 - a**2)))/(2*a**2) - (1/2)*(np.log((x+np.sqrt(abs(x**2-a**2)))/a))

def getDecayGammaBeta(MX, M1, M2, VX=[]):

    if M1 > MX or M2 > MX:
        print('>It is not possible to calculate the momentum for these masses<')
        sys.exit()
    else:
        beta = np.sqrt(np.dot(VX, VX))
        g = 1 / np.sqrt(1 - beta ** 2)
        E1cm = (MX ** 2 + M1 ** 2 - M2 ** 2) / (2 * MX) # calculate the energy of the decaying particle at the center of mass of the mother particle
        if ((E1cm ** 2) - (M1 ** 2)) < 0:
            return beta/np.sqrt(1-(beta**2))
        else:
            P1cm = np.sqrt(E1cm ** 2 - M1 ** 2) # calculate the moment of the decaying particle at the center of mass of the mother particle
            M1t = M1 / g
            if beta == 0:
                V = (P1cm / E1cm)/np.sqrt(1-(P1cm / E1cm)**2)
                return V
            else:
                if P1cm == 0:
                    return beta/np.sqrt(1-(beta**2))
                else:
                    bV1 = (M1t/(2 * beta * P1cm)) * (bintres(E1cm + beta * P1cm, M1t) - bintres(E1cm - beta * P1cm, M1t))

                    return np.around(bV1,4)





def getDecayStdDev(MX, M1, M2, VX=[], gbAvg=None):
    if M1 > MX or M2 > MX:
        print('>It is not possible to calculate the momentum for these masses<')
        sys.exit()


    beta = np.sqrt(np.dot(VX, VX))
    g = 1 / np.sqrt(1 - beta ** 2)
    E1cm = (MX ** 2 + M1 ** 2 - M2 ** 2) / (2 * MX)  # calculate the energy of the decaying particle at the center of mass of the mother particle
    M1t = M1/g
    if ((E1cm ** 2) - (M1 ** 2)) < 0:
        return 0
    else:
        P1cm = np.sqrt(abs(E1cm ** 2 - M1 ** 2)) # calculate the moment of the decaying particle at the center of mass of the mother particle
        if P1cm == 0 or beta == 0:
            return 0

        if gbAvg is None:
            gbAvg = getDecayGammaBeta(MX, M1, M2, VX)

        Ds1b = np.sqrt((((1/(2*beta*P1cm))*((((E1cm+beta*P1cm)**3/(3*M1t**2))-(E1cm+beta*P1cm))-(((E1cm-beta*P1cm)**3/(3*M1t**2))-(E1cm-beta*P1cm))))-(gbAvg ** 2)))

        return np.around(Ds1b,4)

def velocity (M,beta0Process):

    #Function to set a velocity for the mother particle from its mass.

    if beta0Process == 'VH':
        V=-8.57144225e-12*(M**3)+7.63775328e-8*(M**2)-0.000308996381*M+1.00495090
    if beta0Process == 'VGO':
        V = -1.01378244e-11 * (M**3) + 9.33615192e-8 * (M**2) - 0.00035919197 * M + 9.4901006e-1
    if beta0Process == 'VER':
        V  = -8.46530365966e-12 * (M**3) + 7.587104215644e-8 * (M**2) - 0.000309309092734 * M + 1.00793903613707
    if beta0Process == 'VX':
        V = -1.0875078579e-11 * (M**3) + 1.010961291543e-7 * (M**2) - 0.000388407876673 * M + 1.00278035034909
    if beta0Process == 'VUR':
        V = -1.11705558773826e-11 * (M**3) + 1.02524960991756e-7 * (M**2) - 0.000357938654903 * M + 0.938547539837664
    if beta0Process == 'VANA':
        V = -1.10326218522e-11 * (M**3) + 9.915902460305e-8 * (M**2) - 0.000367227071633 * M + 1.00098600163719

    return V
