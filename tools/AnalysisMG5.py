import numpy as np
import pylhe
import sys


def invariant_mass(p1, p2):
    return np.sqrt(
        sum(
            (1 if mu == "e" else -1) * (getattr(p1, mu) + getattr(p2, mu)) ** 2
            for mu in ["e", "px", "py", "pz"]
        )
    )


def CalculoGB(particle):
    modMoment = np.sqrt(particle.px ** 2 + particle.py ** 2 + particle.pz ** 2)
    V = modMoment / particle.e
    g = V / np.sqrt(1 - (V**2))
    Vg = g / np.sqrt(1 + (g**2))
    M = particle.m
    return g

def angle(Moment=[]):

    r = Moment[2] / (np.sqrt(Moment[2]** 2 + Moment[1] ** 2 + Moment[0]** 2))
    #r = np.cos(r)
    return r


def modmome(particle):

    mome = (np.sqrt(particle.px ** 2 + particle.py ** 2 + particle.pz ** 2))

    return mome


def boost(V, p):
    """
    :param V: frame rate
    :param p: 4-moment of the particle being boosted

    :return: 4-moment in the S' frame
    """

    u = np.dot(V, V)
    v4 = [1, V[0], V[1], V[2]]

    if u == 0:
        return p

    if u > 1.0:
        print('>It is not possible to calculate the momentum for this passed velocity<')
        sys.exit()

    g = 1 / np.sqrt(1 - u)

    # Lorentz transformation matrix

    Ce = np.einsum("i,j->ij", V, V) * ((g + 1) / u) + np.identity(3)
    L = [v * g for v in V]
    C = [v * g for v in v4]
    D = np.insert(Ce, 0, L, axis=0)
    M = np.insert(D, 0, C, axis=1)

    P = np.dot(p, M)

    return P





def Analysis(filename, mae, filha1, filha2):

    G35 = []
    G36 = []
    G37 = []
    P21 = []
    P1 = []
    P2 = []
    P3 = []
    P4 = []
    Ang = []
    EF1 = []
    MoM = []
    VMOM =[]
    MOMMass = []
    EF2 = []
    Vfilha1 = []
    VMOM3 = []
    VFilha3 = []
    Massa2 = []
    Ang2 = []
    E2 = []
    PP2 = []
    MI = []
    E2F = []
    P2F = []

    for e in pylhe.readLHE(filename):
        for Mom in e.particles:
            V=0
            if abs(Mom.id) == 21:
                P21.append(Mom.id)
            if abs(Mom.id) == 1:
                P1.append(Mom.id)
            if abs(Mom.id) == 2:
                P2.append(Mom.id)
                if Mom.mother1 == 4:
                    E = Mom.e
                    E2.append(E)
                    P = [Mom.px, Mom.py, Mom.pz]
                    PP2.append(P)


            if abs(Mom.id) == 3:
                P3.append(Mom.id)
            if abs(Mom.id) == 4:
                P4.append(Mom.id)
            if abs(Mom.id) == mae:
                g = CalculoGB(Mom)
                G37.append(g)
                Vx = Mom.px/Mom.e
                Vy = Mom.py/Mom.e
                Vz = Mom.pz/Mom.e
                VMOM.append(np.sqrt(Vx**2+Vz**2+Vy**2))
                VMOM3.append([-Vx,-Vy,-Vz])
                MOMMass.append(Mom.m)
            if abs(Mom.id) == filha1:
                if Mom.mother1 == 3:
                    MTMomento = boost(VMOM3[-2], [Mom.e, Mom.px, Mom.py, Mom.pz])
                    EF1.append(MTMomento[0])
                    Ang.append(angle([MTMomento[1], MTMomento[2], MTMomento[3]]))
                else:
                    MTMomento = boost(VMOM3[-1],[Mom.e,Mom.px,Mom.py,Mom.pz])
                    EF1.append(MTMomento[0])
                    Ang.append(angle([MTMomento[1],MTMomento[2],MTMomento[3]]))
                g = CalculoGB(Mom)
                G36.append(g)
                Vx1 = Mom.px/Mom.e
                Vy1 = Mom.py/Mom.e
                Vz1 = Mom.pz/Mom.e
                VFilha3.append([-Vx1, -Vy1, -Vz1])
                Massa2.append(np.sqrt((Mom.e ** 2) - modmome(Mom) ** 2))
            if abs(Mom.id) == filha2:
                g = CalculoGB(Mom)
                G35.append(g)
                if Mom.mother1 == 4:
                    MT2Momento = boost(VFilha3[-2],[Mom.e,Mom.px,Mom.py,Mom.pz])
                    EF2.append(MT2Momento[0])
                    Ang2.append(angle([MT2Momento[1], MT2Momento[2], MT2Momento[3]]))
                else:
                    MT2Momento = boost(VFilha3[-1],[Mom.e,Mom.px,Mom.py,Mom.pz])
                    EF2.append(MT2Momento[0])
                    Ang2.append(angle([MT2Momento[1], MT2Momento[2], MT2Momento[3]]))




    return G37,G36,G35,P1,P2,P21,P3,P4,Ang,EF1,MoM,VMOM,MOMMass,EF2,Ang2,Massa2,E2,PP2

