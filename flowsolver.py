from math import *
from scipy import optimize
from general import *


class P2C:
    def __init__(self, diam2, diamt, Cd, P1, Pi2, T1, T02, p0, pur1, pur2, R, GAMMA, output=False):
        self.output = output
        self.diam2 = diam2
        self.diamt = diamt
        self.At = pi * self.diamt ** 2 / 4
        self.A2 = pi * self.diam2 ** 2 / 4
        self.Cd = Cd
        self.P1 = P1
        self.Pi2 = Pi2
        self.Pt = P1  # placeholder
        self.Ps2 = Pi2  # placeholder
        self.Pr2 = Pi2  # placeholder
        self.T02 = T02
        self.T1 = T1
        self.p0 = p0
        self.R = R
        self.GAMMA = GAMMA
        self.X1 = P1**G17(self.GAMMA)
        self.Xt = self.X1 # a real value, P. 138
        self.Xi2 = Pi2**G17(self.GAMMA)
        self.Xr2 = 1  # placeholder
        self.Xs2 = 1  # placeholder
        self.T01 = self.T1/(self.X1**2)
        self.d02 = self.p0/(self.R*self.T02)
        self.d0t = self.d02  # from 2.17.3
        self.a02 = sqrt(self.GAMMA*self.R*self.T02)  # 2.17.4
        self.a0t = self.a02  # 2.17.4
        self.d01 = self.p0/(self.R*self.T01) #SHOULD BE T01??
        self.a01 = sqrt(self.GAMMA*self.R*self.T01) #SHOULD BE T01??
        self.dt = self.d0t*self.Xt**G5(self.GAMMA)  # placeholder
        self.Tt = (self.a02*self.X1)**2/(self.GAMMA*self.R)  # placeholder
        self.ct = 0  # placeholder
        self.M = 0  # placeholder
        self.Ts2 = self.T02

    def update(self):
        self.dt = self.d0t*self.Xt**G5(self.GAMMA)  # placeholder
        self.Tt = (self.a02*self.X1)**2/(self.GAMMA*self.R)  # placeholder
        self.Pr2 = self.Xr2**G7(self.GAMMA)
        self.Xs2 = (self.Xi2+self.Xr2-1)
        self.Ps2 = self.Xs2**G7(self.GAMMA)
        self.M = self.ct/(self.a01*self.Xt)
        self.Pt = self.Xt**G7(self.GAMMA)
        self.as2 = self.a02*self.Xs2
        self.Ts2 = (self.as2)**2/(self.GAMMA*self.R)
        self.d2s = (self.d02*self.Xs2**G5(self.GAMMA))
        self.c2s = G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2)
        self.mfr2 = self.A2*self.d2s*self.c2s
        self.M2 = self.c2s/(self.as2)
        self.M = self.ct/(self.a02*self.Xt)
        '''if self.M >= 1.0:
            self.M = 1.0
            self.Xt = self.X1*sqrt(G5(self.GAMMA)/(G5(self.GAMMA)+1))
            self.ct = self.a01*self.Xt'''


    def continuity(self):
        return self.d02*self.X1**G5(self.GAMMA)*self.Cd*self.At*self.ct-self.d02*(self.Xi2+self.Xr2-1)**G5(self.GAMMA)*self.A2*G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2)

    def thd(self):
        first = (G5(self.GAMMA)*(self.a02*(self.Xi2+self.Xr2-1))**2-G5(self.GAMMA)*(self.a02*self.X1)**2)
        second = ((G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2))**2-self.ct**2)
        return first+second

    def equations(self, params):
        self.Xr2, self.ct = params
        self.update()
        return self.continuity(), self.thd()

    def solve(self):
        print("P1: {:.3f}".format(self.P1))
        print("T02: {:.3f}".format(self.T02-273))
        print("dt: {:.3f}".format(self.diamt))
        print("d2: {:.3f}".format(self.diam2))
        print("Cd: {:.3f}".format(self.Cd))
        print("Pi2: {:.3f}".format(self.Pi2))

        print()
        array, infodict, msg, asd = optimize.fsolve(self.equations, (self.Xr2, self.ct), full_output=True)
        if self.output:
            print("Pr2: {:.3f}".format(self.Pr2))
            print("Ps2: {:.3f}".format(self.Ps2))
            print("Ts2: {:.1f}".format(self.Ts2-273))
            print("Pt: {:.2f}".format(self.Pt))
            print("Tt: {:.1f}".format(self.Tt-273))
            print("MFR: {:.1f}".format(self.mfr2*1000))
            print("Ct: {:.1f}".format(self.ct))
            print("Mt: {:.3f}".format(self.M))
            print("c2s: {:.1f}".format(self.c2s))
            print("Ms2: {:.3f}".format(self.M2))
            print("a01: {:.1f}".format(self.a02))

class C2P:
    def __init__(self, diam2, diamt, Cd, P1, Pi2, T1, p0, pur1, pur2, R, GAMMA, output=False):
        self.output = output
        self.diam2 = diam2
        self.diamt = diamt
        self.At = pi*self.diamt**2/4
        self.A2 = pi*self.diam2**2/4
        self.Cd = Cd
        self.P1 = P1
        self.Pi2 = Pi2
        self.Pt = P1  # placeholder
        self.Ps2 = Pi2  # placeholder
        self.Pr2 = Pi2  # placeholder
        self.T1 = T1
        self.p0 = p0
        self.R = R
        self.GAMMA = GAMMA
        self.Cp = (self.GAMMA*self.R)/(self.GAMMA-1)
        self.Cv = self.R/(self.GAMMA-1)
        self.pur1 = pur1
        self.pur2 = pur2

        self.X1 = self.P1**G17(self.GAMMA)
        self.Xi2 = self.Pi2**G17(self.GAMMA)
        self.Xr2 = 1  # placeholder
        self.Xs2 = 1  # placeholder
        self.Xt = self.Pt**G17(self.GAMMA)  # not a real value, see self.Pt
        self.T01 = self.T1/(self.X1**2)
        self.T02 = self.T01  # placeholder

        self.d01 = self.p0/(self.R*self.T01)
        self.d02 = self.d01  # placeholder
        self.d0t = self.d01  # from 2.16.1
        self.dt = self.d01*self.Xt**G5(self.GAMMA)

        self.a01 = sqrt(self.GAMMA*self.R*self.T01)  # from 2.16.2
        self.a0t = self.a01  # from 2.16.2
        self.Tt = (self.a01*self.Xt)**2/(self.GAMMA*self.R)
        self.a02 = self.a01  # placeholder
        self.as2 = self.a02  # placeholder
        self.Ts2 = self.T02  # placeholder
        self.ct = 0  # placeholder

    def update(self):
        self.M = self.ct/(self.a01*self.Xt)
        if self.M >= 1.0:
            self.M = 1.0
            self.Xt = self.X1*sqrt(G5(self.GAMMA)/(G5(self.GAMMA)+1))
            self.ct = self.a01*self.Xt
        self.T02 = (self.a02**2)/(self.GAMMA*self.R)
        self.dt = self.d01*self.Xt**G5(self.GAMMA)
        self.Tt = (self.a01*self.Xt)**2/(self.GAMMA*self.R)
        self.d02 = self.p0/(self.R*self.T02)
        self.Xs2 = self.Xi2+self.Xr2-1
        self.Ps2 = self.Xs2**G7(self.GAMMA)
        self.Pr2 = self.Xr2**G7(self.GAMMA)
        self.as2 = self.a02*self.Xs2
        self.Ts2 = (self.as2)**2/(self.GAMMA*self.R)
        self.Pt = self.Xt**G7(self.GAMMA)
        self.d2s = (self.d02*self.Xs2**G5(self.GAMMA))
        self.c2s = G5(self.GAMMA)*self.a02*(self.Xr2-self.Xi2)
        self.mfr2 = self.A2*self.d2s*self.c2s
        self.M2 = self.c2s/(self.as2)

        #for cyl to pipe conditions, object as mesh
        self.c = self.c2s
        self.a0 = self.a02
        self.d0 = self.d02
        self.T = self.Ts2
        self.Pur = self.pur2

    def continuity(self):
        return self.d01*self.Xt**G5(self.GAMMA)*self.Cd*self.At*-abs(self.ct) - self.d02*(self.Xi2+self.Xr2-1)**G5(self.GAMMA)*self.A2*G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2)

    def thd_1(self):
        return G5(self.GAMMA)*(self.a01*self.X1)**2-((G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2))**2+G5(self.GAMMA)*self.a02**2*(self.Xi2+self.Xr2-1)**2)

    def thd_2(self):
        return G5(self.GAMMA)*((self.a01*self.X1)**2-(self.a01*self.Xt)**2)-self.ct**2

    def momentum(self):
        return self.p0*(self.Xt**G7(self.GAMMA)-(self.Xi2+self.Xr2-1)**G7(self.GAMMA))+(self.d02*(self.Xi2+self.Xr2-1)**G5(self.GAMMA)*G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2))*(-abs(self.ct)-G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2))


    def equations(self, param):
        self.Xr2, self.Xt, self.a02, self.ct = param
        self.update()
        return self.continuity(), self.thd_1(), self.thd_2(), self.momentum()

    def mach_equations(self, param):
        self.Xr2, self.a02 = param
        self.update()
        return self.continuity(), self.thd_1()


    def solve(self):
        array, infodict, msg, asd = optimize.fsolve(self.equations, (self.Xr2, self.Xt, self.a02, self.ct), full_output=True)

        # dirty fix for sonic solution...
        if self.M == 1:
            array, infodict, msg, asd = optimize.fsolve(self.mach_equations, (self.X1, self.a01), full_output=True)
        if self.output:
            print("Pr2: {:.3f}".format(self.Pr2))
            print("Ps2: {:.3f}".format(self.Ps2))
            print("Ts2: {:.3f}".format(self.Ts2-273))
            print("Pt: {:.3f}".format(self.Xt ** G7(self.GAMMA)))
            print("Tt: {:.1f}".format(self.Tt-273))
            print("a02: {:.1f}".format(self.a02))
            print("ct: {:.1f}".format(self.ct))
            print("Mach @ port: {:.1f}".format(self.M))
            print()
            print("cs2: {:.2f}".format(self.c2s))
            print("Mfr2: {:.2f}".format(self.mfr2*1000))
            print("Mach2: {:.3f}".format(self.M2))

class P2P:
    def __init__(self, diam1, diam2, diamt, Cd, P1, P2, T0, T01, p0, R, GAMMA):
        self.diam1 = diam1
        self.diam2 = diam2
        self.diamt = diamt
        self.A1 = pi*diam1**2/4
        self.A2 = pi*diam2**2/4
        self.At = pi*diamt**2/4
        self.Ar = self.A2/self.A1 #area ratio, needed for Benson's constant pressure approximation theory
        self.Cd = Cd
        self.P1 = P1
        self.P2 = P2
        self.Pr1 = P1 #UNKNOWN placeholder, not a real value..
        self.Pr2 = P2 #UNKNOWN placeholder, not a real value..
        self.Pt = self.P1 #UNKNOWN placeholder, not a real value..
        self.R = R
        self.GAMMA = GAMMA
        self.T0 = T0
        self.T01 = T01
        self.T02 = self.T01 #UNKNOWN placeholder, not a real value...
        self.p0 = p0
        self.Xi1 = self.P1**G17(self.GAMMA)
        self.Xi2 = self.P2**G17(self.GAMMA)
        self.a0 = sqrt(self.GAMMA*self.R*self.T0)
        self.a01 = sqrt(self.GAMMA*self.R*self.T01)
        self.d01 = self.p0/(self.R*self.T01)
        self.d0t = self.d01  # from 2.12.1
        self.a0t = self.a01  # from 2.12.2

        #Unknowns
        self.Xr1 = ((1 - self.Ar) * self.Xi1 + 2 * self.Xi2 * self.Ar) / (1 + self.Ar) #Benson's constant p. approx. NOT A REAL VALUE //To be solved 1
        self.Xr2 = (2 * self.Xi1 - self.Xi2 * (1 - self.Ar)) / (1 + self.Ar) #Benson's constant p. approx. NOT A REAL VALUE //To be solved 2
        self.Xt = self.Pt**G17(self.GAMMA) # NOT A REAL VALUE, see self.Pt //To be solved 3
        self.a02 = self.a01 #placeholder, not a real value... //To be solved 4
        self.ct = 0 #placeholder, not a real value... //To be solved 5

        self.c1 = 0 #placeholder, not a real value...
        self.Tt = (self.a01*self.Xt)**2/(self.R*self.GAMMA) #placeholder, not a real value.. (see self.Xt)
        self.d02 = self.d01 #placeholder, not a real value... (see self.T02)
        self.dt = self.d01*self.Xt**G5(self.GAMMA) #placeholder, not a real value... (see self.Xt)
        self.mfr = self.dt*self.At*self.ct #mass flow rate
        self.mfr1 = self.d01*self.A1*self.c1

    def update(self):
        self.Tt = (self.a01*self.Xt)**2/(self.R*self.GAMMA)
        self.T02 = (self.a02)**2/(self.GAMMA*self.R) #solve T02 from 2.12.2
        self.d02 = self.p0/(self.R*self.T02) #from 2.12.1
        self.dt = self.d01*self.Xt**G5(self.GAMMA)
        self.Pt = self.Xt**G7(self.GAMMA) #Pressure ratio P from amplitude ratio X
        self.Pr1 = self.Xr1**G7(self.GAMMA)
        self.Pr2 = self.Xr2**G7(self.GAMMA)

        self.X1s = self.Xi1+self.Xr1-1
        self.X2s = self.Xi2+self.Xr2-1
        self.c1s = G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1)
        self.c2s = G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2)
        self.d1s = self.d01*self.X1s**G5(self.GAMMA)
        self.d2s = self.d02*self.X2s**G5(self.GAMMA)
        self.a1s = self.a01*self.X1s
        self.a2s = self.a02*self.X2s
        self.at = self.a0t*self.Xt
        self.mfr1 = self.d1s*self.A1*self.c1s
        self.mfr2 = self.d2s*self.A2*self.c2s
        self.mfrt = self.d01*self.At*self.ct

        self.M1 = self.c1s/self.a1s
        self.M2 = abs(self.c2s)/self.a2s
        self.Mt = self.ct/self.at


    #Theory 2.10
    def continuity_eq_2107(self):
        first = (self.d01*(self.Xi1+self.Xr1-1)**G5(self.GAMMA)*self.A1*G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1))
        second = (self.d02*(self.Xi2+self.Xr2-1)**G5(self.GAMMA)*self.A2*G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2))
        return first+second
    def first_law_of_thd_2108(self):
        first = ((G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1))**2+G5(self.GAMMA)*self.a01**2*(self.Xi1+self.Xr1-1)**2)
        second = ((G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2))**2+G5(self.GAMMA)*self.a02**2*(self.Xi2+self.Xr2-1)**2)
        return first-second
    def momentum_eq_2109(self):
        first= (self.p0*self.A2*((self.Xi1+self.Xr1-1)**G7(self.GAMMA)-(self.Xi2+self.Xr2-1)**G7(self.GAMMA)))
        second = (self.d01*(self.Xi1+self.Xr1-1)**G5(self.GAMMA)*self.A1*G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1))
        third = (G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1)+G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2))
        return first+second*third

    #Theory 2.11
    def continuity_eq_2117(self):
        return (self.Xi1+self.Xr1-1)**G5(self.GAMMA)*self.A1*(self.Xi1-self.Xr1)+(self.Xi2+self.Xr2-1)**G5(self.GAMMA)*self.A2*(self.Xi2-self.Xr2)
    def first_law_of_thd_2118(self):
        first = (G5(self.GAMMA)*(self.Xi1-self.Xr1)**2+G5(self.GAMMA)*(self.Xi1+self.Xr1-1)**2)
        second = (G5(self.GAMMA)*(self.Xi2-self.Xr2)**2+G5(self.GAMMA)*(self.Xi2+self.Xr2-1)**2)
        return first-second
    def momentum_eq_2109(self):
        first= (self.p0*self.A2*((self.Xi1+self.Xr1-1)**G7(self.GAMMA)-(self.Xi2+self.Xr2-1)**G7(self.GAMMA)))
        second = (self.d01*(self.Xi1+self.Xr1-1)**G5(self.GAMMA)*self.A1*G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1))
        third = (G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1)+G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2))
        return first+second*third

    #Theory 2.12
    def continuity_eq_2128(self):
        first = self.d01*(self.Xi1+self.Xr1-1)**G5(self.GAMMA)*self.A1*self.a01*(self.Xi1-self.Xr1)
        second = self.d02*(self.Xi2+self.Xr2-1)**G5(self.GAMMA)*self.A2*self.a02*(self.Xi2-self.Xr2)
        return first+second
    def continuity_eq_2129(self):
        return (self.Xi1+self.Xr1-1)**G5(self.GAMMA)*self.A1*self.a01*(self.Xi1-self.Xr1)-self.Xt**G5(self.GAMMA)*self.Cd*self.At*self.ct
    def first_law_of_thd_21210(self):
        first = ((G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1))**2+G5(self.GAMMA)*self.a01**2*(self.Xi1+self.Xr1-1)**2)
        second = ((G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2))**2+G5(self.GAMMA)*self.a02**2*(self.Xi2+self.Xr2-1)**2)
        return first-second
    def first_law_of_thd_21211(self):
        return G5(self.GAMMA)*((self.a01*(self.Xi1+self.Xr1-1))**2-(self.a01*self.Xt)**2)+(G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1))**2-self.ct**2
    def momentum_eq_21212(self):
        first = self.p0*self.A2*(self.Xt**G7(self.GAMMA)-(self.Xi2+self.Xr2-1)**G7(self.GAMMA))
        second = (self.d01*(self.Xi1+self.Xr1-1)**G5(self.GAMMA)*self.A1*G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1))*(self.ct+G5(self.GAMMA)*self.a01*(self.Xi2-self.Xr2))
        return first+second

    def continuity_1(self):
        #m1 - m2 = 0
        m1 = self.d01*self.X1s**G5(self.GAMMA)*self.A1*G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1)
        m2 = self.d02*self.X2s**G5(self.GAMMA)*self.A2*G5(self.GAMMA)*self.a02*(self.Xi2-self.Xr2)
        return m1+m2
    def continuity_2(self):
        #m1-mt = 0
        m1 = self.d01*self.X1s**G5(self.GAMMA)*self.A1*G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1)
        mt = self.Cd**2*self.dt*self.At*self.ct
        return m1-mt
    def thd_1(self):
        return (self.c1s**2+G5(self.GAMMA)*self.a1s**2)-(self.c2s**2+G5(self.GAMMA)*self.a2s**2)
    def thd_2(self):
        return (self.c1s**2+G5(self.GAMMA)*self.a1s**2)-(self.ct**2+G5(self.GAMMA)*self.at**2)
    def momentum(self):
        return self.A2*(self.Pt*self.p0-self.X2s**G7(self.GAMMA)*self.p0)+self.d01*self.X1s**G5(self.GAMMA)*self.A1*G5(self.GAMMA)*self.a01*(self.Xi1-self.Xr1)*(self.ct-self.c2s)


    def theory10(self, param):
        self.Xr1, self.Xr2, self.a02 = param
        self.update() #new values for d02, dt, T02 etc...
        if self.M1 > 1:
            self.M1 = 1
            self.Xr1 = (1+G4(self.GAMMA)*self.Xi1)/G6(self.GAMMA)
        return self.continuity_eq_2107(), self.first_law_of_thd_2108(), self.momentum_eq_2109()
    def theory11(self, param):
        self.Xr1, self.Xr2 = param
        self.update() #new values for d02, dt, T02 etc...
        if self.M2 > 1:
            self.M2 = 1
            self.Xr2 = (1+G4(self.GAMMA)*self.Xi2)/G6(self.GAMMA)
        return self.continuity_eq_2117(), self.first_law_of_thd_2118()
    def theory12(self, param):
        self.Xr1, self.Xr2, self.Xt, self.ct, self.a02 = param
        self.update() #new values for d02, dt, T02 etc...
        if self.Mt > 1:
            self.Mt = 1
            self.ct = self.a01*self.Xt
        return self.continuity_1(), self.continuity_2(), self.thd_1(), self.thd_2(), self.momentum_eq_21212()

    def solve(self):
        '''
        if self.A1<self.A2 and self.Cd == 1.0:
            print("T10")
            print()
            array  = optimize.fsolve(self.theory10, (self.Xr1, self.Xr2, self.a01))
        elif self.A1>self.A2 and self.Cd == 1.0:
            print("T11")
            print()
            array = optimize.fsolve(self.theory11, (self.Xr1, self.Xr2))
        else:
            print("T12")
            print()
            array  = optimize.fsolve(self.theory12, (self.Xr1, self.Xr2, self.Xr1, 0, self.a01))
        '''

        array, infodict, msg, asd = optimize.fsolve(self.theory12, (self.Xr1, self.Xr2, self.Xr1, 0, self.a01), full_output=True)

        print("Pr1: {:.4f}".format(self.Pr1))
        print("Pr2: {:.4f}".format(self.Pr2))
        print("T02: {:.1f}".format(self.T02))
        print("MFR1: {:.2f}".format(self.mfr1*1000))
        print("MFRt: {:.2f}".format(self.mfrt*1000))
        print("MFR2: {:.2f}".format(self.mfr2*1000))
        print()
        print("M1: {:.2f}".format(self.M1))
        print("M2: {:.2f}".format(self.M2))
        print("Mt: {:.2f}".format(self.Mt))
        print()
        for val in infodict['fvec']:
            print(val)
        return array

class M2M:
    def __init__(self, T0a, T0b, GAMMAa, GAMMAb, Ra, Rb, X1, X2):
        self.T0a = T0a
        self.T0b = T0b
        self.GAMMAa = GAMMAa
        self.GAMMAb = GAMMAb
        self.Ra = Ra
        self.Rb = Rb
        self.a0a = sqrt(self.GAMMAa*self.Ra*self.T0a)
        self.a0b = sqrt(self.GAMMAb*self.Rb*self.T0b)
        self.X1 = X1
        self.X2 = X2
        self.X1r = 1  # placeholder
        self.X2r = 1  # placeholder

    def momentum(self):
        return (self.a0a*G5(self.GAMMAa))/(self.a0b*G5(self.GAMMAb))*(self.X1-self.X2r)-(self.X1r-self.X2)

    def continuity(self):
        return (self.X1+self.X2r-1)**G7(self.GAMMAa)-(self.X1r+self.X2-1)**G7(self.GAMMAb)

    def equations(self, param):
        self.X1r, self.X2r = param
        return self.continuity(), self.momentum()

    def solve(self):
        array, infodict, msg, asd = optimize.fsolve(self.equations, (1, 1), full_output=True)
        return array