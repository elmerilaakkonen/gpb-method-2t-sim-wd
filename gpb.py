from math import *
from general import *
import scipy.optimize as optimize
import numpy as np
np.set_printoptions(precision=3)
np.set_printoptions(suppress=False)



class Duct:
    def __init__(self, d1, d2, rho0, T0, p0, a0, R, G, L, n, inHandle, outHandle):
        self.d1 = d1/1000       # diameter start
        self.d2 = d2/1000       # diameter end
        self.rho0 = rho0        # reference density
        self.T0 = T0            # reference temperature KELVIN
        self.p0 = p0            # reference pressure
        self.a0 = a0            # reference acoustic velocity
        self.R = R              # reference gas constant (R)
        self.G = G              # reference gamma for gas
        self.L = L/1000         # total length of duct (in meters)

        self.dt = 0             # init dt

        self.inHandle = inHandle
        self.outHandle = outHandle

        self.meshes = self.meshify(n)
        self.calculate_dt()

    def meshify(self, n):
        meshes = []
        mesh_len = self.L / n
        delta_d = (self.d2-self.d1)/(n-1)   # change in diameter between meshes


        for i in range(n):
            d = self.d1+delta_d*i           # mesh diameter = d1 + change * mesh_num
            meshes.append(Mesh(d, mesh_len, self.rho0, self.T0, self.p0, self.a0, self.R, self.G))

        # Default open ended duct
        meshes[0].first = 2
        meshes[-1].last = 2

        return meshes

    def print_info(self):
        print()
        print("    {:6} | {:6} | {:6} | {:6}".format("XR", "XL", "XR1", "XL1"))
        for i in range(len(self.meshes)):
            mesh = self.meshes[i]
            print("{:2}. {:.4f} | {:.4f} | {:.4f} | {:.4f}".format(i, mesh.XR, mesh.XL, mesh.XR1, mesh.XL1))
        print()

    def calculate_dt(self):
        dt = self.meshes[0].dt
        for mesh in self.meshes:
            if dt == 0 or mesh.dt < dt: dt = mesh.dt
        self.dt = dt

    def step(self):
        for i in range(0, len(self.meshes)):
            mesh = self.meshes[i]
            mesh.friction(self.dt)

        for i in range(0, len(self.meshes)-1):
            mesh = self.meshes[i]
            mesh.boundaries(self.meshes[i+1])
            #exit()

        for mesh in self.meshes:
            mesh.update()

        self.calculate_dt()

class Port:
    def __init__(self, d, Cd, R=0, height=0, width=0, fromtdc=0, cyl=None, eng=None, type='exhaust'):
        self.A = pi*(d/1000)**2/4
        self.Cd = Cd
        self.d = d/1000
        self.Mach = False
        self.openArea = 0
        self.R = R / 1000
        self.height = height / 1000
        self.width = width / 1000
        self.fromtdc = fromtdc / 1000  # port top from TDC in mm (if measured from barrel top, take in account the deck height)
        self.cyl = cyl
        self.eng = eng
        self.type = type

        if cyl:
            self.fromtdc -= cyl.deck / 1000

    def open_deg(self):
        eng = self.eng
        l = eng.L_cr
        r = eng.L_st / 2
        x = r + l - self.fromtdc
        A = -1 / 2 * (l**2 - r**2 - x**2) / (r * x)
        return degrees(acos(A))

    def close_deg(self):
        return self.open_deg()+self.duration()

    def open_width(self, x):

        # TODO: FIX, Bottom width not the same as top (w/ same radius)??
        # TODO: for intake port also

        atdc = self.fromtdc                # port roof below tdc in mm
        h = self.height
        R = self.R
        w = self.width

        if self.type != 'intake':
            if x < atdc:
                # port not open, no width
                width = 0
            elif atdc <= x <= atdc+R:
                # port just opening, in the area with corner radius
                width = w - 2*R + 2*np.sqrt(R**2 - round(atdc+R-x, 9)**2)
            elif atdc+h > x > atdc+h-R:
                width = w - 2*R + 2*np.sqrt(R**2 - round(atdc+h-R-x, 9)**2)
            elif atdc+h > x:
                width = w
            else:
                width = 0  # fallback, returning 0 even if port is 100% open for correct area calculation
        else:
            width = 0

        return width

    def area(self, theta):
        eng = self.eng
        od = self.open_deg()  # opening degrees


        # Area is sum of slices with different widths and even height of 0.01mm

        open_x = self.fromtdc
        to_x = eng.piston_position(theta) / 1000
        open_height = to_x - open_x
        A = 0

        if open_height > 0:
            N = int(open_height / (0.01 / 1000))
            slice_height = open_height / N

            A = 0
            for i in range(N+1):
                cur_x = round(open_x + i * slice_height, 9)
                cur_w = self.open_width(cur_x)
                A += cur_w*slice_height

        return A

    def update(self):
        eng = self.eng
        self.A = self.area(eng.angle)

    def duration(self):
        cyl = self.cyl
        eng = self.eng
        stroke = eng.L_st
        conrod = eng.L_cr

        if self.type != 'intake':
            R = stroke / 2
            L = conrod
            C = cyl.deck
            E = self.fromtdc
            T = R + L + C - E
            duration = (180-degrees(acos((T**2+R**2-L**2)/(2*R*T))))*2
        else:
            R = stroke / 2
            L = conrod
            C = cyl.deck
            H = eng.piston_intake_h
            F = self.fromtdc + self.height
            P = R + L + H + C - F
            duration = degrees(acos((P**2+R**2-L**2)/(2*P*R)))*2


        return duration



class Plenum:
    def __init__(self, rho0, T0, p0, a0, R, G, T, p, ports = False, volume = 0, deck = 0):
        self.rho0 = rho0  # reference density of gas in the mesh
        self.T0 = T0  # reference temperature
        self.p0 = p0  # reference pressure
        self.a0 = a0  # reference acoustic velocity
        self.R = R  # gas constant
        self.G = G  # gamma
        self.T = T  # gas temperature
        self.p = p  # pressure in the plenum

        self.eng = None
        self.deck = deck / 1000
        self.volume = volume

        if ports != False:
            self.ports = ports
            for portX in self.ports:
                self.ports[portX].cyl = self
                self.ports[portX].cyl.type = portX

    def boundaries(self, mesh, closed=False):

        if closed:
            mesh.first = 1
            mesh.XR = mesh.Xqf
        else:
            P_cyl = self.p / self.p0  # pressure ratio of the gas in the cylinder
            P_pipe = mesh.Xqf**G7(mesh.G)  # pressure ratio of the gas in the pipe / duct (Pi2 in page. 141 tables)


            if P_cyl > P_pipe:
                mesh.first = 2
                # flowing to pipe
                exhaustHandle = flowSolver(self, mesh, self.ports['exhaust'])
                Xr2, Xt, a02 = exhaustHandle.cylinder_outflow()
                mesh.XR = Xr2
            else:
                # flowing to cyl
                exhaustHandle = flowSolver(self, mesh, self.ports['exhaust'])
                x = exhaustHandle.cylinder_inflow()


class Mesh:
    def __init__(self, d, L, rho0, T0, p0, a0, R, G, XR=1, XL=1, XR1=1, XL1=1, Cd = 1.0, dth = False):
        self.d = d              # diameter
        self.L = L              # length of mesh
        self.rho0 = rho0        # reference density of gas in the mesh
        self.T0 = T0            # reference temperature
        self.p0 = p0            # reference pressure
        self.a0 = a0            # reference acoustic velocity
        self.R = R              # gas constant
        self.G = G              # gamma
        self.L = L              # length
        self.Cd = Cd            # flow coefficient factor
        self.dth = dth
        self.Macht = False      # Mach placeholder
        self.ct = 0             # Mach ct

        if Cd and not dth:
            self.dth = d

        # handle information (0 = false, 1 closed end, 2 open end (atmos), 3 bell mouth inflow, 4 cylinder)
        self.first = 0
        self.last = 0

        self.XR = XR            # pressure superposition ratio at the right end of the mesh
        self.XL = XL            # pressure superposition ratio at the left end of the mesh
        self.XR1 = XR1
        self.XL1 = XL1

        self.cs = 0             # reserved
        self.Tw = T0            # reserved, mesh wall temperature (important for exhaust duct/pipe)
        self.Xpf = XL           # reserved
        self.Xqf = XR1          # reserved
        self.dQf = 0            # reserved
        self.dQh = 0            # reserved

        self.Vj = pi/4 * self.d**2 * self.L                             # volume of the mesh

        # averages
        self.Xj = ((self.XR+self.XL-1)+(self.XR1+self.XL1-1)) / 2
        self.pj = self.p0*self.Xj**(G7(self.G))
        self.rhoj = self.rho0*self.Xj**(G5(self.G))
        self.Tj = self.T0*self.Xj**2

        # acoustic velocities for dt
        self.aR = self.a0*(G6(self.G)*self.XR-G4(self.G)*self.XL-1)
        self.aL = self.a0*(G6(self.G)*self.XL-G4(self.G)*self.XR-1)
        self.aR1 = self.a0*(G6(self.G)*self.XR1-G4(self.G)*self.XL1-1)
        self.aL1 = self.a0*(G6(self.G)*self.XL1-G4(self.G)*self.XR1-1)

        self.dt = min([self.L / abs(self.aR), self.L / abs(self.aR1), self.L / abs(self.aL), self.L / abs(self.aL1)])

        self.a0 = a0
        self.aj = a0*self.Xj                                            # acoustic velocity
        self.mj = self.rhoj*self.Vj                                     # mass of the gas in the mesh

    def update(self):
        # acoustic velocities for dt
        self.aR = self.a0*(G6(self.G)*self.XR-G4(self.G)*self.XL-1)
        self.aL = self.a0*(G6(self.G)*self.XL-G4(self.G)*self.XR-1)
        self.aR1 = self.a0*(G6(self.G)*self.XR1-G4(self.G)*self.XL1-1)
        self.aL1 = self.a0*(G6(self.G)*self.XL1-G4(self.G)*self.XR1-1)

        self.dt = min([self.L / abs(self.aR), self.L / abs(self.aR1), self.L / abs(self.aL), self.L / abs(self.aL1)])

        self.Xj = ((self.XR + self.XL - 1) + (self.XR1 + self.XL1 - 1)) / 2
        self.pj = self.p0 * self.Xj ** (G7(self.G))
        self.rhoj = self.rho0 * self.Xj ** (G5(self.G))
        self.Tj = self.T0 * self.Xj ** 2

        self.aj = self.a0 * self.Xj  # acoustic velocity
        self.mj = self.rhoj * self.Vj  # mass of the gas in the mesh


    def friction(self, dt):
        E = self.a0*dt / self.L
        A = E*(self.XR1-self.XR)
        B = E*(self.XL-self.XL1)
        if A != 0: C = self.XR1 / A
        else: C = 0
        if B != 0: D = self.XL / B
        else: D = 0


        if self.XR == self.XR1 and self.XL1 == self.XL:
            Xq = self.XL1
            Xp = self.XR1
        elif self.XR == self.XR1:
            Xp = self.XR1
            Xq = (1 + D + G4(self.G)*Xp)/(G6(self.G) + 1/B)
        elif self.XL == self.XL1:
            Xq = self.XL1
            Xp = (1 + C + G4(self.G)*Xq)/(G6(self.G) + 1/A)
        else:
            FR = (G6(self.G) + 1 / A) / G4(self.G)
            FL = (G6(self.G) + 1 / B) / G4(self.G)
            Xp = (1 + D + FL + FL*C)/(G4(self.G)*(FR*FL-1))
            Xq = (1 + C + FR + FR*D)/(G4(self.G)*(FR*FL-1))

        # 2.18.20 / 2.18.21, page 149
        # XR1_new = Xp + {± friction effects ± area change effects}
        # XL_new = Xq + {± friction effects ± area change effects}

        # New superposition values for particles during a step:
        Xs = Xp + Xq - 1
        ps = self.p0 * Xs**G7(self.G)
        rhos = self.rho0*Xs**G5(self.G)
        Ts = self.T0 * Xs**2                        # ! Check this... 2.18.25, not a superposition value?
        cs = abs(G5(self.G)*self.a0*(Xp-Xq))        # Value only, no sign convention

        if cs == 0:                                 # No movement, no friction, no heat...
            return True
        # Changes due to friction during a computation step:
        Aj = pi*self.d*self.L                       # surface area of the mesh "wall"

        mu = 7.457 * 10**(-6) + 4.1547*10**(-8) * Ts - 7.4793 * 10**(-12) * Ts**2 # [kg / ms]
        Re = (rhos*self.d*cs)/(mu*Ts)               # Reynolds number
        Cf = 0.0791 / Re**(1/4)
        self.dQf = abs((Aj*Cf*rhos*cs**3*dt)/2)     # Work done by friction. Must be positive... no negative friction :)
        dpf = (2*Cf*rhos*cs**3*dt) / self.d         # Pressure loss caused by friction


        if Xs > 1:                                  # For a compression wave
            psf = ps-dpf                            # Pressure decreases by the amount of dpf
        else:                                       # For a expansion wave
            psf = ps+dpf                            # Pressure increases by the amount of dpf

        Xsf = (psf/self.p0)**G17(self.G)            # New Xs, taking in account the pressure loss due to friction
        csf = cs + (psf-ps)/(rhos*cs)



        self.Xpf = (1 + Xsf + csf/(G5(self.G)*self.a0))/2
        self.Xqf = 1 + Xsf-self.Xpf

        self.cs = cs

        # Changes due to heat transfer during a computation step:
        Ck = 6.1944 * 10**(-3) + 7.3814*10**(-5) * Ts - 1.2491 * 10**(-8) * Ts**2 # [W / mK]
        Ch = (Ck*Cf*Re)/(2*self.d)
        self.dQh = Aj * Ch * (self.Tw - Ts) * dt    # Heat transfer from/to the wall of the mesh

    # Area change effects
    # Xpf / Xqf should be used (friction applied with move_gas())
    def boundaries(self, nd):
        st = self

        # Inside a mesh waves propagate with friction (generally true). New values at boundaries now are:
        if st.first == 0:
            st.XL = st.Xqf
            st.XR1 = st.Xpf

        # first mesh
        elif st.first == 2:
            st.XR1 = st.Xpf
            st.XL = st.Xqf
            st.XL1 = nd.Xqf
        elif st.first == 1:
            st.XR1 = st.Xpf
            st.XL1 = nd.Xqf


        # last mesh
        if st.last == 1:
            # restricted pipe
            st.XR1 = st.Xpf
            nd.XL = nd.Xqf

            # TODO: Reflected pressures, st.XL1 and nd.XR, section 2.12
            # Fig 2.8: Xp1 == Xi1, Xq2 == Xi2

        elif nd.last == 2:
            # atmosphere
            '''
            nd.XL = 2-st.Xpf
            nd.XR1 = st.Xpf

            st.XL1 = nd.XL
            st.XL = st.Xqf
            st.XR1 = st.Xpf'''

            nd.XR1 = 2-nd.Xpf
            nd.XL1 = 1.0
            nd.XL = 2-nd.Xpf
            nd.XR = st.Xpf

            st.XL1 = nd.XL
        elif st.last == 3:
            # branch
            pass

        else:
            # We still need XL1 and XR for the next time step:
            # Different gas properties -> discontinuity -> flowSolver needed
            if st.d == nd.d:
                boundaryHandler = flowSolver(st, nd)
                x = boundaryHandler.intermesh_boundary()
                nd.XR = x[0]
                st.XL1 = x[1]
            elif st.Cd == 1.0 and st.d < nd.d:
                # Expanding pipe. Using Benson's constant pressure criteria for initial guess
                boundaryHandler = flowSolver(st, nd)
                x = boundaryHandler.pipe_expansion()
            elif st.Cd == 1.0 and st.d > nd.d:
                # Contracting pipe. Using Benson's constant pressure criteria for initial guess
                boundaryHandler = flowSolver(st, nd)
                x = boundaryHandler.pipe_contraction()
            else:
                # Restricted flow, coefficient Cd involved
                boundaryHandler = flowSolver(st, nd, st)
                x = boundaryHandler.pipe_restricted()


# TODO: FIX MACH > 1 !!!!!!!
class flowSolver():
    def __init__(self, mesh1, mesh2, throttle=False):
        self.a = mesh1
        self.b = mesh2
        if throttle != False: self.c = throttle


    def intermesh_boundary(self):
        a = self.a
        b = self.b
        x, info, ier, mesg = optimize.fsolve(eq_255_257, np.array([a.Xpf, b.Xqf]), args=(a, b), full_output=True)
        return x

    def pipe_restricted(self):
        mesh1 = self.a
        mesh2 = self.b
        port = self.c
        A1 = pi*(mesh1.d/2)**2
        A2 = pi*(mesh2.d/2)**2
        Ar = A2/A1  # area ratio

        Xi1 = mesh1.Xpf
        Xi2 = mesh2.Xqf
        a02 = mesh1.a0
        ct = mesh1.cs

        # by Benson constant pressure criteria:
        Xr1 = ((1-Ar)*Xi1+2*Xi2*Ar)/(1+Ar)
        Xr2 = (2*Xi1-Xi2*(1-Ar))/(1+Ar)
        Xt = (Xr1+Xr2)/2

        #print("{:4} {:4} {:4} {:4} {:4} {:4} {:4}".format("d1", "d2", "dt", "Cd", "Ar", "Pi1", "Pi2"))
        #print("{:4} {:4} {:4} {:4} {:4} {:.2f} {:.2f}".format(mesh1.d * 1000, mesh2.d * 1000, mesh1.dth * 1000, mesh1.Cd, Ar, Xi1 ** G7(mesh1.G), Xi2 ** G7(mesh2.G)))


        sol = optimize.root(restricted_pipe_flow, np.array([Xr1, Xr2, Xt, a02, ct]), args=(mesh1, mesh2, port), method='hybr')

        Xr1 = sol['x'][0]
        Xr2 = sol['x'][1]
        Xt = sol['x'][2]
        a02 = sol['x'][3]
        ct = sol['x'][4]

        # TODO: FIX MACH > 1 !!!!!!!

        Pr1 = Xr1**G7(mesh1.G)
        Pr2 = Xr2**G7(mesh2.G)
        T02 = a02**2 / (mesh1.G * mesh1.R)
        #print()
        #print("Pr1: {:.4f}, Pr2: {:.4f}, T02: {:.2f}".format(Pr1, Pr2, T02))




    def pipe_contraction(self):
        mesh1 = self.a
        mesh2 = self.b
        A1 = pi*(mesh1.d/2)**2
        A2 = pi*(mesh2.d/2)**2
        Ar = A2/A1  # area ratio

        Xi1 = mesh1.Xpf
        Xi2 = mesh2.Xqf

        # considering the flow as isentropic, which means T01 = T01 and a01 = a02
        # also, there should be no gain in entropy

        # by Benson constant pressure criteria:
        Xr1 = ((1-Ar)*Xi1+2*Xi2*Ar)/(1+Ar)
        Xr2 = (2*Xi1-Xi2*(1-Ar))/(1+Ar)

        sol = optimize.root(isentropic_contraction, np.array([Xr1, Xr2]), args=(mesh1, mesh2), method='hybr')

        Xr1 = sol['x'][0]
        Xr2 = sol['x'][1]

        # TODO: Xr1 and Xr2 to meshes...

        #print("{:4} {:4} {:4} {:4} {:4} {:4} {:4}".format("d1", "d2", "dt", "Cd", "Ar", "Pi1", "Pi2"))
        #print("{:4} {:4} {:4} {:4} {:4} {:.2f} {:.2f}".format(mesh1.d * 1000, mesh2.d * 1000, mesh1.d * 1000, 1.0, Ar,
        #                                                      Xi1 ** G7(mesh1.G), Xi2 ** G7(mesh2.G)))
        #print("{:6} {:6}".format("P1", "P2"))
        #print("{:.4f} {:.4f}".format(Xr1**G7(mesh1.G), Xr2**G7(mesh2.G)))


    def pipe_expansion(self):
        mesh1 = self.a
        mesh2 = self.b
        A1 = pi*(mesh1.d/2)**2
        A2 = pi*(mesh2.d/2)**2
        Ar = A2/A1  # area ratio

        Xi1 = mesh1.Xpf
        Xi2 = mesh2.Xqf

        # by Benson constant pressure criteria:
        Xr1 = ((1-Ar)*Xi1+2*Xi2*Ar)/(1+Ar)
        Xr2 = (2*Xi1-Xi2*(1-Ar))/(1+Ar)
        T02 = mesh1.T0

        sol = optimize.root(non_isentropic_expansion, np.array([Xr1, Xr2, T02]), args=(mesh1, mesh2))

        if mesh1.Macht:
            # Means mesh1 had Mach number over 0.65 -> flow separation should occur
            sol = optimize.root(separated_flow_expansion, np.array([Xr1, Xr2, T02]), args=(mesh1, mesh2))



        #print("{:4} {:4} {:4} {:4} {:4} {:4} {:4}".format("d1", "d2", "dt", "Cd", "Ar", "Pi1", "Pi2"))
        #print("{:4} {:4} {:4} {:4} {:4} {:.2f} {:.2f}".format(mesh1.d*1000, mesh2.d*1000, mesh1.d*1000, 1.0, Ar, Xi1**G7(mesh1.G), Xi2**G7(mesh2.G)))

        Xr1 = sol['x'][0]
        Xr2 = sol['x'][1]
        T02 = sol['x'][2]

        # TODO: Xr1 and Xr2, T02 to meshes...

        #print()
        #print("{:6} {:6} {:4}".format("P1", "P2", "T02"))
        #print("{:.4f} {:.4f} {:.3f}".format(Xr1**G7(mesh1.G), Xr2**G7(mesh2.G), T02))

    def cylinder_inflow(self):
        cyl = self.a
        pipe = self.b
        port = self.c

        # gas properties are those of the gas in the pipe, not the cyl (pipe gas is flowing)
        G = pipe.G
        R = pipe.R

        X1 = (cyl.p / cyl.p0)**G17(G)
        Xi2 = pipe.Xqf

        T01 = cyl.T / X1 ** 2
        rho02 = pipe.p0 / (R*pipe.T0)

        a01 = np.sqrt(G*R*T01)
        a02 = np.sqrt(G*R*pipe.T0)

        Xr2 = 1
        ct = G5(G)*a02*(Xi2-X1)

        P1 = cyl.p / cyl.p0
        T02 = pipe.T0
        dt = port.d
        d2 = pipe.d
        Cd = port.Cd
        Pi2 = pipe.Xqf ** G7(G)

        print("P1: {:.3f}".format(P1))
        print("T02: {:.3f}".format(T02-273))
        print("dt: {:.3f}".format(dt))
        print("d2: {:.3f}".format(d2))
        print("Cd: {:.3f}".format(Cd))
        print("Pi2: {:.3f}".format(Pi2))
        sol = optimize.root(cylinder_inflow, np.array([Xr2, ct]), args=(cyl, pipe, port))

        Xr2 = sol['x'][0]
        ct = sol['x'][1]
        Pr2 = Xr2**G7(G)
        Ps2 = (Xi2+Xr2-1)**G7(G)
        Ts2 = (Xi2+Xr2-1)**2*pipe.T0-273
        Pt = P1
        Tt = (a02*X1)**2 / (G*R) - 273
        cs2 = G5(G) * a02 * (Xi2 - Xr2)
        Xs2 = (Xi2+Xr2-1)

        print("Pr2: {:.3f}  Ps2: {:.3f} Ts2: {:.1f} Pt: {:.2f} Tt: {:.1f} ct: {:.1f} cs2: {:.1f} a01: {:.1f} a02: {:.1f}".format(Pr2, Ps2, Ts2, Pt, Tt, ct, cs2, a01, a02))
        print("Pr2: {:.3f}, Pi2: {:.3f}, Ps2: {:.3f}".format(Xr2 ** G7(G), Xi2 ** G7(G), Xs2 ** G7(G)))



    def cylinder_outflow(self):
        cyl = self.a
        pipe = self.b
        port = self.c

        X1 = (cyl.p / cyl.p0) ** G17(cyl.G)
        T01 = cyl.T / X1 ** 2

        Xr2 = (cyl.p/cyl.p0)**G17(cyl.G)
        Xt = Xr2
        a01 = np.sqrt(cyl.G * cyl.R * T01)
        a02 = a01
        X1 = (cyl.p / cyl.p0) ** G17(cyl.G)
        Xi2 = pipe.Xqf
        ct = G5(cyl.G)*a01*(X1-Xi2)

        sol = optimize.root(cylinder_outflow, np.array([Xr2, Xt, a02, ct]), args=(cyl, pipe, port))

        P1 = (cyl.p / cyl.p0)
        T1 = T01 * X1 ** 2 - 273
        dth = port.d
        d2 = pipe.d
        Cd = port.Cd
        Pi2 = Xi2 ** G7(cyl.G)

        x = sol['x']
        Xr2 = x[0]
        Xt = x[1]
        a02 = x[2]
        ct = abs(x[3])

        if port.Mach == True:
            sol = optimize.root(cylinder_outflow_sonic, np.array([Xr2, Xt, a02]), args=(cyl, pipe, port))
            x = sol['x']
            Xr2 = x[0]
            Xt = x[1]
            a02 = x[2]
            ct = a01*Xt

        Pr2 = Xr2 ** G7(cyl.G)
        Ps2 = (Xr2 + Xi2 - 1) ** G7(cyl.G)
        Mach = ct / (a01*Xt)
        T02 = a02**2 / (cyl.G*cyl.R)
        Ts2 = (Xi2+Xr2-1)**2*T02-273

        Pt = Xt ** G7(cyl.G)
        Tt = (a01 * Xt) ** 2 / (cyl.G * cyl.R) - 273
        cs2 = G5(cyl.G) * a02 * (Xr2 - Xi2)

        if Xr2 < 1: Xr2 = 1

        return Xr2, Xt, a02

        #print("{:4} {:4} {:4} {:4} {:4} {:4}".format("P1", "T1", "dt", "d2", "Cd", "Pi2"))
        #print("{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}".format(P1, T1, dth, d2, Cd, Pi2))
        #print()

        #print("Pr2: {:.4f} Ps2: {:.4f} Ts2: {:.1f} Pt: {:.4f} Tt: {:.1f} ct: {:.2f} cs2: {:.2f} a01: {:.1f} a02:{:.1f} Mach: {:.2f}".format(Pr2, Ps2, Ts2, Pt, Tt, ct, cs2, a01, a02, Mach))
        #print(sol)