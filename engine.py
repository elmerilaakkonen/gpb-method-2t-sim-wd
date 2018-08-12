from math import *

class Engine:

    def __init__(self, bore, stroke, conrod, trapped_stroke, compression_ratio, n, RPM, intake, exhaust, cylinder, crankcase, atmos, debug = False, piston_intake_h = 0):
        self.d_bo = bore/1000                                   # Bore of cylinder in millimeters
        self.L_st = stroke/1000                                 # Stroke length in millimeters
        self.L_ts = trapped_stroke/1000                         # Trapped stroke in millimeters
        self.L_cr = conrod/1000                                 # Connecting rod length in millimeters
        self.n = n                                              # number of cylinders

        self.V_sv = self.n*(pi/4)*self.d_bo**2*self.L_st        # Swept volume of engine
        self.V_ts = self.n*(pi/4)*self.d_bo**2*self.L_ts        # Trapped swept volume

        self.V_cc = 512                                         # Crankcase clearance volume
        self.V_cv = 9                                           # Clearance volume (piston at tdc)

        self.CR_cc = (self.V_cc + self.V_sv)/self.V_cc          # Crankcase compression ratio
        self.CR_g = (self.V_sv+self.V_cv)/self.V_cv             # Geometric compression ratio of the engine
        self.CR_t = (self.V_ts+self.V_cv)/self.V_cv             # Trapped compression ratio of the engine

        self.RPM = RPM                                          # Revolutions per minute
        self.RPS = RPM / 60                                     # Revolutions per second

        self.exhaust = exhaust
        self.intake = intake
        self.cylinder = cylinder
        self.crankcase = crankcase
        self.atmos = atmos
        self.debug = debug
        self.angle = 0                                          # Current angle of crankshaft in degrees, TDC = 0

        self.piston_intake_h = piston_intake_h/1000

        self.single = False

        self.cylinder.eng = self
        for portX in self.cylinder.ports:
            self.cylinder.ports[portX].eng = self


    # Fig 1.10
    # Theta is crankshaft angle after tdc
    def piston_position(self, theta):
        """Returns the piston position in mm below TDC for given crankshaft angle theta in degrees ATDC"""
        theta = radians(theta)
        L_ct = self.L_st / 2
        alpha = pi/2-theta
        G = L_ct * sin(alpha)
        E = L_ct * cos(alpha)
        beta = asin(E/self.L_cr)
        F = cos(beta)*self.L_cr

        return (self.L_cr+L_ct - G - F) * 1000

    def crank_angle(self, piston_pos):
        """Returns crankshaft angle in degrees ATDC (0-90) for given piston position in mm below TDC"""
        l = self.L_cr
        r = self.L_st / 2
        A = -1 / 2 * (l**2 - r**2 - piston_pos**2) / (r * piston_pos)
        return degrees(acos(A))


    def time_to_angle(self, time):
        """Returns crankshaft angle the engine spins in the specified time (s)"""
        return time*self.RPS*360

    def step(self):
        """Spins the engine one step forward (step == dt, from the fastest propagation in a mesh)"""

        cyl = self.cylinder
        exhaust = self.exhaust

        # TODO: dt from all ducts
        dt = exhaust.dt
        dangle = self.time_to_angle(dt)
        #if self.debug: print("One step turns crankshaft {:.2f} degrees.".format(dangle))

        if self.angle + dangle < 360:
            self.angle += dangle
        else:
            self.angle = dangle - (360 - self.angle)


        for single_port in cyl.ports:
            port = cyl.ports[single_port]
            port.update()
            if single_port == 'exhaust':
                if port.A > 0 and self.single == False:
                    # port open
                    # print("open")
                    if single_port == 'exhaust':
                        cyl.boundaries(mesh=exhaust.meshes[0])
                else:
                    # print("closed")
                    cyl.boundaries(mesh=exhaust.meshes[0], closed=True)
                    if self.angle > 180:
                        self.single = True

            if self.intake:
                # TODO: CALCULATE PISTON PORTED INTAKE OPENING / CLOSING
                pass

        #print(self.angle)

        exhaust.step()