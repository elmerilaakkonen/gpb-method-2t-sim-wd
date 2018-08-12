from math import *
import numpy as np


def G3(gamma):
    return (4-2*gamma)/(gamma-1)

def G4(gamma):
    return (3-gamma)/(gamma-1)

def G5(gamma):
    return (2)/(gamma-1)

def G6(gamma):
    return (gamma+1)/(gamma-1)

def G7(gamma):
    return (2*gamma)/(gamma-1)

def G17(gamma):
    return (gamma-1)/(2*gamma)

def G35(gamma):
    return (gamma)/(gamma-1)

def G67(gamma):
    return (gamma+1)/(2*gamma)


# equations 2.5.5 and 2.5.7 combined to a system of equations to be solved
def eq_255_257(data, a, b):
    X1d, X2d = data
    # using the nomenclature by Blair in page 87:
    X1 = a.Xpf
    X2 = b.Xqf
    return G5(a.G)*a.a0*(X1 - X2d) - G5(b.G)*b.a0*(X1d - X2), (X1+X2d-1)**G7(a.G)-(X1d+X2-1)**G7(b.G)

# equations 2.10.7, 2.10.8 and 2.10.9 combined to a system of equations to be solved
def eq_2107_2108_2109(data, a, b):
    Xr1, Xr2, T02 = data
    Xi1 = a.Xpf
    Xi2 = b.Xqf
    rho01 = a.p0 / (a.R*a.T0)
    rho01 = a.p0 / (a.R*a.T0)
    rho02 = b.p0 / (b.R*T02)

def restricted_pipe_flow(data, mesh1, mesh2, port):
    Xr1, Xr2, Xt, a02, ct = data
    A1 = pi * (mesh1.d / 2) ** 2
    A2 = pi * (mesh2.d / 2) ** 2
    Xi1 = mesh1.Xpf
    Xi2 = mesh2.Xqf

    G = mesh1.G
    R = mesh1.R

    rho01 = mesh1.p0/(R*mesh1.T0)
    rho0t = rho01
    rho02 = mesh2.p0/(R*mesh2.T0)
    a01 = np.sqrt(G*R*mesh1.T0)


    rho01 = mesh1.p0 / (G*R*mesh1.T0)
    rho02 = mesh2.p0 / (G*R*mesh2.T0)

    At = pi * (port.dth / 2)**2
    Cd = port.Cd

    # sonic
    Macht = ct / (a01*Xt)
    if Macht > 1:
        mesh1.Macht = True
        mesh1.ct = a01 * Xt
        ct = a01 * Xt
    else:
        mesh1.Macht = False


    conti1 = rho01*(Xi1+Xr1-1)**G5(G)*A1*a01*(Xi1-Xr1) + rho02*(Xi2+Xr2-1)**G5(G)*A2*a02*(Xi2-Xr2)
    conti2 = (Xi1+Xr1-1)**G5(G)*A1*a01*(Xi1-Xr1)-Xt**G5(G)*Cd*At*ct
    thd1 = ((G5(G)*a01*(Xi1-Xr1))**2 + G5(G)*a01**2*(Xi1+Xr1-1)**2) - ((G5(G)*a02*(Xi2-Xr2))**2+G5(G)*a02**2*(Xi2+Xr2-1)**2)
    thd2 = G5(G)*((a01*(Xi1+Xr1-1))**2-(a01*Xt)**2)+(G5(G)*a01*(Xi1-Xr1))**2 - ct**2
    moment = mesh1.p0 * A2 * (Xt**G7(G) - (Xi2+Xr2-1)**G7(G)) + (rho01*(Xi1+Xr1-1)**G5(G)*A1*G5(G)*a01*(Xi1-Xr1)) * (ct+G5(G)*a01*(Xi2-Xr2))

    return conti1, conti2, thd1, thd2, moment


def isentropic_contraction(data, mesh1, mesh2):
    Xr1, Xr2 = data
    A1 = pi * (mesh1.d / 2) ** 2
    A2 = pi * (mesh2.d / 2) ** 2
    Xi1 = mesh1.Xpf
    Xi2 = mesh2.Xqf

    G = mesh1.G
    R = mesh1.R

    conti = (Xi1+Xr1-1)**G5(G)*A1*(Xi1-Xr1) + (Xi2+Xr2-1)**G5(G)*A2*(Xi2-Xr2)
    thd = ((Xi1-Xr1)**2+(Xi1+Xr1-1)**2)-((Xi2-Xr2)**2+(Xi2+Xr2-1)**2)

    return conti, thd

def separated_flow_expansion(data, mesh1, mesh2):
    print("Mach > 0.65, flow separation...")
    Xr1, Xr2, T02 = data
    A1 = pi * (mesh1.d / 2) ** 2
    A2 = pi * (mesh2.d / 2) ** 2

    G = mesh1.G
    R = mesh1.R

    Xi1 = mesh1.Xpf
    Xi2 = mesh2.Xqf

    Xs1 = (Xi1+Xr1-1)
    Xs2 = (Xi2+Xr2-1)

    a01 = np.sqrt(G * R * mesh1.T0)
    a02 = np.sqrt(G * R * T02)

    rho01 = mesh1.p0 / (G * R * mesh1.T0)
    rho02 = mesh2.p0 / (G * R * T02)

    # Check for sonic conditions (not possible)
    Ms1 = (G5(G) * (Xi1 - Xr1)) / (Xi1 + Xr1 - 1)
    if Ms1 > 1.00:
        Xr1 = (1 + G4(G) * Xi1) / G6(G)

    if Ms1 > 0.65:
        mesh1.Macht = True
    else:
        mesh1.Macht = False

    conti = rho01 * (Xi1 + Xr1 - 1) ** G5(G) * A1 * G5(G) * a01 * (Xi1 - Xr1) \
            + rho02 * (Xi2 + Xr2 - 1) ** G5(G) * A2 * G5(G) * a02 * (Xi2 - Xr2)

    thd = ((G5(G) * a01 * (Xi1 - Xr1)) ** 2 + G5(G) * a01 ** 2 * (Xi1 + Xr1 - 1) ** 2) \
          - ((G5(G) * a02 * (Xi2 - Xr2)) ** 2 + G5(G) * a02 ** 2 * (Xi2 + Xr2 - 1) ** 2)

    constant_pressure = Xs1**G7(G) - Xs2**G7(G)

    return conti, thd, constant_pressure


def non_isentropic_expansion(data, mesh1, mesh2):
    Xr1, Xr2, T02 = data
    A1 = pi * (mesh1.d / 2) ** 2
    A2 = pi * (mesh2.d / 2) ** 2

    G = mesh1.G
    R = mesh1.R

    Xi1 = mesh1.Xpf
    Xi2 = mesh2.Xqf

    a01 = np.sqrt(G*R*mesh1.T0)
    a02 = np.sqrt(G*R*T02)

    rho01 = mesh1.p0 / (G*R*mesh1.T0)
    rho02 = mesh2.p0 / (G*R*T02)

    # Check for sonic conditions (not possible)
    Ms1 = (G5(G) * (Xi1 - Xr1)) / (Xi1 + Xr1 - 1)
    if Ms1 > 1.00:
        Xr1 = (1 + G4(G) * Xi1) / G6(G)

    if Ms1 > 0.65:
        mesh1.Macht = True
    else:
        mesh1.Macht = False

    conti = rho01*(Xi1+Xr1-1)**G5(G)*A1*G5(G)*a01*(Xi1-Xr1) \
        + rho02*(Xi2+Xr2-1)**G5(G)*A2*G5(G)*a02*(Xi2-Xr2)

    thd = ((G5(G)*a01*(Xi1-Xr1))**2+G5(G)*a01**2*(Xi1+Xr1-1)**2) \
        - ((G5(G)*a02*(Xi2-Xr2))**2+G5(G)*a02**2*(Xi2+Xr2-1)**2)

    moment = mesh2.p0*A2*((Xi1+Xr1-1)**G7(G)-(Xi2+Xr2-1)**G7(G)) \
        + (rho01*(Xi1+Xr1-1)**G5(G)*A1*G5(G)*a01*(Xi1-Xr1))\
        * (G5(G)*a01*(Xi1-Xr1)+G5(G)*a02*(Xi2-Xr2))

    return conti, thd, moment

def cylinder_inflow(data, cyl, pipe, port):
    Xr2, ct = data
    R = pipe.R
    G = pipe.G

    X1 = (cyl.p/cyl.p0)**G17(G)
    Xi2 = pipe.Xqf
    T02 = pipe.T0
    rho02 = cyl.p0 / (R * T02)
    a02 = np.sqrt(G*R*T02)
    At = pi * (port.d/2)**2
    A2 = pi * (pipe.d/2)**2
    Cd = port.Cd

    conti = rho02*X1**G5(G)*Cd*At*ct-rho02*(Xi2+Xr2-1)**G5(G)*A2*G5(G)*a02*(Xi2-Xr2)
    first = (G5(G)*(a02*(Xi2+Xr2-1))**2-G5(G)*(a02*X1)**2)
    second = ((G5(G)*a02*(Xi2-Xr2))**2-ct**2)
    thd = first+second

    return conti, thd

def cylinder_outflow_sonic(data, cyl, pipe, port):
    Xr2, Xt, a02 = data
    R = cyl.R
    G = cyl.G

    X1 = (cyl.p / cyl.p0) ** G17(G)
    Xi2 = pipe.Xqf
    T01 = cyl.T / X1**2

    T02 = (a02 ** 2) / (G * R)

    rho01 = cyl.p0 / (R * T01)
    rho02 = cyl.p0 / (R * T02)
    a01 = np.sqrt(G * R * T01)

    # Sonic particle velocity
    ct = a01*Xt

    Cd = port.Cd
    At = pi * (port.d/2)**2
    A2 = pi * (pipe.d/2)**2


    conti = rho01 * Xt ** G5(G) * Cd * At * ct - rho02*(Xi2+Xr2-1)**G5(G)*A2*G5(G)*a02*(Xi2-Xr2)
    thd1 = G5(G)*(a01*X1)**2 - ((G5(G)*a02*(Xi2-Xr2))**2+G5(G)*a02**2*(Xi2+Xr2-1)**2)
    thdmom = G5(G)*((a01*X1)**2-(a01*Xt)**2)-(a01*Xt)**2

    return conti, thd1, thdmom


def cylinder_outflow(data, cyl, pipe, port):
    Xr2, Xt, a02, ct = data
    R = cyl.R
    G = cyl.G

    X1 = (cyl.p / cyl.p0) ** G17(G)
    Xi2 = pipe.Xqf
    T01 = cyl.T / X1**2
    T02 = (a02 ** 2) / (G * R)


    rho01 = cyl.p0 / (R * T01)
    rho02 = cyl.p0 / (R * T02)
    a01 = np.sqrt(G * R * T01)
    a1 = a01 * X1
    rhot = rho01 * Xt ** G5(G)
    Tt = (a01 * Xt)**2 / (G*R)

    Cd = port.Cd
    At = port.A
    A2 = pi * (pipe.d/2)**2


    conti = rho01 * Xt ** G5(G) * Cd * At * ct - rho02*(Xi2+Xr2-1)**G5(G)*A2*G5(G)*a02*(Xi2-Xr2)
    thd1 = G5(G)*(a01*X1)**2 - ((G5(G)*a02*(Xi2-Xr2))**2+G5(G)*a02**2*(Xi2+Xr2-1)**2)
    thd2 = G5(G)*((a01*X1)**2-(a01*Xt)**2)-ct**2
    moment = cyl.p0 * (Xt**G7(G)-(Xi2+Xr2-1)**G7(G))+(rho02*(Xi2+Xr2-1)**G5(G)*G5(G)*a02*(Xi2-Xr2))*(ct-G5(G)*a02*(Xi2-Xr2))

    Mt = ct/(a01*Xt)
    if Mt > 1.0:
        port.Mach = True
    else:
        port.Mach = False

    return conti, thd1, thd2, moment