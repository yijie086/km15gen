# distutils: language = c++
from gepard.fits import th_KM15
import gepard as g
import numpy as np
from utils.const import *
from utils.physics import *
from scipy.special import spence
from libcpp.vector cimport vector

cpdef double nu(double xB, double Q2):
  return Q2/(2*M*xB)

cpdef double y(double xB, double Q2):
  return nu(xB, Q2)/10.604

cpdef double W2(double xB, double Q2):
  return M*M+2.0*M*nu(xB, Q2)-Q2

cpdef double tmin2(double xB, double Q2):
  return -0.5*((Q2/xB-Q2)*(Q2/xB-np.sqrt((Q2/xB)**2+4.*M*M*Q2))+2.*M*M*Q2)/W2(xB, Q2)

cpdef double tmax2(double xB, double Q2):
  return -0.5*((Q2/xB-Q2)*(Q2/xB+np.sqrt((Q2/xB)**2+4.*M*M*Q2))+2.*M*M*Q2)/W2(xB, Q2)

cpdef double genEvent(double xBmin, double xBmax, 
  double Q2min, double Q2max, double tmin, double tmax,
  double ymin, double ymax, int rad = 0, double Ed = 10.604):

  cdef int elPold  = 2*np.random.randint(2) - 1
  cdef double xBd     = xBmin + (xBmax - xBmin) * np.random.rand()
  cdef double Q2d     = Q2min + (Q2max - Q2min) * np.random.rand()
  cdef double yd      = y(xBd, Q2d)

  if (yd < ymin) or (yd > ymax):
    vec.push_back(0)
    return 0

  cdef double td      = tmin + (tmax - tmin) * np.random.rand()
  
  if (td < -tmin2(xBd, Q2d)) or (td > -tmax2(xBd, Q2d)):
    vec.push_back(0)
    return 0

  cdef double phigd     = np.random.rand()*2.*np.pi
  cdef double phield    = np.random.rand()*2.*np.pi

  cdef vector[double] xsvec = bhdvcs(xBd, Q2d, td, phigd, rad, Ed)


cpdef vector[double] bhdvcs(double xBd, double Q2d, double td, double phigd, int rad = 0, double Ed = 10.604):

  cdef vector[double] vec
  cdef double costel, afac, dE1, nud, Esc, dE2, Eprime_e, E_e, eta, deltaEs, delta_vertex, delta_vac, delta_R, delta_vvr, deld  
  cdef double rho, aks, delta_1, Eprime_p, pprime_p, delta_2, xs_born, xs       

  if rad:
    # Follow procedure of M. Vanderhaeghen et al., PHYSICAL REVIEW C 62 025501
    afac = alpha/np.pi * (np.log(Q2d/me**2) - 1.)
    dE1 = np.random.rand()**(1/afac) * Ed
    E_e       = Ed  - dE1                            #A71

    nud = Q2d / 2. /M / xBd     
    Esc = Ed  - nud                
    costel = 1 - Q2d/(2*E_e*Esc)

    Q2d = Q2d * E_e/Ed
    nud     = nud - dE1
    xBd     = Q2d/ 2. / M / nud

    xs_born = printKM(xBd, Q2d, td, phigd, pol = elPold, mode = 5, E = Ed - dE1)

    afac = alpha/np.pi * (np.log(Q2d/me**2) - 1.)
    dE2 = np.random.rand()**(1/afac) * Esc
    Eprime_e  = Esc - dE2

    # costel = 1-Q2d/2/Eprime_e/E_e

    # eta       = E_e/Esc               #82
    # deltaEs   = eta * dE2             #82
    # delta_vvr = delta_vac + delta_vertex + delta_R #A71
    delta_vertex = alpha/np.pi * ( 1.5*np.log(Q2d/me**2) - 2 - 1./2. * np.log(Q2d/me**2)**2  + np.pi**2/6.)
    delta_vac    = alpha/np.pi * 2./3. * ( -5./3. + 1. * np.log(Q2d/me**2))
    delta_R      = alpha/np.pi * ( - 0.5 * np.log(E_e/Eprime_e)**2 
                                  + 0.5 * np.log(Q2d/me**2)**2- np.pi**2./3. + spence(1-((1+costel)/2.)) )

    delta_vvr = alpha/np.pi * ( (3./2.+2/3.)*np.log(Q2d/me**2)
                 -28./9. - 0.5 * np.log(E_e/Eprime_e)**2 - np.pi**2/6. + spence(1-((1+costel)/2.)) ) #cos^2(theta/2) = (1+cos theta)/2

    # rad dominated by electron side
    # deld      = np.sqrt(td)                 #A75
    # rho       = np.sqrt(td + 4*M**2)        #A75
    # aks       = (deld + rho)**2/4/M**2         #A75
    # delta_1   = 2*alpha/np.pi * (np.log(4*dE1*dE2 /td/aks) * np.log(eta) + spence(1-(1-eta/aks)) - spence(1-(1-1/eta/aks)) ) #A74
    # Eprime_p  = M + td/2/M
    # pprime_p  = np.sqrt(Eprime_p**2 - M**2)
    # delta_2   = alpha/np.pi * ( np.log(4*dE1*dE2 /M**2)*(Eprime_p/pprime_p*np.log(aks)-1) + 1 
    #                         + Eprime_p/pprime_p * (-0.5*np.log(aks)**2 - np.log(aks)*np.log(rho**2/M**2) + np.log(aks)
    #                         - spence( 1 - (1-1/aks/aks)) +2*spence( 1 - (-1/aks)) + np.pi**2/6 ) ) #A76
    print(dE1, dE2, delta_vertex, delta_vac, delta_R, np.exp(delta_vertex + delta_R)/(1-delta_vac * 0.5)**2, delta_vvr)
    xs      = xs_born* np.exp(delta_vertex + delta_R)/(1-delta_vac * 0.5)**2

    vec.push_back(elPold)
    vec.push_back(xBd)
    vec.push_back(Q2d)
    vec.push_back(yd)
    vec.push_back(td)
    vec.push_back(phigd)
    vec.push_back(phield)
    vec.push_back(xs)
    vec.push_back(xs_born)
    vec.push_back(dE1)
    vec.push_back(dE2)
    return vec

  else:
    xs  = printKM(xBd, Q2d, td, phigd, pol = elPold, mode = 5, E = Ed)
    vec.push_back(elPold)
    vec.push_back(xBd)
    vec.push_back(Q2d)
    vec.push_back(yd)
    vec.push_back(td)
    vec.push_back(phigd)
    vec.push_back(phield)
    vec.push_back(xs)
    return vec

def getphoton(xBd, Q2d, td, phigd, phield):
  nud = Q2d / 2 /M / xBd
  Esc = Ed  - nud
  yb = nud/Ed
  costel = 1 - Q2d/(2*Ed*Esc)
  sintel = np.sqrt(1-costel**2)
  cosphe = np.cos(phield)
  sinphe = np.sin(phield)

  V3k1 = 0
  V3k2 = 0
  V3k3 = Ed

  V3l1 = Esc*sintel*cosphe
  V3l2 = Esc*sintel*sinphe
  V3l3 = Esc*costel

  # non rad case!
  V3q1 = V3k1 - V3l1
  V3q2 = V3k2 - V3l2
  V3q3 = V3k3 - V3l3

  Ep = M + td/2/M
  Egam = nud -  td/2/M

  qmod = np.sqrt(nud**2 + Q2d)
  costVq = (Ed - Esc * costel)/qmod
  sintVq = np.sqrt(1 - costVq**2)

  costgg = (2*Egam*(M + nud) + Q2d - 2*M*nud)/(2*Egam*qmod)
  sintgg = np.sqrt(1 - costgg**2)
  Vgx = Egam * sintgg * np.cos(phigd)
  Vgy = Egam * sintgg * np.sin(phigd)
  Vgz = Egam*costgg

  V3gam1 = Vgx*costVq*cosphe - Vgz*sintVq*cosphe - Vgy*sinphe
  V3gam2 = Vgx*costVq*sinphe - Vgz*sintVq*sinphe + Vgy*cosphe
  V3gam3 = Vgx*sintVq        + Vgz*costVq

  V3p1 = V3q1 - V3gam1
  V3p2 = V3q2 - V3gam2
  V3p3 = V3q3 - V3gam3
