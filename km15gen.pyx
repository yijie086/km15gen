from gepard.fits import th_KM15
import gepard as g
import numpy as np
from utils.const import *
from utils.physics import *

cpdef double nu(double xB, double Q2):
  return Q2/(2*M*xB)

cpdef double y(double xB, double Q2):
  return nu(xB, Q2)/10.604

cpdef double W2(double xB, double Q2):
  return M*M+2.0*M*nu(xB, Q2)-Q2

cpdef double tmin2(double xB, double Q2):
  return -0.5*((Q2/xB-Q2)*(Q2/xB-np.sqrt((Q2/xB)**2+4*M*M*Q2))+2*M*M*Q2)/W2(xB, Q2)

cpdef double tmax2(double xB, double Q2):
  return -0.5*((Q2/xB-Q2)*(Q2/xB+np.sqrt((Q2/xB)**2+4*M*M*Q2))+2*M*M*Q2)/W2(xB, Q2)

cpdef double bhdvcs(double xBmin, double xBmax, 
  double Q2min, double Q2max, double tmin, double tmax,
  double ymin, double ymax, int xs_only = 0):
  
  cdef int elPold  = 2*np.random.randint(2) - 1
  cdef double xBd     = xBmin + (xBmax - xBmin) * np.random.rand()
  cdef double Q2d     = Q2min + (Q2max - Q2min) * np.random.rand()
  cdef double yd      = y(xBd, Q2d)
  if (yd < ymin) or (yd>ymax):
    return 0
  cdef double td      = tmin + (tmax - tmin) * np.random.rand()
  if (td < -tmin2(xBd, Q2d)) or (td > -tmax2(xBd, Q2d)):
    return 0
  cdef double phigd    = np.random.rand()*2*np.pi
  cdef double xs = printKM(xBd, Q2d, td, phigd, pol = elPold, mode = 6)
  if xs_only:
    return xs
  return elPold, xBd, Q2d, yd, td, phigd, xs

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
