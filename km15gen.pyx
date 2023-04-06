from gepard.fits import th_KM15
import gepard as g
import numpy as np
import argparse
from utils.const import *
from utils.physics import *

def nu(xB, Q2, t=0, phi=0):
  return Q2/(2*M*xB)

def y(xB, Q2, t=0, phi=0):
  return nu(xB, Q2, t, phi)/10.604

def W2(xB, Q2, t = 0, phi = 0):
  return M*M+2.0*M*nu(xB, Q2, t, phi)-Q2

def tmin2(xB, Q2, t=0, phi=0):
  return -0.5*((Q2/xB-Q2)*(Q2/xB-np.sqrt((Q2/xB)**2+4*M*M*Q2))+2*M*M*Q2)/W2(xB, Q2, t, phi)

def tmax2(xB, Q2, t=0, phi=0):
  return -0.5*((Q2/xB-Q2)*(Q2/xB+np.sqrt((Q2/xB)**2+4*M*M*Q2))+2*M*M*Q2)/W2(xB, Q2, t, phi)

def bhdvcs(xBmin, xBmax, Q2min, Q2max, tmin, tmax, xs_only = 0):
  elPold  = 2*np.random.randint(2) - 1
  xBd     = xBmin + (xBmax - xBmin) * np.random.rand()
  Q2d     = Q2min + (Q2max - Q2min) * np.random.rand()
  yd      = y(xBd, Q2d)
  if (yd < ymin) or (yd>ymax):
    return
  td      = tmin + (tmax - tmin) * np.random.rand()
  if (td < -tmin2(xBd, Q2d)) or (td > -tmax2(xBd, Q2d)):
    return
  phigd    = np.random.rand()*2*np.pi
  xs = printKM(xBd, Q2d, td, phigd, pol = elPold, mode = 6)
  if xs_only:
    return xs
  return elPold, xBd, Q2d, yd, td, phigd, xs

def getphoton(xBd, Q2d, td, phigd, phield):
  nud = Q2d / 2 /M / xBd
  Esc = Ed  - nud
  yb = nud/Ed
  costel = 1 - Q2d/(2*Ed*Esc)
  sintel = np.sqrt(1-costel**2)

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

  cosphe = np.cos(phield)
  sinphe = np.sin(phield)
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


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("-Ed", "--Ed", type = float, default = 10.604)
  parser.add_argument("-trig", "--trig", type = int, default = 1)
  parser.add_argument("-xBmin", "--xBmin", type = float, default = 0.05)
  parser.add_argument("-xBmax", "--xBmax", type = float, default = 0.75)
  parser.add_argument("-Q2min", "--Q2min", type = float, default = 0.9)
  parser.add_argument("-Q2max", "--Q2max", type = float, default = 11)
  parser.add_argument("-tmin", "--tmin", type = float, default = 0.085)
  parser.add_argument("-tmax", "--tmax", type = float, default = 1.79)
  parser.add_argument("-ymin", "--ymin", type = float, default = 0.19)
  parser.add_argument("-ymax", "--ymax", type = float, default = 0.85)
  parser.add_argument("-find_max", "--find_max", action = 'store_true')
  args = parser.parse_args()


  M = 0.938272081 # target mass
  Ed    = args.Ed
  xBmin = args.xBmin
  xBmax = args.xBmax
  Q2min = args.Q2min
  Q2max = args.Q2max
  tmin  = args.tmin
  tmax  = args.tmax
  ymin  = args.ymin
  ymax  = args.ymax

  if args.find_max:
    xs_array = []
    while len(xs_array) < args.trig:
      xs = bhdvcs(xBmin, xBmax, Q2min, Q2max, tmin, tmax, xs_only = 1)
      if xs:
        xs_array.append(xs)

    print(np.max(xs_array))