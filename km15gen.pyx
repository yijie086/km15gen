# distutils: language = c++
import numpy as np
from scipy.special import spence
from libcpp.vector cimport vector
import gepard as g
from gepard.fits import th_KM15

cdef double M = 0.938272081 # target mass
cdef double me = 0.5109989461 * 0.001 # electron mass
cdef double ycolcut = 0.000001 #0.0005 # P1 cut
cdef double alpha = 1/137.036



cpdef double printKM(double xB, double Q2, double t, double phi, int pol = 0, double E = 10.604, str model = 'km15'):
    phi = np.pi - phi
    pt1 = g.DataPoint(xB=xB, t=-t, Q2=Q2, phi=phi,
            process='ep2epgamma', exptype='fixed target', frame ='trento',
            in1energy= E, in1charge=-1, in1polarization=pol)
    pt1.prepare()
    if model == 'km15':
      return th_KM15.XS(pt1)
    elif model == 'km15_bh':
      return th_KM15.PreFacSigma(pt1)*th_KM15.TBH2unp(pt1)


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

# cpdef double xi(double xB, double Q2, double t):
#     return xB*(1-t/2/Q2)/(2-xB+xB*(-t/Q2))

cpdef double del2(double t):
    return -t

cpdef double del2q2(double t, double Q2):
    return del2(t)/Q2

cpdef double eps(double xB, double Q2):
    return 2*xB*M/np.sqrt(Q2)

cpdef double eps2(double xB, double Q2):
    return eps(xB, Q2)**2

cpdef double qeps2(double xB, double Q2):
    return 1+eps2(xB, Q2)

cpdef double sqeps2(double xB, double Q2):
    return np.sqrt(qeps2(xB, Q2))

cpdef double y1eps(double xB, double Q2):
    cdef double yd = y(xB, Q2)
    return 1 - yd - yd*yd*eps2(xB, Q2)/4

cpdef double Kfac2(double xB, double Q2, double t):
    cdef double tmind = -tmin2(xB, Q2)
    cdef double eps2d  = eps2(xB, Q2)
    return (-del2q2(t, Q2))*(1 - xB)*y1eps(xB, Q2)*(1 - (tmind)/t)*(np.sqrt(1 + eps2d) + 
            ((4*xB*(1 - xB) + eps2d)/(4*(1 - xB)))*(-(t - (tmind))/Q2))

cpdef double Kfac(double xB, double Q2, double t):
    return np.sqrt(Kfac2(xB, Q2, t))

cpdef double Jfac(double xB, double Q2, double t):
    cdef double yd = y(xB, Q2)
    cdef double eps2d  = eps2(xB, Q2)
    cdef double del2q2d = del2q2(t, Q2)
    return (1 - yd - yd*eps2d/2)*(1 + del2q2d) - (1 - xB)*(2 - yd)*del2q2d

cpdef double P1(double xB, double Q2, double t, double phi):    
    cdef double yd = y(xB, Q2)
    cdef double eps2d  = eps2(xB, Q2)
    return -(Jfac(xB, Q2, t) + 2*Kfac(xB, Q2, t)*np.cos(np.pi-np.radians(phi)))/(yd*(1 + eps2d))


cpdef double getScale(double xBmin, double xBmax, 
  double Q2min, double Q2max, double tmin, double tmax,
  double ymin, double ymax, double w2min, int rad = 0, double Ed = 10.604):

  cdef int nx=40
  cdef int nq=20
  cdef int nt=40

  cdef double dx    = (xBmax-xBmin)/nx
  cdef double dq    = (Q2max-Q2min)/nq
  cdef double dt

  cdef double xBd 
  cdef double Q2d 
  cdef double yd  
  cdef double w2d
  cdef double td

  cdef int elPold
  cdef double phigd
  cdef double dstot   = 0
  cdef double xs
  for elPold in [-1, 1]:
    for ix in range(1, nx+1):
      xBd = xBmin + ix * dx
      for iq in range(1, nq+1):
        Q2d = Q2min + iq * dq
        yd   = y(xBd, Q2d)
        w2d  = W2(xBd, Q2d)
        if (yd < ymin) or (yd > ymax):
          continue
        if (w2d < w2min):
          continue
        for it in range(1, nt+1):
          tmax = min(-tmax2(xBd, Q2d), tmax)
          tmin = max(-tmin2(xBd, Q2d), tmin)
          if (tmax<tmin):
            continue
          dt = (tmax-tmin)/nt
          td = tmin + it*dt
          if (abs(P1(xBd, Q2d, td, 0)) < ycolcut):
            continue
          xs = printKM(xBd, Q2d, td, 0, pol = elPold, E = Ed)
          if xs > dstot:
            dstot = xs
  return dstot


cpdef str genOneEvent(double xBmin, double xBmax, 
  double Q2min, double Q2max, double tmin, double tmax,
  double ymin, double ymax, double w2min, double xsmax, int rad = 0, double Ed = 10.604, str filename = "km15gen", str model = "km15"):

  cdef vector[double] kine
  cdef double cl_be, costel, afac, dE1, nud, Esc, dE2, Eprime_el_e, E_el_e, eta, deltaEs, delta_vertex, delta_vac, delta_R, delta_vvr, deld  
  cdef double rho, aks, delta_1, Eprime_p, pprime_p, delta_2, xs_born, xs_part, xs, xBd_tr, Q2d_tr, nud_tr
  cdef double V3l1, V3l2, V3l3, V3gam1, V3gam2, V3gam3, V3p1, V3p2, V3p3
  cdef double sintel, cosphe, sinphe
  cdef double vx, vy, vz
  cdef double radQ2, radxB, radEd
  cdef str result

  result = ""

  vx = 0.025*(np.random.rand() - 0.5)
  vy = 0.025*(np.random.rand() - 0.5)
  vz = np.random.rand()*5. # for now there's no offset for simulating the external radiator

  cl_be = Ed

  cdef int elPold  = 2*np.random.randint(2) - 1

  cdef double xBd     = xBmin + (xBmax - xBmin) * np.random.rand()
  cdef double Q2d     = Q2min + (Q2max - Q2min) * np.random.rand()
  cdef double yd      = y(xBd, Q2d)
  cdef double w2d     = W2(xBd, Q2d)
  cdef double td      = tmin + (tmax - tmin) * np.random.rand()
  cdef double phigd     = np.random.rand()*2.*np.pi
  cdef double phield    = np.random.rand()*2.*np.pi
  radQ2d = Q2d
  radxBd = xBd

  if (yd < ymin) or (yd > ymax):
    return result
  if (w2d < w2min):
    return result

  nud = Q2d / 2. /M / xBd     
  Esc = Ed  - nud
  costel = 1 - Q2d/(2*Ed*Esc)          
  if (td < -tmin2(xBd, Q2d)) or (td > -tmax2(xBd, Q2d)):
    return result
  if (-P1(xBd, Q2d, td, phigd) < ycolcut):
    return result
  xs_born = printKM(xBd, Q2d, td, phigd, pol = elPold, E = Ed, model = model)
  if np.isnan(xs_born):
    return result
  if np.isinf(xs_born):
    return result
  if xs_born == 0 :
    return result

  if rad:
    # external radiator
    Ed = Ed - Ed * np.random.rand()**(3./4./ (vz/929. + 0.003/8.897))
    # LH2 x0 = 929 cm, https://github.com/JeffersonLab/JPsiGen/blob/eb40dd934bb9f022873414a57e0dad9d1ccbcbdf/include/KinFunctions.h
    # Al  x0 = 8.897 cm https://pdg.lbl.gov/2022/AtomicNuclearProperties/HTML/aluminum_Al.html
    # Follow procedure of M. Vanderhaeghen et al., PHYSICAL REVIEW C 62 025501
    nud = nud - (cl_be - Ed)
    if nud <= 0:
      return result
    Q2d = Q2d * Ed/cl_be
    xBd = Q2d / 2. /M /nud

    if (td < -tmin2(xBd, Q2d)) or (td > -tmax2(xBd, Q2d)):
      return result
    if (-P1(xBd, Q2d, td, phigd) < ycolcut):
      return result
    afac = alpha/np.pi * (np.log(Q2d/me**2) - 1.)
    dE1 = np.random.rand()**(1/afac) * Ed
    E_el_e       = Ed  - dE1                            #A71 

    if dE1 >= nud:
      return result

    nud_tr  = nud - dE1

    afac = alpha/np.pi * (np.log(Q2d/me**2) - 1.)
    dE2 = np.random.rand()**(1/afac) * Esc
    Eprime_el_e = Esc + dE2
    if dE1 + dE2 >= nud:
      return result
    nud_tr  = nud_tr - dE2
    Q2d_tr  = Q2d * Eprime_el_e/Esc * E_el_e/Ed
    xBd_tr  = Q2d_tr/ 2. / M / nud_tr

    if (xBd_tr > 1) or (xBd_tr < 0):
      return result

    if (td < -tmin2(xBd_tr, Q2d_tr)) or (td > -tmax2(xBd_tr, Q2d_tr)):
      return result
    if (-P1(xBd_tr, Q2d_tr, td, phigd) < ycolcut):
      return result
    xs_part = printKM(xBd_tr, Q2d_tr, td, phigd, pol = elPold, E = Ed - dE1, model = model)
    if np.isnan(xs_part):
      return result
    if np.isinf(xs_part):
      return result
    if xs_part == 0 :
      return result

    # costel = 1-Q2d/2/Eprime_e/E_el_e

    # eta       = E_el_e/Esc               #82
    # deltaEs   = eta * dE2             #82
    # delta_vvr = delta_vac + delta_vertex + delta_R #A71
    delta_vertex = alpha/np.pi * ( 1.5*np.log(Q2d/me**2) - 2 - 1./2. * np.log(Q2d/me**2)**2  + np.pi**2/6.)
    delta_vac    = alpha/np.pi * 2./3. * ( -5./3. + 1. * np.log(Q2d/me**2))
    delta_R      = alpha/np.pi * ( - 0.5 * np.log(Ed/Esc)**2 
                                  + 0.5 * np.log(Q2d/me**2)**2- np.pi**2./3. + spence(1-((1+costel)/2.)) ) # the first term of delta_R used for the soft photon emission already.

    delta_vvr = alpha/np.pi * ( (3./2.+2/3.)*np.log(Q2d/me**2)
                 -28./9. - 0.5 * np.log(Ed/Esc)**2 - np.pi**2/6. + spence(1-((1+costel)/2.)) ) #cos^2(theta/2) = (1+cos theta)/2

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
    # print(dE1, dE2, delta_vertex, delta_vac, delta_R, np.exp(delta_vertex + delta_R)/(1-delta_vac * 0.5)**2, delta_vvr)
    
    xs      = xs_part* np.exp(delta_vertex + delta_R)/(1-delta_vac * 0.5)**2
    if np.isnan(xs):
      return result
    if np.isinf(xs):
      return result
    if xs == 0 :
      return result

  else:
    if (td < -tmin2(xBd, Q2d)) or (td > -tmax2(xBd, Q2d)):
      return result
    if (abs(P1(xBd, Q2d, td, phigd)) < ycolcut):
      return result
    xBd_tr   = xBd
    Q2d_tr   = Q2d
    xs      = xs_born

  # print(xsvec[0], xsmax, xsmax * np.random.rand())
  # if xs > xsmax * np.random.rand(): # this was for the rejection sampling but it's too slow for now.
  vz = vz - 5.5
  if rad:
    kine    = getphoton(xBd_tr, Q2d_tr, td, phigd, phield, Ed = Ed - dE1)
    V3l1, V3l2, V3l3, V3gam1, V3gam2, V3gam3, V3p1, V3p2, V3p3, costgg = kine
    if costgg>1:
      return result
    sintel = np.sqrt(1-costel**2)
    cosphe = np.cos(phield)
    sinphe = np.sin(phield)
    V3l1 = V3l1 - dE2*sintel*cosphe
    V3l2 = V3l2 - dE2*sintel*sinphe
    V3l3 = V3l3 - dE2*costel
    El     = np.sqrt(V3l1**2 + V3l2**2 + V3l3**2)
    Ep     = np.sqrt(V3p1**2 + V3p2**2 + V3p3**2 + M**2)
    Egam   = np.sqrt(V3gam1**2 + V3gam2**2 + V3gam3**2)

    with open("{}.dat".format(filename), "a") as file_out:
      if (dE1>=0.1 ) and (dE2>=0.1): #both s and p
        result = result + "5   1       1  0.0{:>4}   11   {:.3f}   1       1      {:6f}\n".format(elPold,cl_be,xs)
        result = result + "1  {: .4f}  1   11   0    4   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(radxBd, V3l1, V3l2, V3l3, radQ2d, td, vx, vy, vz)
        result = result + "2  {: .4f}  1 2212   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(phigd, V3p1, V3p2, V3p3, xBd_tr, Q2d_tr, vx, vy, vz)
        result = result + "3  {: .4f}      1   22   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(Ed, V3gam1, V3gam2, V3gam3, Egam, xs_born, vx, vy, vz)
        result = result + "4   1.      1   22   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(0, 0, dE1, dE1, 0, vx, vy, vz)
        result = result + "5   1.      1   22   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(dE2*sintel*cosphe, dE2*sintel*sinphe, dE2*costel, dE2, 0, vx, vy, vz)
      elif (dE1>=0.1) and (dE2<0.1): # s peak only
        result = result + "4   1       1  0.0{:>4}   11   {:.3f}   1       1      {:6f}\n".format(elPold,cl_be,xs)
        result = result + "1  {: .4f}  1   11   0    2   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(radxBd, V3l1, V3l2, V3l3, radQ2d, td, vx, vy, vz)
        result = result + "2  {: .4f}  1 2212   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(phigd, V3p1, V3p2, V3p3, xBd_tr, Q2d_tr, vx, vy, vz)
        result = result + "3  {: .4f}      1   22   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(Ed, V3gam1, V3gam2, V3gam3, Egam, xs_born, vx, vy, vz)
        result = result + "4   1.      1   22   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(0, 0, dE1, dE1, 0, vx, vy, vz)
      elif (dE1<0.1) and (dE2>=0.1): # p peak only
        result = result + "4   1       1  0.0{:>4}   11   {:.3f}   1       1      {:6f}\n".format(elPold,cl_be,xs)
        result = result + "1  {: .4f}  1   11   0    3   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(radxBd, V3l1, V3l2, V3l3, radQ2d, td, vx, vy, vz)
        result = result + "2  {: .4f}  1 2212   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(phigd, V3p1, V3p2, V3p3, xBd_tr, Q2d_tr, vx, vy, vz)
        result = result + "3  {: .4f}      1   22   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(Ed, V3gam1, V3gam2, V3gam3, Egam, xs_born, vx, vy, vz)
        result = result + "4   1.      1   22   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(dE2*sintel*cosphe, dE2*sintel*sinphe, dE2*costel, dE2, 0, vx, vy, vz)
      elif (dE1 < 0.1 ) and (dE2 < 0.1): # nonrad
        result = result + "3   1       1  0.0{:>4}   11   {:.3f}   1       1      {:6f}\n".format(elPold,cl_be,xs)
        result = result + "1  {: .4f}  1   11   0    1   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(radxBd, V3l1, V3l2, V3l3, radQ2d, td, vx, vy, vz)
        result = result + "2  {: .4f}  1 2212   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(phigd, V3p1, V3p2, V3p3, xBd_tr, Q2d_tr, vx, vy, vz)
        result = result + "3  {: .4f}      1   22   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  {: .4f}\n".format(Ed, V3gam1, V3gam2, V3gam3, Egam, xs_born, vx, vy, vz)

  else:
    kine    = getphoton(xBd, Q2d, td, phigd, phield)
    V3l1, V3l2, V3l3, V3gam1, V3gam2, V3gam3, V3p1, V3p2, V3p3, costgg = kine
    sintel = np.sqrt(1-costel**2)
    cosphe = np.cos(phield)
    sinphe = np.sin(phield)
    El     = np.sqrt(V3l1**2 + V3l2**2 + V3l3**2)
    Ep     = np.sqrt(V3p1**2 + V3p2**2 + V3p3**2 + M**2)
    Egam   = np.sqrt(V3gam1**2 + V3gam2**2 + V3gam3**2)
    xs     = xs_born
    with open("{}.dat".format(filename), "a") as file_out:
      result = result + "3   1       1  0.0{:>4}   11   {:.3f}   1       1      {:6f}\n".format(elPold,cl_be,xs)
      result = result + "1  {: .4f}  1     11   0    1   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  \n".format(xBd, V3l1, V3l2, V3l3, Q2d, td, vx, vy, vz)
      result = result + "2  {: .4f}  1   2212   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  \n".format(phigd, V3p1, V3p2, V3p3, xBd_tr, Q2d_tr, vx, vy, vz)
      result = result + "3  {: .4f}      1   22   0    0   {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f} {: .4f}  \n".format(Ed, V3gam1, V3gam2, V3gam3, Egam, xs_born, vx, vy, vz)
  return result

cpdef vector[double] bhdvcs(double xBd, double Q2d, double td, double phigd, int elPold, int rad = 0, double Ed = 10.604):

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

    xs_born = printKM(xBd, Q2d, td, phigd, pol = elPold, E = Ed - dE1)

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
    # print(dE1, dE2, delta_vertex, delta_vac, delta_R, np.exp(delta_vertex + delta_R)/(1-delta_vac * 0.5)**2, delta_vvr)
    
    xs      = xs_born* np.exp(delta_vertex + delta_R)/(1-delta_vac * 0.5)**2

    vec.push_back(xs)
    vec.push_back(xs_born)
    vec.push_back(xBd)
    vec.push_back(Q2d)
    vec.push_back(dE1)
    vec.push_back(dE2)
    return vec

  else:
    xs  = printKM(xBd, Q2d, td, phigd, pol = elPold, E = Ed)
    vec.push_back(xs)
    return vec

cpdef vector[double] getphoton(double xBd, double Q2d, double td, double phigd, double phield, double Ed = 10.604):
  cdef double nud = Q2d / 2 /M / xBd
  cdef double Esc = Ed  - nud
  cdef double yb = nud/Ed
  cdef double costel = 1 - Q2d/(2*Ed*Esc)
  cdef double sintel = np.sqrt(1-costel**2)
  cdef double cosphe = np.cos(phield)
  cdef double sinphe = np.sin(phield)

  cdef double V3k1 = 0
  cdef double V3k2 = 0
  cdef double V3k3 = Ed

  cdef double V3l1 = Esc*sintel*cosphe
  cdef double V3l2 = Esc*sintel*sinphe
  cdef double V3l3 = Esc*costel

  cdef double V3q1 = V3k1 - V3l1
  cdef double V3q2 = V3k2 - V3l2
  cdef double V3q3 = V3k3 - V3l3

  cdef double Ep = M + td/2/M
  cdef double Egam = nud -  td/2/M

  cdef double qmod = np.sqrt(V3q1**2 + V3q2**2 + V3q3**2)
  cdef double costVq = (Ed - Esc * costel)/qmod
  cdef double sintVq = np.sqrt(1 - costVq**2)

  cdef double costgg = (2*Egam*(M + nud) + Q2d - 2*M*nud)/(2*Egam*qmod)
  cdef double sintgg = np.sqrt(1 - costgg**2)
  cdef double Vgx = Egam * sintgg * np.cos(phigd)
  cdef double Vgy = Egam * sintgg * np.sin(phigd)
  cdef double Vgz = Egam*costgg

  cdef double V3gam1 = Vgx*costVq*cosphe - Vgz*sintVq*cosphe - Vgy*sinphe
  cdef double V3gam2 = Vgx*costVq*sinphe - Vgz*sintVq*sinphe + Vgy*cosphe
  cdef double V3gam3 = Vgx*sintVq        + Vgz*costVq

  cdef double V3p1 = V3q1 - V3gam1
  cdef double V3p2 = V3q2 - V3gam2
  cdef double V3p3 = V3q3 - V3gam3

  cdef vector[double] vec

  vec.push_back(V3l1) 
  vec.push_back(V3l2) 
  vec.push_back(V3l3) 
  vec.push_back(V3gam1) 
  vec.push_back(V3gam2) 
  vec.push_back(V3gam3) 
  vec.push_back(V3p1) 
  vec.push_back(V3p2) 
  vec.push_back(V3p3)
  vec.push_back(costgg)

  return vec
