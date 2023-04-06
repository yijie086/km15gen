#!/usr/bin/env python3
"""
Modules help pandas algebra without using ROOT.
"""
import numpy as np
import pandas as pd
from utils.const import *
from gepard.fits import th_KM15
import gepard as g
import os, subprocess
from functools import reduce
# from scipy import stats
from scipy.stats import chi2

def dot(vec1, vec2):
	# dot product of two 3d vectors
    return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]

def mag(vec1):
	# L2 norm of vector
	return np.sqrt(dot(vec1, vec1))

def mag2(vec1):
	# square of L2 norm
	return	dot(vec1, vec1)

def cosTheta(vec1, vec2):
	# cosine angle between two 3d vectors
    return dot(vec1,vec2)/np.sqrt(mag2(vec1) * mag2(vec2))

def angle(vec1, vec2):
	# angle between two 3d vectors
	return 180/np.pi*np.arccos(np.minimum(1, cosTheta(vec1, vec2)))

def cross(vec1, vec2):
	# cross product of two 3d vectors
    return [vec1[1]*vec2[2]-vec1[2]*vec2[1], vec1[2]*vec2[0]-vec1[0]*vec2[2], vec1[0]*vec2[1]-vec1[1]*vec2[0]]

def vecAdd(gam1, gam2):
	# add two 3d vectors
	return [gam1[0]+gam2[0], gam1[1]+gam2[1], gam1[2]+gam2[2]]

def pi0Energy(gam1, gam2):
	# reconstructed pi0 energy of two 3d photon momenta
	return mag(gam1)+mag(gam2)

def pi0InvMass(gam1, gam2):
	# pi0 invariant mass of two 3d photon momenta
	pi0mass2 = pi0Energy(gam1, gam2)**2-mag2(vecAdd(gam1, gam2))
	pi0mass2 = np.where(pi0mass2 >= 0, pi0mass2, 10**6)
	pi0mass = np.sqrt(pi0mass2)
	pi0mass = np.where(pi0mass > 100, -1000, pi0mass)
	return pi0mass

def getPhi(vec1):
	# azimuthal angle of one 3d vector
	return 180/np.pi*np.arctan2(vec1[1], vec1[0])

def getTheta(vec1):
	# polar angle of one 3d vector
	return 180/np.pi*np.arctan2(np.sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]), vec1[2])

def getEnergy(vec1, mass):
	# for taken 3d momenta p and mass m, return energy = sqrt(p**2 + m**2)
	return np.sqrt(mag2(vec1)+mass**2)

def rotateDCHitPosition(x, y, sec):
    ang = np.radians((sec-1) * sect_angle_coverage)
    x1_rot = y * np.sin(ang) + x * np.cos(ang)
    y1_rot = y * np.cos(ang) - x * np.sin(ang)
    return x1_rot, y1_rot

def determineSector(x, y):
    phi = getPhi([x, y])
    sector_cond = [(phi < 30) & (phi >= -30), (phi < 90) & (phi >= 30), (phi < 150) & (phi >= 90), (phi >= 150) | (phi < -150), (phi < -90) & (phi >= -150), (phi < -30) & (phi >= -90)]
    return np.select(sector_cond, [1, 2, 3, 4, 5, 6])


def thetaphifromhit(x, y, z):
    theta = getTheta([x, y, z])
    phi = getPhi([x, y, z])
    sector_cond = [(phi < 30) & (phi >= -30), (phi < 90) & (phi >= 30), (phi < 150) & (phi >= 90), (phi >= 150), (phi < -150), (phi < -90) & (phi >= -150), (phi < -30) & (phi >= -90)]
    return theta, phi + np.select(sector_cond, [0, -60, -120, -180, 180, 120, 60])

def e_DC_fiducial_cut_XY(x_rot, sec, region, minparams, maxparams):
    # if (pid==11) pid_ind =0
    # else if (pid==2212) pid_ind=1
    # else if (pid==211)  pid_ind=2
    # else if (pid==-211) pid_ind=3
    # else if (pid==321)  pid_ind=4
    # else if (pid==-321) pid_ind=5
    calc_min = np.array(minparams)[0, sec-1, region, 0] + np.array(minparams)[0, sec-1, region, 1] * x_rot;
    calc_max = np.array(maxparams)[0, sec-1, region, 0] + np.array(maxparams)[0, sec-1, region, 1] * x_rot;
    return calc_min, calc_max

def p_DC_fiducial_cut_XY(x_rot, sec, region, minparams, maxparams):
    # if (pid==11) pid_ind =0
    # else if (pid==2212) pid_ind=1
    # else if (pid==211)  pid_ind=2
    # else if (pid==-211) pid_ind=3
    # else if (pid==321)  pid_ind=4
    # else if (pid==-321) pid_ind=5
    calc_min = np.array(minparams)[1, sec-1, region, 0] + np.array(minparams)[1, sec-1, region, 1] * x_rot;
    calc_max = np.array(maxparams)[1, sec-1, region, 0] + np.array(maxparams)[1, sec-1, region, 1] * x_rot;
    return calc_min, calc_max

def p_DC_fiducial_cut_thetaphi(theta_DC, sec, region, minparams, maxparams):
    calc_phi_min = np.array(minparams)[1, sec-1, region, 0] + np.array(minparams)[1, sec-1, region, 1] * np.log(theta_DC) + np.array(minparams)[1, sec-1, region, 2] * theta_DC + np.array(minparams)[1, sec-1, region, 3] * theta_DC * theta_DC;
    calc_phi_max = np.array(maxparams)[1, sec-1, region, 0] + np.array(maxparams)[1, sec-1, region, 1] * np.log(theta_DC) + np.array(maxparams)[1, sec-1, region, 2] * theta_DC + np.array(maxparams)[1, sec-1, region, 3] * theta_DC * theta_DC;
    return calc_phi_min, calc_phi_max

def nu(xB, Q2, t, phi):
    return Q2/(2*M*xB)

def y(xB, Q2, t, phi):
    return nu(xB, Q2, t, phi)/10.604

def xi(xB, Q2, t, phi):
    return xB*(1-t/2/Q2)/(2-xB+xB*(-t/Q2))

def del2(xB, Q2, t, phi):
    return -t

def del2q2(xB, Q2, t, phi):
    return del2(xB, Q2, t, phi)/Q2

def eps(xB, Q2, t, phi):
    return 2*xB*M/np.sqrt(Q2)

def eps2(xB, Q2, t, phi):
    return eps(xB, Q2, t, phi)**2

def qeps2(xB, Q2, t, phi):
    return 1+eps2(xB, Q2, t, phi)

def sqeps2(xB, Q2, t, phi):
    return np.sqrt(qeps2(xB, Q2, t, phi))

def y1eps(xB, Q2, t, phi):
    return 1 - y(xB, Q2, t, phi) - y(xB, Q2, t, phi)*y(xB, Q2, t, phi)*eps2(xB, Q2, t, phi)/4

def tmin(xB, Q2, t, phi):
    return -Q2*(2*(1-xB)*(1-sqeps2(xB, Q2, t, phi))+eps2(xB, Q2, t, phi))/(4*xB*(1-xB)+eps2(xB, Q2, t, phi))

def tmax(xB, Q2, t, phi):
    return -Q2*(2*(1-xB)*(1+sqeps2(xB, Q2, t, phi))+eps2(xB, Q2, t, phi))/(4*xB*(1-xB)+eps2(xB, Q2, t, phi))

def W2(xB, Q2, t, phi):
    return M*M+2.0*M*nu(xB, Q2, t, phi)-Q2

def W(xB, Q2, t, phi):
    return np.sqrt(W2(xB, Q2, t, phi))

def tmin2(xB, Q2, t, phi):
    return -0.5*((Q2/xB-Q2)*(Q2/xB-np.sqrt((Q2/xB)**2+4*M*M*Q2))+2*M*M*Q2)/W2(xB, Q2, t, phi)

def tmax2(xB, Q2, t, phi):
    return -0.5*((Q2/xB-Q2)*(Q2/xB+np.sqrt((Q2/xB)**2+4*M*M*Q2))+2*M*M*Q2)/W2(xB, Q2, t, phi)

def Kfac2(xB, Q2, t, phi):
    return (-del2q2(xB, Q2, t, phi))*(1 - xB)*y1eps(xB, Q2, t, phi)*(1 - np.abs(tmin(xB, Q2, t, phi))/t)*(np.sqrt(1 + eps2(xB, Q2, t, phi)) + 
            ((4*xB*(1 - xB) + eps2(xB, Q2, t, phi))/(4*(1 - xB)))*(-(t - np.abs(tmin(xB, Q2, t, phi)))/Q2))
def Kfac(xB, Q2, t, phi):
    return np.sqrt(Kfac2(xB, Q2, t, phi))

def Jfac(xB, Q2, t, phi):
    return (1 - y(xB, Q2, t, phi) - y(xB, Q2, t, phi)*eps2(xB, Q2, t, phi)/2)*(1 + del2q2(xB, Q2, t, phi)) - (1 - xB)*(2 - y(xB, Q2, t, phi))*del2q2(xB, Q2, t, phi)

def P1(xB, Q2, t, phi):    
    return -(Jfac(xB, Q2, t, phi) + 2*Kfac(xB, Q2, t, phi)*np.cos(np.pi-np.radians(phi)))/(y(xB, Q2, t, phi)*(1 + eps2(xB, Q2, t, phi)))

def P2(xB, Q2, t, phi):
    return 1 + del2q2(xB, Q2, t, phi) - P1(xB, Q2, t, phi)

def printKMarray(xBarray, Q2array, tarray, phiarray, **kwargs):
    BHarray = []
    if isinstance(xBarray, pd.core.series.Series):
        xBarray = xBarray.to_numpy()
        Q2array = Q2array.to_numpy()
        tarray = tarray.to_numpy()
        phiarray = phiarray.to_numpy()
        
    for xB, Q2, t, phi in zip(xBarray, Q2array, tarray, phiarray):
        BHarray.append(printKM(xB, Q2, t, phi, **kwargs))
    return np.array(BHarray)

def printKM(xB, Q2, t, phi, frame = 'trento', pol = 0, mode = 5):
    phi = np.pi - phi
    pt1 = g.DataPoint(xB=xB, t=-t, Q2=Q2, phi=phi,
    				process='ep2epgamma', exptype='fixed target', frame =frame,
    				in1energy=10.604, in1charge=-1, in1polarization=pol)
    pt1.prepare()
    if mode == 0:
        try:
            return th_KM15.PreFacSigma(pt1)*pt1.P1P2*th_KM15.TBH2unp(pt1)
        except:
            print(xB, Q2, t, phi)
            return 0
    if mode == 1:
        try:
            return th_KM15.PreFacSigma(pt1)*th_KM15.TBH2unp(pt1)
        except:
            print(xB, Q2, t, phi)
            return 0
    if mode == 2:
        try:
            return th_KM15.PreFacSigma(pt1)*th_KM15.TINTunp(pt1)
        except:
            print(xB, Q2, t, phi)
            return 0
    if mode == 3:
        try:
            return th_KM15.PreFacSigma(pt1)*th_KM15.TDVCS2unp(pt1)
        except:
            print(xB, Q2, t, phi)
            return 0
    if mode == 4:
        try:
            return th_KM15.PreFacSigma(pt1)*(th_KM15.TINTunp(pt1)+th_KM15.TDVCS2unp(pt1))
        except:
            print(xB, Q2, t, phi)
            return 0
    if mode == 5:
        try:
            return th_KM15.PreFacSigma(pt1)*(th_KM15.TBH2unp(pt1)+th_KM15.TINTunp(pt1)+th_KM15.TDVCS2unp(pt1))
        except:
            print(xB, Q2, t, phi)
            return 0
    else:
        try:
            return th_KM15.XS(pt1)
        except:
            print(xB, Q2, t, phi)
            return 0

def printVGGarray(xBarray, Q2array, tarray, phiarray, **kwargs):
    VGGarray = []
    if isinstance(xBarray, pd.core.series.Series):
        xBarray = xBarray.to_numpy()
        Q2array = Q2array.to_numpy()
        tarray = tarray.to_numpy()
        phiarray = phiarray.to_numpy()
        
    for xB, Q2, t, phi in zip(xBarray, Q2array, tarray, phiarray):
        VGGarray.append(printVGG(xB, Q2, t, phi, **kwargs))
    return np.array(VGGarray)

def printVGG(xB, Q2, t, phi, globalfit = True, pol = 0, local = False):
    my_env = os.environ.copy()
    path = "/home/sangbaek/printDVCSBH/"
    if local:
        path = "/Users/sangbaek/CLAS12/dvcs/print/"
    my_env["PATH"] = "{}:".format(path) + my_env["PATH"]
    my_env["CLASDVCS_PDF"] = "{}".format(path)
    if globalfit:
        dstot = subprocess.check_output(['{}/dvcsgen'.format(path), '--beam', '10.604', '--x', str(xB), str(xB), '--q2', str(Q2), str(Q2),'--t', str(t), str(t), '--bh', '3', '--phi', str(phi), '--gpd', '101', '--globalfit'], env = my_env)
    else:
        dstot = subprocess.check_output(['{}/dvcsgen'.format(path), '--beam', '10.604', '--x', str(xB), str(xB), '--q2', str(Q2), str(Q2),'--t', str(t), str(t), '--bh', '3', '--phi', str(phi), '--gpd', '101'], env = my_env)
    try:
        if pol == 0:
            i = 0
        if pol == 1:
            i = 2
        if pol == -1:
            i = 1
        dstot = float(dstot.splitlines()[-1-i].decode("utf-8"))
        return dstot
    except:
        print(xB, Q2, t, phi)
        return 0

def printBHarray(xBarray, Q2array, tarray, phiarray, **kwargs):
    BHarray = []
    if isinstance(xBarray, pd.core.series.Series):
        xBarray = xBarray.to_numpy()
        Q2array = Q2array.to_numpy()
        tarray = tarray.to_numpy()
        phiarray = phiarray.to_numpy()
        
    for xB, Q2, t, phi in zip(xBarray, Q2array, tarray, phiarray):
        BHarray.append(printBHonly(xB, Q2, t, phi, **kwargs))
    return np.array(BHarray)

def printBHonly(xB, Q2, t, phi, globalfit = True, local = False):
    path = "/home/sangbaek/printDVCSBH/"
    if local:
        path = "/Users/sangbaek/CLAS12/dvcs/print/"
    if globalfit:
        dstot = subprocess.check_output(['{}/dvcsgen'.format(path), '--beam', '10.604', '--x', str(xB), str(xB), '--q2', str(Q2), str(Q2),'--t', str(t), str(t), '--bh', '1', '--phi', str(phi), '--globalfit'])
    else:
        dstot = subprocess.check_output(['{}/dvcsgen'.format(path), '--beam', '10.604', '--x', str(xB), str(xB), '--q2', str(Q2), str(Q2),'--t', str(t), str(t), '--bh', '1', '--phi', str(phi)])
    try:
        dstot = float(dstot.splitlines()[-1].decode("utf-8"))
        return dstot
    except:
        print(xB, Q2, t, phi)
        return 0

def nphistmean(hist, bins):
    s=0
    for i in range(len(hist)):
        s += hist[i] * ((bins[i] + bins[i+1]) / 2) 
    mean = s / np.sum(hist)
    return mean

def createBinEdges(binCenters):
    start = binCenters[0] - np.diff(binCenters)[0]/2
    end = binCenters[-1] + np.diff(binCenters)[-1]/2
    middle = binCenters[:-1] + np.diff(binCenters)/2
    return np.array([start, *middle, end])

def makeReduced(df):
    columns_needed = ["helicity", "polarity", "config", "beamCurrent", "xB", "Q2", "t1", "phi1"]
    return df.loc[:, columns_needed]

def readReduced(parent, jobNum, polarity, beamCurrent):
    df = pd.read_pickle(parent + "{}.pkl".format(jobNum))
    df.loc[:, "polarity"] = polarity
    df.loc[:, "beamCurrent"] = beamCurrent
    columns_needed = ["polarity", "config", "beamCurrent", "xB", "Q2", "t1", "phi1"]
    return df.loc[:, columns_needed]

def divideHist(df1, df2, threshold = 0):
	return np.divide(df1, df2, where = (df2!=0) & (df1>threshold), out = np.zeros(df2.shape, dtype = float))

def inverseHist(df1):
	return np.divide(np.ones(df1.shape), df1, where = df1!=0, out = np.zeros_like(df1))

def binVolumes(xBbin, Q2bin, tbin, finehist, k=0):
	xBbins  = collection_xBbins[k]
	Q2bins  = collection_Q2bins[k]
	tbins   = collection_tbins [k]
	phibins = collection_phibins[k]
	fineVols = []
	for phibin in range(len(phibins)-1):
		fineVol = finehist[6*xBbin:6*(xBbin+1), 6*Q2bin:6*(Q2bin+1), 6*tbin:6*(tbin+1), 6*phibin:6*(phibin+1)].flatten()
		fineVols.append(fineVol)
	return (np.sum(fineVols, axis = 1)/6**4)*np.diff(np.radians(phibins))*np.diff(xBbins)[xBbin]*np.diff(Q2bins)[Q2bin]*np.diff(tbins)[tbin]

def getCFFarrays(xB, Q2, t, phi):
    phigd = np.pi - np.radians(phi)
    if isinstance(xB, np.ndarray):
        H1_RE, H1_IM, E1_RE, E1_IM, H1T_RE, H1T_IM, E1T_RE, E1T_IM = ([], [], [], [], [], [], [], [])
        for i in range(len(xB)):
            pt1 = g.DataPoint(xB=xB[i], t=-t[i], Q2=Q2[i], phi=phigd[i],
                          process='ep2epgamma', exptype='fixed target', frame ='trento',
                          in1energy=10.604, in1charge=-1, in1polarization=0)
            H1_RE.append(th_KM15.ReH(pt1))
            H1_IM.append(th_KM15.ImH(pt1))
            E1_RE.append(th_KM15.ReE(pt1))
            E1_IM.append(th_KM15.ImE(pt1))
            H1T_RE.append(th_KM15.ReHt(pt1))
            H1T_IM.append(th_KM15.ImHt(pt1))
            E1T_RE.append(th_KM15.ReEt(pt1))
            E1T_IM.append(th_KM15.ImEt(pt1))
        H1_RE, H1_IM, E1_RE, E1_IM, H1T_RE, H1T_IM, E1T_RE, E1T_IM = (np.array(H1_RE), np.array(H1_IM), np.array(E1_RE), np.array(E1_IM), np.array(H1T_RE), np.array(H1T_IM), np.array(E1T_RE), np.array(E1T_IM))
    else:
        pt1 = g.DataPoint(xB=xB, t=-t, Q2=Q2, phi=phigd,
                      process='ep2epgamma', exptype='fixed target', frame ='trento',
                      in1energy=10.604, in1charge=-1, in1polarization=0)
        H1_RE, H1_IM, E1_RE, E1_IM, H1T_RE, H1T_IM, E1T_RE, E1T_IM = (th_KM15.ReH(pt1), th_KM15.ImH(pt1), th_KM15.ReE(pt1), th_KM15.ImE(pt1), th_KM15.ReHt(pt1), th_KM15.ImHt(pt1), th_KM15.ReEt(pt1), th_KM15.ImEt(pt1))
    return H1_RE, H1_IM, E1_RE, E1_IM, H1T_RE, H1T_IM, E1T_RE, E1T_IM

def getBHDVCS(xB, Q2, t, phi, Ed = 10.604, mode = 1):
    coeff = 10**9*hc2*alpha**3
    del2  = -t
    Phi_gb = np.pi-np.radians(phi)
    phigd = Phi_gb
    nu  = Q2/(2*M*xB)
    yb = nu/Ed
    eps = 2*xB*M/np.sqrt(Q2)
    eps2=eps*eps
    qeps2=1 + eps2
    sqeps2=np.sqrt(qeps2)
    ds = 2*np.pi*0.001*coeff*(xB*yb**2/(16*np.pi**2*Q2**2))/sqeps2
    del2min=-Q2*(2*(1-xB)*(1-sqeps2)+eps2)
    del2min=del2min/(4*xB*(1-xB)+eps2)
    tau  = del2/(4*M**2)
    taum1=1-tau
    xtau=xB*xB/tau
    del2q2=del2/Q2
    del2q4=del2q2*del2q2
    del2q2m1=1-del2q2
    phip = phigd 
    GE_p, GM_p, GE_n, GM_n, dipol = ffs( - t)
    delm = del2/(2*M)**2
    F1 = (GE_p - delm*GM_p)/(1-delm)
    F2 = (GM_p - GE_p)/(1-delm)
    y1eps=1 - yb - yb*yb*eps2/4
    sqy1eps=np.sqrt(y1eps)
    Kfac = np.sqrt((-del2q2)*(1 - xB)*y1eps*
        (1 - del2min/del2)*(np.sqrt(1 + eps2) + 
        ((4*xB*(1 - xB) + eps2)/(4*(1 - xB)))*((del2 - del2min)/Q2)))
    Jfac = (1 - yb - yb*eps2/2)*(1 + del2q2) - (1 - xB)*(2 - yb)*del2q2
    P1 = -(Jfac + 2*Kfac*np.cos(phip))/(yb*(1 + eps2))
    P2 = 1 + del2q2 - P1
    F1_M_F = F1 + tau*F2
    F1_M_F2 = F1**2 - tau*F2**2
    F1_P_F = F1 + F2
    F1_P_F2 = F1_P_F*F1_P_F
    c01_BH = 8*Kfac**2*((2 + 3*eps2)*Q2*F1_M_F2/del2 + 2*xB**2*F1_P_F2)
    c02_BH = (2 - yb)**2*((2 + eps2)*F1_M_F2*((2*xB*M)**2*(1 + del2q2)**2/del2 + 
        4*(1 - xB)*(1 + xB*del2q2)) +
        4*xB**2*F1_P_F2*(xB + (1 - xB + eps2/2)*(del2q2m1)**2 - xB*(1 - 2*xB)*(del2q2)**2))
    c03_BH = 8*(1 + eps2)*(1 - yb - yb*yb*eps2/4)*(2*eps2*(1 - del2/(4*M**2))*F1_M_F2 - xB**2*(del2q2m1)**2*F1_P_F2)
    c0_BH = c01_BH + c02_BH + c03_BH
    c1_BH = 8*Kfac*(2 - yb)*(F1_M_F2*(4*(xB*M)**2/del2 - 2*xB - eps2) + F1_P_F2*2*xB**2*(1 - (1 - 2*xB)*del2q2))
    c2_BH = 8*(xB*Kfac)**2*(F1_M_F2*4*M**2/del2 + 2*F1_P_F2)
    BHfact = ds/((xB*yb*(1 + eps2))**2*del2*P1*P2)
    hc0BH=c0_BH*BHfact
    hc1BH=c1_BH*BHfact
    hc2BH=c2_BH*BHfact
    if mode == 0:
        return (c0_BH-c1_BH+c2_BH)*BHfact*P1*P2
    if mode == 1:
        return hc0BH +hc1BH*np.cos(Phi_gb)+hc2BH*np.cos(2*Phi_gb)


    H1_RE, H1_IM, E1_RE, E1_IM, H1T_RE, H1T_IM, E1T_RE, E1T_IM = getCFFarrays(xB, Q2, t, phi)

    deldel    = 1 - del2min/del2
    deldel_sq = np.sqrt(deldel)
    del2m2    = -del2/M**2
    del2m4    = -del2m2/4
    delm2_sq  = np.sqrt(del2m2)
    cy2   = 2 - 2*yb + yb**2

    C_I_re = F1*H1_RE + xB/(2-xB)*(F1+F2)*H1T_RE - del2m4*F2*E1_RE
    C_I_im = F1*H1_IM + xB/(2-xB)*(F1+F2)*H1T_IM - del2m4*F2*E1_IM
    RE2    = xB/(2-xB)*(H1_RE + E1_RE) + H1T_RE
    C_I_re_eff = -xB*C_I_re
    C_I_im_eff = -xB*C_I_im
    b1 = (2 - xB)*(1 - yb) 
    b1= b1 + (2-yb)**2/(1-yb)*Kfac*Kfac/del2*Q2
    b2 = (1 - yb)*xB*(F1 + F2)
    c0_I = -8*(2 - yb)*( b1*C_I_re - b2*RE2 )
    c1_I = -8*cy2*C_I_re
    s1_I =  8*yb*(2 - yb)*C_I_im
    c2_I = -16*((2 - yb)/(2 - xB))*C_I_re_eff
    s2_I =  16*(yb/(2 - xB))*C_I_im_eff
    Intfac=ds/(xB*yb**3*P1*P2*(del2))
    hs2Iunp=Kfac*Kfac*s2_I*Intfac
    hs1Iunp=Kfac*s1_I*Intfac
    hc2Iunp=Kfac*Kfac*c2_I*Intfac
    hc1Iunp=Kfac*c1_I*Intfac
    hc0Iunp=del2/Q2*c0_I*Intfac

    a1 = H1_RE**2 + H1_IM**2 + H1T_RE**2 + H1T_IM**2
    a2 = 2.*( H1_RE*E1_RE  +  H1_IM*E1_IM +  H1T_RE*E1T_RE + H1T_IM*E1T_IM )
    a3 =  E1_RE**2 +  E1_IM**2
    a4 = E1T_RE**2 + E1T_IM**2 
    C_DVCS = ( 4*(1-xB)*a1 - a2*xB**2 - (xB**2 + (2-xB)**2*del2m4)*a3 - xB**2*del2m4*a4 )/(2 - xB)**2
    C_DVCS_eff = -xB*C_DVCS
    c0_DVCS = 2*cy2*C_DVCS
    c1_DVCS = 8*((2 - yb)/(2 - xB))*C_DVCS_eff
    DVCSfac=ds/(yb**2*Q2)
    hc0dvcs=c0_DVCS*DVCSfac
    hc1dvcs=Kfac*c1_DVCS*DVCSfac

    if mode == 2:
        return hc0Iunp+hc1Iunp*np.cos(Phi_gb)+hc2Iunp*np.cos(2*Phi_gb)
    if mode == 3:
        return hc0dvcs+hc1dvcs*np.cos(Phi_gb)
    if mode == 4:
        return hc0dvcs+hc1dvcs*np.cos(Phi_gb) + hc0Iunp+hc1Iunp*np.cos(Phi_gb)+hc2Iunp*np.cos(2*Phi_gb)
    if mode == 5:
        return hc0BH +hc1BH*np.cos(Phi_gb)+hc2BH*np.cos(2*Phi_gb)+hc0dvcs+hc1dvcs*np.cos(Phi_gb) + hc0Iunp+hc1Iunp*np.cos(Phi_gb)+hc2Iunp*np.cos(2*Phi_gb)

def ffs(del2):
    a0 = [0.239163298067, 0.264142994136, 0.048919981379, 0.257758326959]
    a1 = [-1.109858574410, -1.095306122120, -0.064525053912, -1.079540642058]
    a2 = [1.444380813060, 1.218553781780, -0.240825897382, 1.182183812195]
    a3 = [0.479569465603, 0.661136493537, 0.392108744873, 0.711015085833]
    a4 = [-2.286894741870, -1.405678925030, 0.300445258602, -1.348080936796]
    a5 = [1.126632984980, -1.356418438880, -0.661888687179, -1.662444025208]
    a6 = [1.250619843540, 1.447029155340, -0.175639769687, 2.624354426029]
    a7 = [-3.631020471590, 4.235669735900, 0.624691724461, 1.751234494568]
    a8 = [4.082217023790, -5.334045653410, -0.077684299367, -4.922300878888]
    a9 = [0.504097346499, -2.916300520960, -0.236003975259, 3.197892727312]
    a10 = [-5.085120460510, 8.707403067570, 0.090401973470, -0.712072389946]
    a11 = [3.967742543950, -5.706999943750, 0.0, 0.0]
    a12 = [-0.981529071103, 1.280814375890, 0.0, 0.0]

    Mv = 0.843
    kp = 1.79285
    mu_n = -1.91
    tcut = 4*0.13957**2
    t0 = -0.7


    dipol = 1/(1 - del2/Mv**2)**2

    z = (np.sqrt(tcut-del2)-np.sqrt(tcut-t0))/(np.sqrt(tcut-del2)+np.sqrt(tcut-t0)) 

    GE_p = (a0[0]*z**0+a1[0]*z**1+a2[0]*z**2+a3[0]*z**3+a4[0]*z**4+a5[0]*z**5+a6[0]*z**6+a7[0]*z**7+a8[0]*z**8+a9[0]*z**9+a10[0]*z**10+a11[0]*z**11+a12[0]*z**12)
    GM_p = (1 + kp)*(a0[1]*z**0+a1[1]*z**1+a2[1]*z**2+a3[1]*z**3+a4[1]*z**4+a5[1]*z**5+a6[1]*z**6+a7[1]*z**7+a8[1]*z**8+a9[1]*z**9+a10[1]*z**10+a11[1]*z**11+a12[1]*z**12)
    GE_n = (a0[2]*z**0+a1[2]*z**1+a2[2]*z**2+a3[2]*z**3+a4[2]*z**4+a5[2]*z**5+a6[2]*z**6+a7[2]*z**7+a8[2]*z**8+a9[2]*z**9+a10[2]*z**10+a11[2]*z**11+a12[2]*z**12)
    GM_n = mu_n*(a0[3]*z**0+a1[3]*z**1+a2[3]*z**2+a3[3]*z**3+a4[3]*z**4+a5[3]*z**5+a6[3]*z**6+a7[3]*z**7+a8[3]*z**8+a9[3]*z**9+a10[3]*z**10+a11[3]*z**11+a12[3]*z**12)

    dipol = 1/(1 - del2/Mv**2)**2

    delm = del2/(2*M)**2

    F1new = (GE_p - delm*GM_p)/(1-delm)
    F2new = (GM_p - GE_p)/(1-delm)         
    F1old = (dipol - delm*(1+kp)*dipol)/(1-delm)
    F2old = ((1+kp)*dipol - dipol)/(1-delm)         

    return GE_p, GM_p, GE_n, GM_n, dipol