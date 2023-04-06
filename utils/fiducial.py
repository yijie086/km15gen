from utils.const import *
from utils.physics import *

def electronFiducial(df_electronRec, pol = "inbending", mc = False, fidlevel = 'mid'):
	df_electronRec.loc[:, "EFid"] = 1

	# #PCAL dead wires
	# exclusion1_1 = (df_electronRec.EcalW1 > 74) & (df_electronRec.EcalW1 < 79.8)
	# exclusion1_2 = (df_electronRec.EcalW1 > 83.6) & (df_electronRec.EcalW1 < 92.2)
	# exclusion1_3 = (df_electronRec.EcalW1 > 212.5) & (df_electronRec.EcalW1 < 230)
	# exclusion1 = exclusion1_1 | exclusion1_2 | exclusion1_3
	# df_electronRec.loc[(df_electronRec.Esector == 1) & exclusion1, "EFid"] = 0
	# exclusion2_1 = (df_electronRec.EcalW1 < 14)
	# exclusion2_2 = (df_electronRec.EcalU1 > 111.2) & (df_electronRec.EcalU1 < 119.3)
	# exclusion2_3 = (df_electronRec.EcalV1 > 113) & (df_electronRec.EcalV1 < 118.7)
	# exclusion2 = exclusion2_1 | exclusion2_2 | exclusion2_3
	# df_electronRec.loc[(df_electronRec.Esector == 2) & exclusion2, "EFid"] = 0
	# exclusion3 = df_electronRec.EcalW1 < 14
	# df_electronRec.loc[(df_electronRec.Esector == 3) & exclusion3, "EFid"] = 0
	# exclusion4_1 = (df_electronRec.EcalV1 < 14)
	# exclusion4_2 = (df_electronRec.EcalV1 > 229.4) & (df_electronRec.EcalV1 < 240.7)
	# exclusion4_3 = (df_electronRec.EcalW1 > 135) & (df_electronRec.EcalW1 < 150)
	# exclusion4 = exclusion4_1 | exclusion4_2 | exclusion4_3
	# df_electronRec.loc[(df_electronRec.Esector == 4) & exclusion4, "EFid"] = 0
	# exclusion6 = (df_electronRec.EcalW1 > 170) & (df_electronRec.EcalW1 < 192)
	# df_electronRec.loc[(df_electronRec.Esector == 6) & exclusion6, "EFid"] = 0

	# passElectronTrackQualityCut (pass)
	sector_cond = [df_electronRec.Esector ==1, df_electronRec.Esector ==2, df_electronRec.Esector ==3, df_electronRec.Esector ==4, df_electronRec.Esector ==5, df_electronRec.Esector ==6]

	# passElectronSamplingFractionCut
	ecal_e_sampl_mu_0 = np.select(sector_cond, ecal_e_sampl_mu[0])
	ecal_e_sampl_mu_1 = np.select(sector_cond, ecal_e_sampl_mu[1])
	ecal_e_sampl_mu_2 = np.select(sector_cond, ecal_e_sampl_mu[2])
	ecal_e_sampl_sigm_0 = np.select(sector_cond, ecal_e_sampl_sigm[0])
	ecal_e_sampl_sigm_1 = np.select(sector_cond, ecal_e_sampl_sigm[1])
	ecal_e_sampl_sigm_2 = np.select(sector_cond, ecal_e_sampl_sigm[2])

	if mc:
		ecal_e_sampl_mu_0 = np.select(sector_cond, ecal_e_sampl_mu_mc[0])
		ecal_e_sampl_mu_1 = np.select(sector_cond, ecal_e_sampl_mu_mc[1])
		ecal_e_sampl_mu_2 = np.select(sector_cond, ecal_e_sampl_mu_mc[2])
		ecal_e_sampl_sigm_0 = np.select(sector_cond, ecal_e_sampl_sigm_mc[0])
		ecal_e_sampl_sigm_1 = np.select(sector_cond, ecal_e_sampl_sigm_mc[1])
		ecal_e_sampl_sigm_2 = np.select(sector_cond, ecal_e_sampl_sigm_mc[2])
	mean = ecal_e_sampl_mu_0 + ecal_e_sampl_mu_1/1000*pow(df_electronRec.Ep-ecal_e_sampl_mu_2,2)
	sigma = ecal_e_sampl_sigm_0 + ecal_e_sampl_sigm_1/(10*(df_electronRec.Ep-ecal_e_sampl_sigm_2))
	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - e_sampl_sigma_range*sigma, "EFid"]  = 0
		df_electronRec.loc[df_electronRec.ESamplFrac > mean + e_sampl_sigma_range*sigma, "EFid"]  = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - (e_sampl_sigma_range-0.5)*sigma, "EFid"]  = 0
		df_electronRec.loc[df_electronRec.ESamplFrac > mean + (e_sampl_sigma_range-0.5)*sigma, "EFid"]  = 0


	#passElectronNpheCut
	df_electronRec.loc[df_electronRec.Enphe <= min_nphe, "EFid"] = 0

	#passElectronVertexCut
	if pol == 'inbending':
		min_vz = vz_min_inb
		max_vz = vz_max_inb
	if pol == 'outbending':
		min_vz = vz_min_outb
		max_vz = vz_max_outb
	df_electronRec.loc[df_electronRec.Evz <= min_vz, "EFid"] = 0
	df_electronRec.loc[df_electronRec.Evz >= max_vz, "EFid"] = 0

	# passElectronPCALFiducialCut
	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.EcalV1 <= min_v, "EFid"] = 0
		df_electronRec.loc[df_electronRec.EcalW1 <= min_w, "EFid"] = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.EcalV1 <= min_v+5, "EFid"] = 0
		df_electronRec.loc[df_electronRec.EcalW1 <= min_w+5, "EFid"] = 0

	#passElectronPCALEdepCut
	df_electronRec.loc[df_electronRec.Eedep1 <= min_pcal_dep, "EFid"] = 0

	#passElectronDCR1
	if pol == 'inbending':
		minparams = e_dc_minparams_in
		maxparams = e_dc_maxparams_in
	if pol == 'outbending':
		minparams = e_dc_minparams_out
		maxparams = e_dc_maxparams_out

	dcsec = determineSector(df_electronRec.EDc1Hitx, df_electronRec.EDc1Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc1Hitx, df_electronRec.EDc1Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 0, minparams, maxparams)
	df_electronRec.loc[y_rot<=calc_min, "EFid"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid"] = 0
	#passElectronDCR2
	dcsec = determineSector(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 1, minparams, maxparams)
	df_electronRec.loc[y_rot<=calc_min, "EFid"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid"] = 0

	#passElectronDCR3
	dcsec = determineSector(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 2, minparams, maxparams)
	df_electronRec.loc[y_rot<=calc_min, "EFid"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid"] = 0

	# #passElectronAntiPionCut
	df_electronRec.loc[(df_electronRec.Ep>4.5)&(-df_electronRec.Eedep1/df_electronRec.Ep + anti_pion_threshold > df_electronRec.Eedep2/df_electronRec.Ep), "EFid"] = 0
	return df_electronRec.loc[df_electronRec.EFid==1, :]

def gammaFiducial(df_gammaRec):
	df_gammaRec.loc[:, "GFid"] = 1
	#passGammaPCALFiducialCut
	df_gammaRec.loc[(df_gammaRec.GcalV1 <= g_min_v) & (df_gammaRec.Gsector<7), "GFid"] = 0
	df_gammaRec.loc[(df_gammaRec.GcalW1 <= g_min_w) & (df_gammaRec.Gsector<7), "GFid"] = 0
	#passGammaBetaCut
	df_gammaRec.loc[df_gammaRec.Gbeta <= min_Gbeta, "GFid"] = 0
	df_gammaRec.loc[df_gammaRec.Gbeta >= max_Gbeta, "GFid"] = 0

	return df_gammaRec.loc[df_gammaRec.GFid==1, :]

def protonFiducial(df_protonRec, pol = 'inbending'):
	df_protonRec.loc[:, "PFid"] = 1

	dcsec = determineSector(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity)
	if pol == 'inbending':
		minparams = p_dc_minparams_in
		maxparams = p_dc_maxparams_in

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity, df_protonRec.PDc1Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 0, minparams, maxparams)
		df_protonRec.loc[(phi_DC<=phi_DC_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(phi_DC>=phi_DC_max) & (df_protonRec.Psector<7), "PFid"] = 0

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity, df_protonRec.PDc2Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 1, minparams, maxparams)
		df_protonRec.loc[(phi_DC<=phi_DC_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(phi_DC>=phi_DC_max) & (df_protonRec.Psector<7), "PFid"] = 0

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity, df_protonRec.PDc3Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 2, minparams, maxparams)
		df_protonRec.loc[(phi_DC<=phi_DC_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(phi_DC>=phi_DC_max) & (df_protonRec.Psector<7), "PFid"] = 0

	if pol == 'outbending':
		minparams = p_dc_minparams_out
		maxparams = p_dc_maxparams_out

		dcsec = determineSector(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 0, minparams, maxparams)
		df_protonRec.loc[(y_rot<=calc_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(y_rot>=calc_max) & (df_protonRec.Psector<7), "PFid"] = 0
		#passElectronDCR2
		dcsec = determineSector(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 1, minparams, maxparams)
		df_protonRec.loc[(y_rot<=calc_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(y_rot>=calc_max) & (df_protonRec.Psector<7), "PFid"] = 0

		#passElectronDCR3
		dcsec = determineSector(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 2, minparams, maxparams)
		df_protonRec.loc[(y_rot<=calc_min) & (df_protonRec.Psector<7), "PFid"] = 0
		df_protonRec.loc[(y_rot>=calc_max) & (df_protonRec.Psector<7), "PFid"] = 0

	return df_protonRec.loc[df_protonRec.PFid==1, :]


def electronFiducialCounting(df_electronRec, pol = "inbending", mc = False, fidlevel = 'mid'):
	df_electronRec.loc[:, "EFid_dw"] = 1
	df_electronRec.loc[:, "EFid_sf"] = 1
	df_electronRec.loc[:, "EFid_vz"] = 1
	df_electronRec.loc[:, "EFid_pcal1"] = 1
	df_electronRec.loc[:, "EFid_vz"] = 1
	df_electronRec.loc[:, "EFid_edep"] = 1
	df_electronRec.loc[:, "EFid_dc"] = 1
	df_electronRec.loc[:, "EFid_ap"] = 1

	#PCAL dead wires
	exclusion1_1 = (df_electronRec.EcalW1 > 74) & (df_electronRec.EcalW1 < 79.8)
	exclusion1_2 = (df_electronRec.EcalW1 > 83.6) & (df_electronRec.EcalW1 < 92.2)
	exclusion1_3 = (df_electronRec.EcalW1 > 212.5) & (df_electronRec.EcalW1 < 230)
	exclusion1 = exclusion1_1 | exclusion1_2 | exclusion1_3
	df_electronRec.loc[(df_electronRec.Esector == 1) & exclusion1, "EFid_dw"] = 0
	exclusion2_1 = (df_electronRec.EcalW1 < 14)
	exclusion2_2 = (df_electronRec.EcalU1 > 111.2) & (df_electronRec.EcalU1 < 119.3)
	exclusion2_3 = (df_electronRec.EcalV1 > 113) & (df_electronRec.EcalV1 < 118.7)
	exclusion2 = exclusion2_1 | exclusion2_2 | exclusion2_3
	df_electronRec.loc[(df_electronRec.Esector == 2) & exclusion2, "EFid_dw"] = 0
	exclusion3 = df_electronRec.EcalW1 < 14
	df_electronRec.loc[(df_electronRec.Esector == 3) & exclusion3, "EFid_dw"] = 0
	exclusion4_1 = (df_electronRec.EcalV1 < 14)
	exclusion4_2 = (df_electronRec.EcalV1 > 229.4) & (df_electronRec.EcalV1 < 240.7)
	exclusion4_3 = (df_electronRec.EcalW1 > 135) & (df_electronRec.EcalW1 < 150)
	exclusion4 = exclusion4_1 | exclusion4_2 | exclusion4_3
	df_electronRec.loc[(df_electronRec.Esector == 4) & exclusion4, "EFid_dw"] = 0
	exclusion6 = (df_electronRec.EcalW1 > 170) & (df_electronRec.EcalW1 < 192)
	df_electronRec.loc[(df_electronRec.Esector == 6) & exclusion6, "EFid_dw"] = 0

	# passElectronTrackQualityCut (pass)
	sector_cond = [df_electronRec.Esector ==1, df_electronRec.Esector ==2, df_electronRec.Esector ==3, df_electronRec.Esector ==4, df_electronRec.Esector ==5, df_electronRec.Esector ==6]

	# passElectronSamplingFractionCut
	ecal_e_sampl_mu_0 = np.select(sector_cond, ecal_e_sampl_mu[0])
	ecal_e_sampl_mu_1 = np.select(sector_cond, ecal_e_sampl_mu[1])
	ecal_e_sampl_mu_2 = np.select(sector_cond, ecal_e_sampl_mu[2])
	ecal_e_sampl_sigm_0 = np.select(sector_cond, ecal_e_sampl_sigm[0])
	ecal_e_sampl_sigm_1 = np.select(sector_cond, ecal_e_sampl_sigm[1])
	ecal_e_sampl_sigm_2 = np.select(sector_cond, ecal_e_sampl_sigm[2])

	if mc:
		ecal_e_sampl_mu_0 = np.select(sector_cond, ecal_e_sampl_mu_mc[0])
		ecal_e_sampl_mu_1 = np.select(sector_cond, ecal_e_sampl_mu_mc[1])
		ecal_e_sampl_mu_2 = np.select(sector_cond, ecal_e_sampl_mu_mc[2])
		ecal_e_sampl_sigm_0 = np.select(sector_cond, ecal_e_sampl_sigm_mc[0])
		ecal_e_sampl_sigm_1 = np.select(sector_cond, ecal_e_sampl_sigm_mc[1])
		ecal_e_sampl_sigm_2 = np.select(sector_cond, ecal_e_sampl_sigm_mc[2])
	mean = ecal_e_sampl_mu_0 + ecal_e_sampl_mu_1/1000*pow(df_electronRec.Ep-ecal_e_sampl_mu_2,2)
	sigma = ecal_e_sampl_sigm_0 + ecal_e_sampl_sigm_1/(10*(df_electronRec.Ep-ecal_e_sampl_sigm_2))
	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - e_sampl_sigma_range*sigma, "EFid_sf"]  = 0
		df_electronRec.loc[df_electronRec.ESamplFrac > mean + e_sampl_sigma_range*sigma, "EFid_sf"]  = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.ESamplFrac < mean - (e_sampl_sigma_range-0.5)*sigma, "EFid_sf"]  = 0
		df_electronRec.loc[df_electronRec.ESamplFrac > mean + (e_sampl_sigma_range-0.5)*sigma, "EFid_sf"]  = 0


	# #passElectronNpheCut
	# df_electronRec.loc[df_electronRec.Enphe <= min_nphe, "EFid"] = 0

	#passElectronVertexCut
	if pol == 'inbending':
		min_vz = vz_min_inb
		max_vz = vz_max_inb
	if pol == 'outbending':
		min_vz = vz_min_outb
		max_vz = vz_max_outb
	df_electronRec.loc[df_electronRec.Evz <= min_vz, "EFid_vz"] = 0
	df_electronRec.loc[df_electronRec.Evz >= max_vz, "EFid_vz"] = 0

	# passElectronPCALFiducialCut
	if fidlevel == 'mid':
		df_electronRec.loc[df_electronRec.EcalV1 <= min_v, "EFid_pcal1"] = 0
		df_electronRec.loc[df_electronRec.EcalW1 <= min_w, "EFid_pcal1"] = 0
	elif fidlevel == 'tight':
		df_electronRec.loc[df_electronRec.EcalV1 <= min_v+10, "EFid_pcal1"] = 0
		df_electronRec.loc[df_electronRec.EcalW1 <= min_w+10, "EFid_pcal1"] = 0

	#passElectronPCALEdepCut
	df_electronRec.loc[df_electronRec.Eedep1 <= min_pcal_dep, "EFid_edep"] = 0

	#passElectronDCR1
	if pol == 'inbending':
		minparams = e_dc_minparams_in
		maxparams = e_dc_maxparams_in
	if pol == 'outbending':
		minparams = e_dc_minparams_out
		maxparams = e_dc_maxparams_out

	dcsec = determineSector(df_electronRec.EDc1Hitx, df_electronRec.EDc1Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc1Hitx, df_electronRec.EDc1Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 0, minparams, maxparams)
	df_electronRec.loc[y_rot<=calc_min, "EFid_dc"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid_dc"] = 0
	#passElectronDCR2
	dcsec = determineSector(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc2Hitx, df_electronRec.EDc2Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 1, minparams, maxparams)
	df_electronRec.loc[y_rot<=calc_min, "EFid_dc"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid_dc"] = 0

	#passElectronDCR3
	dcsec = determineSector(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity)
	x_rot, y_rot = rotateDCHitPosition(df_electronRec.EDc3Hitx, df_electronRec.EDc3Hity, dcsec)
	calc_min, calc_max = e_DC_fiducial_cut_XY(x_rot, dcsec, 2, minparams, maxparams)
	df_electronRec.loc[y_rot<=calc_min, "EFid_dc"] = 0
	df_electronRec.loc[y_rot>=calc_max, "EFid_dc"] = 0

	# #passElectronAntiPionCut
	df_electronRec.loc[(df_electronRec.Ep>4.5)&(-df_electronRec.Eedep1/df_electronRec.Ep + anti_pion_threshold > df_electronRec.Eedep2/df_electronRec.Ep), "EFid_ap"] = 0

	df_electronRec.loc[:, "EFid_pcal"] = df_electronRec.loc[:, "EFid_pcal1"] * df_electronRec.loc[:, "EFid_dw"]
	df_electronRec.loc[:, "EFid"] = df_electronRec.EFid_dw*df_electronRec.EFid_sf*df_electronRec.EFid_vz*df_electronRec.EFid_pcal1*df_electronRec.EFid_vz*df_electronRec.EFid_edep*df_electronRec.EFid_dc*df_electronRec.EFid_ap
	return df_electronRec

def gammaFiducialCounting(df_gammaRec):
	df_gammaRec.loc[:, "GFid_beta"] = 1
	df_gammaRec.loc[:, "GFid_Pcal1"] = 1
	df_gammaRec.loc[:, "GFid_Pcal2"] = 0
	df_gammaRec.loc[df_gammaRec.Gsector>7, "GFid_Pcal2"] = 1
	df_gammaRec.loc[:, "GFid_FT"] = 1
	#passGammaPCALFiducialCut
	df_gammaRec.loc[(df_gammaRec.GcalV1 <= g_min_v) & (df_gammaRec.Gsector<7), "GFid_Pcal1"] = 0
	df_gammaRec.loc[(df_gammaRec.GcalW1 <= g_min_w) & (df_gammaRec.Gsector<7), "GFid_Pcal1"] = 0
	#passGammaBetaCut
	df_gammaRec.loc[df_gammaRec.Gbeta <= min_Gbeta, "GFid_beta"] = 0
	df_gammaRec.loc[df_gammaRec.Gbeta >= max_Gbeta, "GFid_beta"] = 0

	#photon FD fiducial cuts by F.X. Girod
	#apply photon fiducial cuts
	sector_cond = [df_gammaRec.Gsector ==1, df_gammaRec.Gsector ==2, df_gammaRec.Gsector ==3, df_gammaRec.Gsector ==4, df_gammaRec.Gsector ==5, df_gammaRec.Gsector ==6]
	psplit = np.select(sector_cond, [87, 82, 85, 77, 78, 82])
	tleft = np.select(sector_cond, [58.7356, 62.8204, 62.2296, 53.7756, 58.2888, 54.5822])
	tright = np.select(sector_cond, [58.7477, 51.2589, 59.2357, 56.2415, 60.8219, 49.8914])
	sleft = np.select(sector_cond, [0.582053, 0.544976, 0.549788, 0.56899, 0.56414, 0.57343])
	sright = np.select(sector_cond, [-0.591876, -0.562926, -0.562246, -0.563726, -0.568902, -0.550729])
	rleft = np.select(sector_cond, [64.9348, 64.7541, 67.832, 55.9324, 55.9225, 60.0997])
	rright = np.select(sector_cond, [65.424, 54.6992, 63.6628, 57.8931, 56.5367, 56.4641])
	qleft = np.select(sector_cond, [0.745578, 0.606081, 0.729202, 0.627239, 0.503674, 0.717899])
	qright = np.select(sector_cond, [-0.775022, -0.633863, -0.678901, -0.612458, -0.455319, -0.692481])
	#first condition
	ang = np.radians((df_gammaRec.loc[df_gammaRec.Gsector<7, "Gsector"]-1) * 60)
	GcX_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.sin(ang) + df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.cos(ang)
	GcY_rot = df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] * np.cos(ang) - df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] * np.sin(ang)

	df_gammaRec.loc[df_gammaRec.Gsector<7, "GcX"] = GcX_rot
	df_gammaRec.loc[df_gammaRec.Gsector<7, "GcY"] = GcY_rot

	cond1_1 = df_gammaRec.GcX >= psplit
	cond1_2 = df_gammaRec.GcY < sleft * (df_gammaRec.GcX - tleft)
	cond1_3 = df_gammaRec.GcY > sright * (df_gammaRec.GcX - tright)
	cond1_4 = df_gammaRec.Gsector < 7
	cond1 = cond1_1 & cond1_2 & cond1_3 & cond1_4
	df_gammaRec.loc[cond1, "GFid_Pcal2"] = 1
	#second condition else if the first
	# cond2_0 = df_gammaRec.GFid == 0 # not necessary, because cond2_1 rules out the first (S. Lee)
	cond2_1 = df_gammaRec.GcX < psplit
	cond2_2 = df_gammaRec.GcY < qleft * (df_gammaRec.GcX - rleft)
	cond2_3 = df_gammaRec.GcY > qright * (df_gammaRec.GcX - rright)
	cond2_4 = df_gammaRec.Gsector < 7
	cond2 = cond2_1 & cond2_2 & cond2_3 & cond2_4
	df_gammaRec.loc[cond2, "GFid_Pcal2"] = 1

	#FT fiducial cuts
	circleCenterX1 = -8.419
	circleCenterY1 = 9.889
	circleRadius1 = 1.6

	circleCenterX2 = -9.89
	circleCenterY2 = -5.327
	circleRadius2 = 1.6

	circleCenterX3 = -6.15
	circleCenterY3 = -13
	circleRadius3 = 2.3

	circleCenterX4 = 3.7
	circleCenterY4 = -6.5
	circleRadius4 = 2

	circle1 = (df_gammaRec.GcX - circleCenterX1)**2 + (df_gammaRec.GcY - circleCenterY1)**2 < circleRadius1**2
	circle2 = (df_gammaRec.GcX - circleCenterX2)**2 + (df_gammaRec.GcY - circleCenterY2)**2 < circleRadius2**2
	circle3 = (df_gammaRec.GcX - circleCenterX3)**2 + (df_gammaRec.GcY - circleCenterY3)**2 < circleRadius3**2
	circle4 = (df_gammaRec.GcX - circleCenterX4)**2 + (df_gammaRec.GcY - circleCenterY4)**2 < circleRadius4**2

	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle1, "GFid_FT"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle2, "GFid_FT"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle3, "GFid_FT"] = 0
	df_gammaRec.loc[(df_gammaRec.Gsector > 7) & circle4, "GFid_FT"] = 0

	exclusion1_1 = (df_gammaRec.GcalW1 > 74) & (df_gammaRec.GcalW1 < 79.8)
	exclusion1_2 = (df_gammaRec.GcalW1 > 83.6) & (df_gammaRec.GcalW1 < 92.2)
	exclusion1_3 = (df_gammaRec.GcalW1 > 212.5) & (df_gammaRec.GcalW1 < 230)
	exclusion1 = exclusion1_1 | exclusion1_2 | exclusion1_3
	df_gammaRec.loc[(df_gammaRec.Gsector == 1) & exclusion1, "GFid_FT"] = 0
	exclusion2_1 = (df_gammaRec.GcalW1 < 14)
	exclusion2_2 = (df_gammaRec.GcalU1 > 111.2) & (df_gammaRec.GcalU1 < 119.3)
	exclusion2_3 = (df_gammaRec.GcalV1 > 113) & (df_gammaRec.GcalV1 < 118.7)
	exclusion2 = exclusion2_1 | exclusion2_2 | exclusion2_3
	df_gammaRec.loc[(df_gammaRec.Gsector == 2) & exclusion2, "GFid_FT"] = 0
	exclusion3 = df_gammaRec.GcalW1 < 14
	df_gammaRec.loc[(df_gammaRec.Gsector == 3) & exclusion3, "GFid_FT"] = 0
	exclusion4_1 = (df_gammaRec.GcalV1 < 14)
	exclusion4_2 = (df_gammaRec.GcalV1 > 229.4) & (df_gammaRec.GcalV1 < 240.7)
	exclusion4_3 = (df_gammaRec.GcalW1 > 135) & (df_gammaRec.GcalW1 < 150)
	exclusion4 = exclusion4_1 | exclusion4_2 | exclusion4_3
	df_gammaRec.loc[(df_gammaRec.Gsector == 4) & exclusion4, "GFid_FT"] = 0
	exclusion6 = (df_gammaRec.GcalW1 > 170) & (df_gammaRec.GcalW1 < 192)
	df_gammaRec.loc[(df_gammaRec.Gsector == 6) & exclusion6, "GFid_FT"] = 0

	df_gammaRec.loc[:, "GFid_Pcal"] = 	df_gammaRec.loc[:, "GFid_Pcal1"] * 	df_gammaRec.loc[:, "GFid_Pcal2"]
	df_gammaRec.loc[:, "GFid"] = df_gammaRec.GFid_beta * df_gammaRec.GFid_Pcal1 * df_gammaRec.GFid_Pcal2 * df_gammaRec.GFid_FT 
	return df_gammaRec

def protonFiducialCounting(df_protonRec, pol = 'inbending', fidlevel = 'mid'):
	df_protonRec.loc[:, "PFid_dc"] = 1
	df_protonRec.loc[:, "PFid_cvt"] = 0
	df_protonRec.loc[df_protonRec.Psector<7, "PFid_cvt"] = 1 #FD fid done by previous pipeline
	df_protonRec.loc[:, "PFid_chi"] = 1
	df_protonRec.loc[:, "PFid_vz"] = 1

	dcsec = determineSector(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity)
	if pol == 'inbending':
		minparams = p_dc_minparams_in
		maxparams = p_dc_maxparams_in

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity, df_protonRec.PDc1Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 0, minparams, maxparams)
		df_protonRec.loc[(phi_DC<=phi_DC_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(phi_DC>=phi_DC_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity, df_protonRec.PDc2Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 1, minparams, maxparams)
		df_protonRec.loc[(phi_DC<=phi_DC_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(phi_DC>=phi_DC_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0

		theta_DC, phi_DC = thetaphifromhit(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity, df_protonRec.PDc3Hitz)
		phi_DC_min, phi_DC_max = p_DC_fiducial_cut_thetaphi(theta_DC, dcsec, 2, minparams, maxparams)
		df_protonRec.loc[(phi_DC<=phi_DC_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(phi_DC>=phi_DC_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0

		pchi2CD_lb,   pchi2CD_ub   = -5.000, 6.345
		pchi2FD_S1_lb, pchi2FD_S1_ub = -3.296, 3.508
		pchi2FD_S2_lb, pchi2FD_S2_ub = -3.552, 4.000
		pchi2FD_S3_lb, pchi2FD_S3_ub = -3.446, 3.937
		pchi2FD_S4_lb, pchi2FD_S4_ub = -2.747, 3.190
		pchi2FD_S5_lb, pchi2FD_S5_ub = -2.851, 3.418
		pchi2FD_S6_lb, pchi2FD_S6_ub = -3.174, 3.514

		vzdiffCD_lb,    vzdiffCD_ub    = -2.011, 2.314
		vzdiffFD_S1_lb, vzdiffFD_S1_ub = -3.209, 4.017
		vzdiffFD_S2_lb, vzdiffFD_S2_ub = -3.612, 4.139
		vzdiffFD_S3_lb, vzdiffFD_S3_ub = -3.328, 4.287
		vzdiffFD_S4_lb, vzdiffFD_S4_ub = -3.411, 4.108
		vzdiffFD_S5_lb, vzdiffFD_S5_ub = -3.607, 4.246
		vzdiffFD_S6_lb, vzdiffFD_S6_ub = -2.999, 3.927

	if pol == 'outbending':
		minparams = p_dc_minparams_out
		maxparams = p_dc_maxparams_out

		dcsec = determineSector(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc1Hitx, df_protonRec.PDc1Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 0, minparams, maxparams)
		df_protonRec.loc[(y_rot<=calc_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(y_rot>=calc_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		#passElectronDCR2
		dcsec = determineSector(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc2Hitx, df_protonRec.PDc2Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 1, minparams, maxparams)
		df_protonRec.loc[(y_rot<=calc_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(y_rot>=calc_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0

		#passElectronDCR3
		dcsec = determineSector(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity)
		x_rot, y_rot = rotateDCHitPosition(df_protonRec.PDc3Hitx, df_protonRec.PDc3Hity, dcsec)
		calc_min, calc_max = p_DC_fiducial_cut_XY(x_rot, dcsec, 2, minparams, maxparams)
		df_protonRec.loc[(y_rot<=calc_min) & (df_protonRec.Psector<7), "PFid_dc"] = 0
		df_protonRec.loc[(y_rot>=calc_max) & (df_protonRec.Psector<7), "PFid_dc"] = 0

		pchi2CD_lb,   pchi2CD_ub   = -5.592,  6.785
		pchi2FD_S1_lb, pchi2FD_S1_ub = -3.905, 4.088
		pchi2FD_S2_lb, pchi2FD_S2_ub = -3.411, 3.939
		pchi2FD_S3_lb, pchi2FD_S3_ub = -4.042, 5.954
		pchi2FD_S4_lb, pchi2FD_S4_ub = -3.820, 5.065
		pchi2FD_S5_lb, pchi2FD_S5_ub = -3.384, 4.232
		pchi2FD_S6_lb, pchi2FD_S6_ub = -5.077, 5.100


	cut_CD = df_protonRec.Psector > 7
	if fidlevel == 'mid':
		cut_right = cut_CD & (df_protonRec.Ptheta<max_Ptheta)
	elif fidlevel == 'tight':
		cut_right = cut_CD & (df_protonRec.Ptheta<max_Ptheta-5)
	cut_bottom = cut_CD & (df_protonRec.PCvt12theta>44.5)
	cut_sidel = cut_CD & (df_protonRec.PCvt12theta<-2.942 + 1.274*df_protonRec.Ptheta)
	cut_sider = cut_CD & (df_protonRec.PCvt12theta>-3.523 + 1.046*df_protonRec.Ptheta)

	cut_trapezoid = cut_CD & cut_right & cut_bottom & cut_sidel & cut_sider

	cut_gaps1 = ~((df_protonRec.PCvt12phi>-95) & (df_protonRec.PCvt12phi<-80))
	cut_gaps2 = ~((df_protonRec.PCvt12phi>25) & (df_protonRec.PCvt12phi<40))
	cut_gaps3 = ~((df_protonRec.PCvt12phi>143) & (df_protonRec.PCvt12phi<158))
	cut_gaps = cut_CD & cut_gaps1 & cut_gaps2 & cut_gaps3
	cut_total = cut_gaps & cut_trapezoid

	df_protonRec.loc[cut_total, "PFid_cvt"] = 1 #CD fid
            

	df_protonRec.loc[ (df_protonRec.Psector>4000) & ((df_protonRec.Pchi2pid<pchi2CD_lb)   | (df_protonRec.Pchi2pid>pchi2CD_ub)  ), "PFid_chi"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==1)   & ((df_protonRec.Pchi2pid<pchi2FD_S1_lb) | (df_protonRec.Pchi2pid>pchi2FD_S1_ub)), "PFid_chi"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==2)   & ((df_protonRec.Pchi2pid<pchi2FD_S2_lb) | (df_protonRec.Pchi2pid>pchi2FD_S2_ub)), "PFid_chi"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==3)   & ((df_protonRec.Pchi2pid<pchi2FD_S3_lb) | (df_protonRec.Pchi2pid>pchi2FD_S3_ub)), "PFid_chi"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==4)   & ((df_protonRec.Pchi2pid<pchi2FD_S4_lb) | (df_protonRec.Pchi2pid>pchi2FD_S4_ub)), "PFid_chi"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==5)   & ((df_protonRec.Pchi2pid<pchi2FD_S5_lb) | (df_protonRec.Pchi2pid>pchi2FD_S5_ub)), "PFid_chi"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==6)   & ((df_protonRec.Pchi2pid<pchi2FD_S6_lb) | (df_protonRec.Pchi2pid>pchi2FD_S6_ub)), "PFid_chi"] = 0

	vzdiffCD_lb,    vzdiffCD_ub    = -2.737, 2.096
	vzdiffFD_S1_lb, vzdiffFD_S1_ub = -4.435, 3.429
	vzdiffFD_S2_lb, vzdiffFD_S2_ub = -4.646, 2.978
	vzdiffFD_S3_lb, vzdiffFD_S3_ub = -3.922, 3.040
	vzdiffFD_S4_lb, vzdiffFD_S4_ub = -4.646, 3.493
	vzdiffFD_S5_lb, vzdiffFD_S5_ub = -3.901, 3.750
	vzdiffFD_S6_lb, vzdiffFD_S6_ub = -3.846, 3.623

	df_protonRec.loc[ (df_protonRec.Psector>4000) & ((df_protonRec.vzdiff<vzdiffCD_lb)   | (df_protonRec.vzdiff>vzdiffCD_ub)  ), "PFid_vz"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==1)   & ((df_protonRec.vzdiff<vzdiffFD_S1_lb) | (df_protonRec.vzdiff>vzdiffFD_S1_ub)), "PFid_vz"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==2)   & ((df_protonRec.vzdiff<vzdiffFD_S2_lb) | (df_protonRec.vzdiff>vzdiffFD_S2_ub)), "PFid_vz"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==3)   & ((df_protonRec.vzdiff<vzdiffFD_S3_lb) | (df_protonRec.vzdiff>vzdiffFD_S3_ub)), "PFid_vz"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==4)   & ((df_protonRec.vzdiff<vzdiffFD_S4_lb) | (df_protonRec.vzdiff>vzdiffFD_S4_ub)), "PFid_vz"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==5)   & ((df_protonRec.vzdiff<vzdiffFD_S5_lb) | (df_protonRec.vzdiff>vzdiffFD_S5_ub)), "PFid_vz"] = 0
	df_protonRec.loc[ (df_protonRec.Psector==6)   & ((df_protonRec.vzdiff<vzdiffFD_S6_lb) | (df_protonRec.vzdiff>vzdiffFD_S6_ub)), "PFid_vz"] = 0

	df_protonRec.loc[:, "PFid"] = df_protonRec.PFid_dc * df_protonRec.PFid_cvt * df_protonRec.PFid_chi * df_protonRec.PFid_vz 
	return df_protonRec
