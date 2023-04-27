# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 00:40:01 2023

@author: 91956
"""

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import os
import math
import time as clock

########################################
########### basic constants ############
########################################

# os.chdir('/home/antony/mycodes/astron_euler/data')
start_time = clock.time()
dt = 0.00005 #0.00005seconds = 0.05ms = 50us = 50000ns
stim_time = 300 #seconds
N = int(stim_time / dt)

#some constants
# setting up cel constants
volume_ratio = 2
cyt_vol = 0.37e-18
cyt_area = 4.921e-12
pmca_density = 400e12
ava_number = 6.022e23

# release machinery constants
gama_1 = 0.000417
gama_2 = 115
gama_3 = 20#20 in anup's, 8 in kiran's


########################################
############## functions ###############
########################################


# pmca flux
def PMCA_flux(Ca_cyt, m0, m1, m2):
	kf1 = 15
	kf2 = 20
	kf3 = 100
	k1 = 0.6
	kb1 = 20
	#print(Ca_cyt,m0,m1,m2)
	dm0dt = ((-m0 * Ca_cyt * kf1) + (m1 * kb1) + (m2 * kf3)) * dt #+ np.sqrt(dt)*np.random.normal()*m0
	dm1dt = ((-m1 * kb1) + (m0 * Ca_cyt * kf1) + (-m1 * kf2)) * dt #+ np.sqrt(dt)*np.random.normal()*m1
	dm2dt = ((-m2 * kf3) + (m1 * kf2)) * dt# + np.sqrt(dt)*np.random.normal()*m2
	J_pmca = ((-m0 * Ca_cyt * kf1) + (m1 * kb1) + (m0 * k1)) * dt
	#print(J_pmca,dm0dt,dm1dt,dm2dt)
	return J_pmca, dm0dt, dm1dt, dm2dt

# cytosolic buffer
def CB_flux(Ca_cyt, c0, c1):
	kf = 60
	kb = 1200
	dc0 = ((-c0 * Ca_cyt * kf) + (c1 * kb)) * dt
	dc1 = ((-c1 * kb) + (c0 * Ca_cyt * kf)) * dt
	#print(c0+c1)
	J_cytB = ((-Ca_cyt * c0 * kf) + (c1 * kb)) * dt
	return J_cytB, dc0, dc1

# serca flux (hill)
def SERCA_flux(Ca_cyt):
	V_max_serca = 250
	ka = 0.1
	n = 2
	J_serca = (V_max_serca * ((Ca_cyt ** n) / ((Ca_cyt ** n) + (ka ** n)))) * dt
	return J_serca

# er leak flux
def ER_leak_flux(Ca_cyt, Ca_ER):
	k_leak = 0.2
	J_ERleak = (k_leak * (Ca_ER - Ca_cyt)) * dt
	return J_ERleak

# ER buffer flux
def ER_buffer_flux(Ca_ER, B0, B1):
	kfl = 1
	kbl = 80
	# print(B0+B1)
	dB0 = ((-B0 * Ca_ER * kfl) + (B1 * kbl)) * dt
	dB1 = ((-B1 * kbl) + (B0 * Ca_ER * kfl)) * dt
	J_ERB = ((-Ca_ER * B0 * kfl) + (B1 * kbl)) * dt
	return J_ERB, dB0, dB1

# IP3 flux
def IP3_flux(Ca_cyt, Ca_ER, IP3_conc, h):
	V_max_ip3R = 10
	a1 = 400
	a2 = 0.2
	a3 = 400
	a5 = 20
	b1 = 52
	b2 = 0.2
	b3 = 377.36
	b5 = 1.64
	d1 = b1 / a1
	d2 = b2 / a2
	d3 = b3 / a3
	d5 = b5 / a5
	q2 = d2 * ((IP3_conc + d1) / (IP3_conc + d3))
	n = 5
	n_inf = IP3_conc / (IP3_conc + d1)
	m_inf = Ca_cyt / (Ca_cyt + d5)
	alpha_h = (a2 * d2) * ((IP3_conc + d1) / (IP3_conc + d3))
	beta_h = a2 * Ca_cyt
	varience = (alpha_h * (1 - h) + beta_h * h) / n
	# print(Ca_cyt, IP3_conc, varience)
	std_dev = math.sqrt(2*abs(varience)*dt)
	J_IP3R = (V_max_ip3R * (m_inf ** 3) * (n_inf ** 3) * (h ** 3) * (Ca_ER - Ca_cyt)) * dt
	dh = ((alpha_h * (1 - h)) - (beta_h * h)) * dt + np.random.normal()*std_dev #noise: np.random.normal()*std_dev
	return J_IP3R, dh

# IP3 kinase
def IP3_kinase(Ca_cyt, IP3_conc):
	V_max_ip3K = 5
	ka1 = 0.4
	ka2 = 10
	J_IP3K = (V_max_ip3K * ((Ca_cyt ** 4) / ((Ca_cyt ** 4) + (ka1 ** 4))) * (IP3_conc / (IP3_conc + ka2))) * dt
	return J_IP3K

#IP3 phospatase
def IP3_phos(IP3_conc):
	V_max_ip3P = 1.25
	IP3_base = 0.160
	J_IP3P = (V_max_ip3P * (IP3_conc - IP3_base)) * dt
	return J_IP3P

# PLC delta pathway
def IP3_plc_delta(Ca_cyt, IP3_conc):
	V_max_ip3PLC = 0.02
	ka1_plc = 1
	ka2_plc = 1.5
	J_IP3PLC = ((V_max_ip3PLC / (1 + (IP3_conc / ka1_plc))) * ((Ca_cyt ** 2) / ((Ca_cyt ** 2) + (ka2_plc**2)))) * dt
	return J_IP3PLC

# glutamate degradation and flux
def Glutamate_flux(Glu_conc):
	kGlu = 160
	dGlu = (-Glu_conc * kGlu) * dt
	V_max_mGluR = 0.65
	ka_glu = 11
	J_mGluR = (V_max_mGluR * ((Glu_conc ** 2) / ((Glu_conc ** 2) + (ka_glu ** 2)))) * dt
	return J_mGluR, dGlu

#syt4 flux
#added cooperativity term
def syt4_flux(Ca_cyto, s0, s1, s2):
	kf4 = 153 #1.53e8 in anup's code, 153 in kiran's
	kb4 = 3500 #3500 in anup's code
	b4 = 0.25 #0.25 in anup's code
	ds0 = ((s1 * kb4) - (s0 * 2 * kf4 * Ca_cyto)) * dt
	ds1 = ((s0 * 2 * kf4 * Ca_cyto) + (s2 * 2 * kb4 * b4) - (s1 * kb4) - (s1 * kf4 * Ca_cyto)) * dt
	ds2 = ((s1 * 1 * kf4 * Ca_cyto) - (s2 * 2 * kb4 * b4)) * dt
	J_syt4 = ((-s0 * 2 * kf4 * Ca_cyto) + (s1 * kb4) - (s1 * kf4 * Ca_cyto) + (s2 * 2 * kb4 * b4)) * dt
	return J_syt4, ds0, ds1, ds2

#syt7 flux
#added cooperativity term
def syt7_flux(Ca_cyto, y0, y1, y2, y3, y4, y5):
	kf7 = 3.82 #3.82e6 in anup's code, 3.82 in kiran's
	kb7 = 60 #60 in anup's code
	b7 = 0.25 #0.25 in anup's code
	dy0 = ((y1 * kb7) - (y0 * 5 * kf7 * Ca_cyto)) * dt
	dy1 = ((y0 * 5 * kf7 * Ca_cyto) + (y2 * 2 * kb7 * (b7 ** 1)) - (y1 * 1 * kb7 * (b7 ** 0)) - (y1 * 4 * kf7 * Ca_cyto)) * dt
	dy2 = ((y1 * 4 * kf7 * Ca_cyto) + (y3 * 3 * kb7 * (b7 ** 2)) - (y2 * 2 * kb7 * (b7 ** 1)) - (y2 * 3 * kf7 * Ca_cyto)) * dt
	dy3 = ((y2 * 3 * kf7 * Ca_cyto) + (y4 * 4 * kb7 * (b7 ** 3)) - (y3 * 3 * kb7 * (b7 ** 2)) - (y3 * 2 * kf7 * Ca_cyto)) * dt
	dy4 = ((y3 * 2 * kf7 * Ca_cyto) + (y5 * 5 * kb7 * (b7 ** 4)) - (y4 * 4 * kb7 * (b7 ** 3)) - (y4 * 1 * kf7 * Ca_cyto)) * dt
	dy5 = ((y4 * 1 * kf7 * Ca_cyto) - (y5 * 5 * kb7 * (b7 ** 4))) * dt
	#print(y0,y1,y2,y3,y4,y5)
	J_syt7 = ((-5 * y0 * kf7 * Ca_cyto) + (1 * y1 * kb7 * (b7 ** 0)) + (
		-4 * y1 * kf7 * Ca_cyto) + (2 * y2 * kb7 * (b7 ** 1)) + (
		-3 * y2 * kf7 * Ca_cyto) + (3 * y3 * kb7 * (b7 ** 2)) + (
		-2 * y3 * kf7 * Ca_cyto) + (4 * y4 * kb7 * (b7 ** 3)) + (
		-1 * y4 * kf7 * Ca_cyto) + (5 * y5 * kb7 * (b7 ** 4))) * dt
	return J_syt7, dy0, dy1, dy2, dy3, dy4, dy5

### recycling scheme
def recycling_scheme(docked, mobile, doc_released, doc_reacidified, doc_endocytosed, mob_released, mob_reacidified, mob_endocytosed, flag, mem):
	r04 = 0.000417
	rCa4 = 0#115 in anup's code, 115 in kiran's
	r07 = 0.000417
	rCa7 = 0#20 in anup's code, 8 in kiran's
	kef = 0.66
	kaf = 0.6
	kes = 0.16
	kas = 0.052
	k_mob1 = 0.615 #0.615 #make k_mob1=0 for suppressing anything ever going to docked
	k_mob2 = 0.615
	k_doc = 0.75 #1 in anup's code (i think) (rmob2doc=1, but fracmob2doc=0.25), 0.75 in kiran's #make k_doc=0 for suppressing anything ever going to docked. otherwise k_doc=0.75
	#k_doc is 0.75 in kiran's code because in anup's, th term is "(mobile * rmob2doc * (1-fracmob2doc))", which turns out to be 0.75
	if flag == 0:
		rCa4 = 0
		rCa7 = 0
	elif flag == 1:
		rCa4 = 115
		rCa7 = 0
	elif flag == 2:
		rCa4 = 115
		rCa7 = 0
	elif flag == 3:
		rCa4 = 0
		rCa7 = 20#8
	if mem == 0: #no release or normal recycle scheme
		ddocked = ((k_doc * mobile) - ((rCa4 + r04) * docked)) * dt
		dmobile = ((k_mob1 * doc_reacidified) + (k_mob2 * mob_reacidified) - (k_doc * mobile) - ((rCa7 + r07) * mobile)) * dt
		ddoc_released = (((rCa4 + r04) * docked) - (kef * doc_released)) * dt
		ddoc_endocytosed = ((kef * doc_released) - (kaf * doc_endocytosed)) * dt
		ddoc_reacidified = ((kaf * doc_endocytosed) - (k_mob1 * doc_reacidified)) * dt
		dmob_released = (((rCa7 + r07) * mobile) - (kes * mob_released)) * dt
		dmob_endocytosed = ((kes * mob_released) - (kas * mob_endocytosed)) * dt
		dmob_reacidified = ((kas * mob_endocytosed) - (k_mob2 * mob_reacidified)) * dt
	elif mem == 1: #mobile pool < 1
		ddocked = (- ((rCa4 + r04) * docked)) * dt
		dmobile = ((k_mob1 * doc_reacidified) + (k_mob2 * mob_reacidified)) * dt
		ddoc_released = (((rCa4 + r04) * docked) - (kef * doc_released)) * dt
		ddoc_endocytosed = ((kef * doc_released) - (kaf * doc_endocytosed)) * dt
		ddoc_reacidified = ((kaf * doc_endocytosed) - (k_mob1 * doc_reacidified)) * dt
		dmob_released = ( - (kes * mob_released)) * dt
		dmob_endocytosed = ((kes * mob_released) - (kas * mob_endocytosed)) * dt
		dmob_reacidified = ((kas * mob_endocytosed) - (k_mob2 * mob_reacidified)) * dt
	return ddocked, dmobile, ddoc_released, ddoc_reacidified, ddoc_endocytosed, dmob_released, dmob_reacidified, dmob_endocytosed


########################################
########## integration loop ############
########################################


for j in range(0,5):
	# initialising the parameters
	time = np.zeros(N, float)
	J_PMCA = np.zeros(N, float)
	m0 = np.zeros(N, float)
	m1 = np.zeros(N, float)
	m2 = np.zeros(N, float)
	J_cytB = np.zeros(N, float)
	c0 = np.zeros(N, float)
	c1 = np.zeros(N, float)
	J_serca = np.zeros(N, float)
	J_ERleak = np.zeros(N, float)
	J_ERB = np.zeros(N, float)
	b0 = np.zeros(N, float)
	b1 = np.zeros(N, float)
	J_IP3Flux = np.zeros(N, float)
	h = np.zeros(N, float)
	J_IP3Kinase = np.zeros(N, float)
	J_IP3Phos = np.zeros(N, float)
	J_IP3PLC = np.zeros(N, float)
	J_mGluR = np.zeros(N, float)
	Glu_conc = np.zeros(N, float)
	J_syt4 = np.zeros(N, float)
	s0 = np.zeros(N, float)
	s1 = np.zeros(N, float)
	s2 = np.zeros(N, float)
	J_syt7 = np.zeros(N, float)
	y0 = np.zeros(N, float)
	y1 = np.zeros(N, float)
	y2 = np.zeros(N, float)
	y3 = np.zeros(N, float)
	y4 = np.zeros(N, float)
	y5 = np.zeros(N, float)

	Ca_cyto = np.zeros(N, float)
	Ca_ER = np.zeros(N, float)
	IP3_conc = np.zeros(N, float)

	docked = np.zeros(N, float)
	mobile = np.zeros(N, float)
	doc_released = np.zeros(N, float)
	doc_reacidified = np.zeros(N, float)
	doc_endocytosed = np.zeros(N, float)
	mob_released = np.zeros(N, float)
	mob_reacidified = np.zeros(N, float)
	mob_endocytosed = np.zeros(N, float)

	release = np.zeros(N, int)
	spontaneous = np.zeros(N, float)
	synchronous = np.zeros(N, float)
	asynchronous = np.zeros(N, float)

	#setting initial conditions
	# initial conditons for PMCA
	m1[0] = 0.252
	m2[0] = 0.050
	m0[0] = 1.0 * (pmca_density * cyt_area) / (ava_number * cyt_vol * 1000) * 1e6 - (m1[0]+m2[0])
	# initial conditions for cytosolic buffer
	c0[0] = 49.804
	c1[0] = 0.196
	# initial conditions for ER buffer
	b0[0] = 1.75
	b1[0] = 8.25
	# initial conditions for recycling scheme
	docked[0] = 0.800 #0.800
	mobile[0] = 0.250 #make docked=0 and mobile=0.8 to study only asynchronous release. otherwise 0.2. anup said 0.25 for mobile.
	# initial conditions for syt4
	s0[0] = 0.993
	s1[0] = 0.007
	s2[0] = 0
	# initial conditions for syt7
	y0[0] = 0.975
	y1[0] = 0.025
	y2[0] = 0
	y3[0] = 0
	y4[0] = 0
	y5[0] = 0
	#initial conditions for ion concentrations
	IP3_conc[0] = 0.160
	Ca_cyto[0] = 0.0
	Ca_ER[0] = 400
	
	#integration step
	mem = 0 #memory of vesicle release in case of asynchronous release
	for i in range(1, N):
		release_flag = 0
		#increase in time
		time[i] = time[i - 1] + dt
		#pmca flux
		pmca_temp = PMCA_flux(Ca_cyto[i - 1], m0[i - 1], m1[i - 1], m2[i - 1])
		J_PMCA[i] = pmca_temp[0]
		m0[i] = m0[i - 1] + pmca_temp[1]
		m1[i] = m1[i - 1] + pmca_temp[2]
		m2[i] = m2[i - 1] + pmca_temp[3]
		#cytosolic buffer
		cytB_temp = CB_flux(Ca_cyto[i - 1], c0[i - 1], c1[i - 1])
		J_cytB[i] = cytB_temp[0]
		c0[i] = c0[i - 1] + cytB_temp[1]
		c1[i] = c1[i - 1] + cytB_temp[2]
		#serca flux
		J_serca[i] = SERCA_flux(Ca_cyto[i - 1])
		#er leak
		J_ERleak[i] = ER_leak_flux(Ca_cyto[i - 1], Ca_ER[i - 1])
		#er buffer
		ER_buffer_temp = ER_buffer_flux(Ca_ER[i - 1], b0[i - 1], b1[i - 1])
		J_ERB[i] = ER_buffer_temp[0]
		b0[i] = b0[i - 1] + ER_buffer_temp[1]
		b1[i] = b1[i - 1] + ER_buffer_temp[2]
		#ip3 flux
		while True:
			IP3_flux_temp = IP3_flux(Ca_cyto[i - 1], Ca_ER[i - 1], IP3_conc[i - 1], h[i - 1])
			J_IP3Flux[i] = IP3_flux_temp[0]
			h[i] = h[i - 1] + IP3_flux_temp[1]
			if 0 <= h[i] <= 1: #or True:
				break
		#ip3 kinase
		J_IP3Kinase[i] = IP3_kinase(Ca_cyto[i - 1], IP3_conc[i - 1])
		#ip3 phosphatase
		J_IP3Phos[i] = IP3_phos(IP3_conc[i - 1])
		#plc delta pathway
		J_IP3PLC[i] = IP3_plc_delta(Ca_cyto[i - 1], IP3_conc[i - 1])
		#glutamate flux and degradation
		Glutamate_temp = Glutamate_flux(Glu_conc[i - 1])
		J_mGluR[i] = Glutamate_temp[0]
		# if 240 < dt * i < 270 and np.sin(dt*i*2*np.pi*0.5) >0.8:
		if 210 < dt * i < 240:
			Glu_conc[i] = 400
		else:
			Glu_conc[i] = Glu_conc[i - 1] + Glutamate_temp[1]
		syt4_temp = syt4_flux(Ca_cyto[i-1], s0[i-1], s1[i-1], s2[i-1])
		#syt4 flux
		J_syt4[i] = syt4_temp[0]
		s0[i] = (s0[i-1] + syt4_temp[1])
		s1[i] = (s1[i - 1] + syt4_temp[2])
		s2[i] = (s2[i - 1] + syt4_temp[3])
		#syt7 flux
		syt7_temp = syt7_flux(Ca_cyto[i-1], y0[i-1], y1[i-1], y2[i-1], y3[i-1], y4[i-1], y5[i-1])
		J_syt7[i] = syt7_temp[0]
		y0[i] = (y0[i-1] + syt7_temp[1])
		y1[i] = (y1[i - 1] + syt7_temp[2])
		y2[i] = (y2[i - 1] + syt7_temp[3])
		y3[i] = (y3[i - 1] + syt7_temp[4])
		y4[i] = (y4[i - 1] + syt7_temp[5])
		y5[i] = (y5[i - 1] + syt7_temp[6])
		#ion concentrations
		Ca_cyto[i] = Ca_cyto[i - 1] + (J_PMCA[i] + J_cytB[i] - J_serca[i] + J_ERleak[i] + J_IP3Flux[i] + J_syt4[i] + J_syt7[i])
		Ca_ER[i] = Ca_ER[i - 1] + ((J_serca[i] - J_ERleak[i] - J_IP3Flux[i]) * volume_ratio + J_ERB[i])
		IP3_conc[i] = IP3_conc[i - 1] + (J_mGluR[i] - J_IP3Kinase[i] - J_IP3Phos[i] + J_IP3PLC[i])
		
#	 #calculating release probabilities - kiran's scheme
#	 spontaneous_r = gama_1 * s0[i]*y0[i]  *dt
#	 spontaneous[i] = spontaneous_r
#	 synchronous_r = gama_2 * (s2[i]*y0[i] + s2[i]*y1[i] + s2[i]*y2[i] + s2[i]*y3[i] + s2[i]*y4[i] + s2[i]*y5[i])*dt
#	 synchronous[i] = synchronous_r
#	 asynchronous_r = gama_3 * (s0[i]*y5[i] + s1[i]*y5[i] + s2[i]*y5[i])*dt
#	 asynchronous[i] = asynchronous_r
#	 #suhita's suggestion of release probabilities
#	 spontaneous_r = gama_1 * s0[i]  *dt
#	 spontaneous[i] = spontaneous_r
#	 synchronous_r = gama_2 * s2[i] *dt
#	 synchronous[i] = synchronous_r
#	 asynchronous_r = gama_3 * y5[i] *dt
#	 asynchronous[i] = asynchronous_r
	#different scheme for release rates, as taken from anup's code, but in anup's code he calculates vesicle recycling scheme first and then the rates and i calculate rates first so im not sure ho this will work out
	'''	#synchronous release
		s0frac = s0[i-1]/(s0[i-1]+s1[i-1]+s2[i-1])
		s2frac = s2[i-1]/(s0[i-1]+s1[i-1]+s2[i-1])
		synchronous_r = gama_2 * docked[i-1] * s2frac * dt
		#asynchronous release
		y0frac = y0[i-1]/(y0[i-1]+y1[i-1]+y2[i-1]+y3[i-1]+y4[i-1]+y5[i-1])
		y5frac = y5[i-1]/(y0[i-1]+y1[i-1]+y2[i-1]+y3[i-1]+y4[i-1]+y5[i-1])
		asynchronous_r = gama_3 * mobile[i-1] * y5frac * dt
		#spontaneous release
		spontaneous_r = gama_1 * ((docked[i-1] * s0frac) + (mobile[i-1] * y0frac)) * dt
		
		if spontaneous_r > np.random.uniform():
			release_flag = 1
			release[i] = 1
			mem = 1 #mem, which was 0 before starting the integration loop, has now become 1, i.e., synchronous release has happened
			print("Spontaneous release at: ", time[i])
		if ((asynchronous_r > np.random.uniform()) and (mobile[i-1] > 0.9)):
		#if asynchronous_r > np.random.uniform():
			release_flag = 3
			release[i] = 3
			mem = 1 #mem, which was 0 before starting the integration loop, has now become 1, i.e., asynchronous release has happened
			print("Asynchronous release at: ", time[i])
		if ((synchronous_r > np.random.uniform()) and (docked[i-1] > 0.99)):
			release_flag = 2
			release[i] = 2
			mem = 1 #mem, which was 0 before starting the integration loop, has now become 1, i.e., synchronous release has happened
			print("Synchronous release at: ", time[i]
    '''
	# this lower else block had me confused for a very long time. whatever you do, do not uncomment
	# #	 else:
	# #		 release[i] = 0

	# #uncomment the following and remove all the above recycle schemes to get back a singular vesicle release followed by overwriting
	'''	#vesicle recycling scheme
		recycling_temp = recycling_scheme(docked[i-1], mobile[i-1], doc_released[i - 1], doc_reacidified[i - 1], doc_endocytosed[i - 1],mob_released[i - 1], mob_reacidified[i - 1], mob_endocytosed[i - 1], release_flag, mem)
		mobile[i] = mobile[i-1] + recycling_temp[1]
		doc_released[i] = doc_released[i - 1] + recycling_temp[2]
		doc_reacidified[i] = doc_reacidified[i - 1] + recycling_temp[3]
		doc_endocytosed[i] = doc_endocytosed[i - 1] + recycling_temp[4]
		mob_released[i] = mob_released[i - 1] + recycling_temp[5]
		mob_reacidified[i] = mob_reacidified[i - 1] + recycling_temp[6]
		mob_endocytosed[i] = mob_endocytosed[i - 1] + recycling_temp[7]
		docked[i] = docked[i-1] + recycling_temp[0]

		#setting maximum limit = 1 or all 8 variables because these are fractions, not whole numbers
		if docked[i] > 1.5:
			docked[i] = 1.05
		if mobile[i] > 1.05:
			mobile[i] = 1.05
		if doc_released[i] > 1:
			doc_released[i] = 1
		if doc_endocytosed[i] > 1:
			doc_endocytosed[i] = 1
		if doc_reacidified[i] > 1:
			doc_reacidified[i] = 1
		if mob_released[i] > 1:
			mob_released[i] = 1
		if mob_endocytosed[i] > 1:
			mob_endocytosed[i] = 1
		if mob_reacidified[i] > 1:
			mob_reacidified[i] = 1

		#setting the vesicle numbers as 1 or 0 based on type of vesicle release
		if release[i] == 3:  #asynchronous release
			mobile[i] = 0
			mob_released[i] = 1
		elif release[i] == 2:  #synchronous release
			docked[i] = 0
			doc_released[i] = 1
		elif release[i] == 1:  #spontaneous release
			docked[i] = 0
			doc_released[i] = 1

		#updating memory of vesicle release i mobile pool = 1
		if mobile[i] > 0.99:
			mem = 0 # mem, which had been set to 1 upon asynchronous release, has been updated to 0 so normal recycle scheme to be used
    '''

	print("--- %s seconds ---" % (clock.time() - start_time)) #how long it took to run simulation
	#df = pd.DataFrame({"Time": time[4000000:], "Ca_cyto": Ca_cyto[4000000:],"Docked": docked[4000000:], "Released": doc_released[4000000:], "Endocytosed": doc_endocytosed[4000000:], "Reacidified": doc_reacidified[4000000:], "R_type": release[4000000:]}) #make a data frame but for values after 200s (time step 4000000)
	df = pd.DataFrame({"Time": time[0:], "Ca_cyto": Ca_cyto[0:],"Ca_ER": Ca_ER[0:], "IP3_conc": IP3_conc[0:]}) #make a data frame but for values after 200s (time step 4000000)
	df.to_csv(r"C:\semester project\shweta\sims\glu_400_30" + "\output" + str(j) + ".csv", index=False) #save data frame as csv
