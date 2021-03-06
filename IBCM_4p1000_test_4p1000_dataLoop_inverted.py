############################################################
#ICBM 
#Equations from (Andren and Katterer 1997) 
#(https://doi.org/10.1890/1051-0761(1997)007[1226:ITICBM]2.0.CO;2)

import sys
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.optimize import least_squares
import numdifftools as ndt
import math
from random import gauss
import xlrd
import pandas as pd
import time
import datetime
from datetime import datetime

#################
#Type of simulation
##################
optimization_4p1000=True
optimization_param=False
#######################

#Save data
if(optimization_param):
        out1=open("SOC_IBCM_optim6.txt","wb")
if(optimization_4p1000):
        out2=open("SOC_IBCM_4p1000.txt","wb")

##################
#Functions
##################
#########################################################
def AB_NanZeroRemover(site_T0,site_T0_name,iout,ROOTDIR):
#########################################################
        if ( iout > 0 ):
                out1=open(ROOTDIR+"AB.data","wb")
        #######################################################################
        # remove NaNs and ZEROs from ABOVE and BELOW
        #######################################################################
        yy=np.asarray(site_T0['Year']).astype(np.int16)             # array of years
        aa=np.asarray(site_T0['ABOVE']).astype(np.float64)/10.  # array of C_above (kgC/m2/yr)
        bb=np.asarray(site_T0['BELOW']).astype(np.float64)/10.  # array of C_below (kgC/m2/yr)
        aa0=np.where(np.isnan(aa),0,aa)  # replace NaN with zero in "above"
        YEAR=yy[aa0>0]                     # select years where above>0
        abo=aa[aa0>0]                      # select ABOVE>0
        bel=bb[aa0>0]                      # select corresponding BELOW
        if (iout > 0):
                XX=np.stack((YEAR,abo,bel),axis=0)
                np.save(out1,XX)
        #print site_T0_name,': AB_NanZeroRemover --> selected ',len(YEAR),' out of ',len(yy)
        return abo,bel,YEAR
#############################################
def SOC_NanZeroRemover(site_T0,site_T0_name):
#############################################
        BigRelativeError=0.15 # mean percentage variance amongst all sites
        # put year, soc, variance into numpy arrays
        yy0=np.asarray(site_T0['Year']).astype(np.int16)            # array of years
        ss0=np.asarray(site_T0['SOC']).astype(np.float64)           # array of SOCs
        vv0=np.asarray(site_T0['SOC variance']).astype(np.float64)  # array of SOC variances
        ss0=np.where(np.isnan(ss0),0,ss0)      # replace NaN with 0
        sc=ss0[ss0>0]                            # cut away all 0s, sc now corresponds to real measurements
        YEAR=yy0[ss0>0]                          # select the years corresponding to sc
        sc=sc/10.                                # pass to kgC/m2
        vv0=vv0/100.
        std2=np.std(sc)**2                      # square standard deviation of the measurements (use when no error provided - ??)
        if (std2 == 0):
                std2=(BigRelativeError*sc)**2    # <-- check  # if std == 0 use BigRelativeError in spite of std
        vv0=np.where(np.isnan(vv0),std2,vv0)   # Replace NaN in variance array with std2
        vv0=np.where(vv0==0,std2,vv0)           # Replace 0 in variance with std2
        var=vv0[ss0>0]                           # Restrict variance corresponding to the selected SOCs data
        print site_T0_name,': SOC_NanZeroRemover (cleanup of SOC data) --> selected ',len(YEAR), ' years out of ',len(yy0)
        return sc,var,YEAR

###################
#Equations
##################

def Young_SS(litterin,k1,r):
	#litterin_0 = litterin[0]
	litterin_0 = litterin
	Y_0 = litterin_0/(k1*r)	
	return Y_0

def Young_C(n_an,litterin,k1,r,Y_SS):
	t = np.arange(n_an-1)
	Y_t = litterin/(k1*r) + (Y_SS - litterin/(k1*r))*np.exp(-k1*r*t)
	Y_0_t = np.append(Y_SS,Y_t)
	return Y_0_t

def Old_SS(litterin,k2,r,h):
        #litterin_0 = litterin[0]
	litterin_0 = litterin
        O_0 = h*litterin_0/(k2*r)
        return O_0

def Old_C(n_an,litterin,k1,k2,r,h,Y_SS,O_SS):
        t = np.arange(n_an-1)
	O_t = h*litterin/(k2*r) + (O_SS - h*litterin/(k2*r)-h*(k1*r*Y_SS-litterin)/(r*(k2-k1)))*np.exp(-k2*r*t)+(h*(k1*r*Y_SS-litterin)/(r*(k2-k1)))*np.exp(-k1*r*t)
	O_0_t = np.append(O_SS,O_t)
	return O_0_t


############################################################
# J_new ######## OBJECTIVE FUNCTION
############################################################
def J_new(in_new):

        global NEW_ITER, n_an_4p1000, predict_c
        global k1_opt,k2_opt,h_opt,rstraw,hstraw,Y_SS_opt,O_SS_opt
	
	predict_Young_C = Young_C(n_an_4p1000,in_new,k1_opt,rstraw,Y_SS_opt)
	predict_Old_C = Old_C(n_an_4p1000,in_new,k1_opt,k2_opt,rstraw,hstraw,Y_SS_opt,O_SS_opt)
        predict_c = predict_Young_C+predict_Old_C
        J_new = abs(target*n_an_4p1000 - np.sum(predict_c[predict_c.shape[0]-1]-predict_c[0]))
	#print 'OBJ FUN', J_new
	#print 'predict_c',predict_c

        NEW_ITER+=1
        if ( NEW_ITER % 100 == 0 ):
                print "NEW_ITER ",NEW_ITER," in_new=",in_new," J_new=",J_new
        param_opt=in_new
        return J_new

############################################################
# J_param ######## OBJECTIVE FUNCTION
############################################################
def J_param(k1_new_k2_new_h_new):
	
        global NEW_ITER, predict_c, sc, SOC_var, SOC_all_meas
        global rstraw, hstraw, r0
	global Y_0site, O_0site

	k1_new=k1_new_k2_new_h_new[0]
	k2_new=k1_new_k2_new_h_new[1]
	h_new=k1_new_k2_new_h_new[2]
		
	litter_in = SITE_litter_inc[j]
	n_an = SITE_n_an[j]

	Y_0site = Young_SS(litter_in,k1_new,r0)
	O_0site = Old_SS(litter_in,k2_new,r0,h_new)

	
        predict_Young_C = Young_C(n_an,litter_in,k1_new,rstraw,Y_0site)
        predict_Old_C = Old_C(n_an,litter_in,k1_new,k2_new,rstraw,h_new,Y_0site,O_0site)
        predict_c = predict_Young_C+predict_Old_C

	
	ss0=np.where(np.isnan(SOC_all_meas),0,SOC_all_meas)      # replace NaN with 0
        predict_c=predict_c[ss0>0]

	Delta = predict_c - sc
	Delta2w = Delta**2/SOC_var
	J_param = np.sum(Delta2w)	


        NEW_ITER+=1
        if ( NEW_ITER % 100 == 0 ):
                print "NEW_ITER ",NEW_ITER," k2_new=",k2_new," J_new=",J_new
		print "k1_new=",k1_new
		print "k2_new=",k2_new
		print "h_new=",h_new

	        #print 'Delta', Delta
        	#print 'OBJ FUN', J_param

        	print 'sc',sc
        	print 'predict_c',predict_c

	        #print 'k2_new,h_new,litter_in'
        	#print k2_new,h_new,litter_in


        param_opt=k1_new_k2_new_h_new

        return J_param

################################################
# OPTIMIZATION CONSTRAINTS
#k1>k2
################################################
def constr1(x):
        return x[0]-x[1]

con1={'type':'ineq','fun':constr1}
cons=[con1]

###################
#Parameters (test)
##################
#Decomposition rates
k1 = 0.8
k2 = 0.00605
########################
#humification coefficient
#and enviromental parameter
#Steady state
h0 = 0.125
r0 = 1.
#Straw (without N)
hstraw = 0.125
rstraw = 1.22

#Priors
k1_k2_h_priors = np.array([k1,k2,hstraw])
bnds_param = [(0.,20.),(0.,10.),(0.,1.)] 


if(optimization_4p1000):
        ##################
        #Import optimized params
        ##################
        n=0 # contatore
        inp=open('k1_k2_h_optim6.txt','rb')
        while inp.read(1):
                inp.seek(-1,1)
                XX_optimized_param=np.load(inp)
                n+=1
	#Number of parameters optimized
	P_opt = len(XX_optimized_param[1])

###################
#Data

ROOTDIR='/Users/ebruni/Desktop/DOTTORATO/'
loc_exp = ROOTDIR+'SIMULAZIONI_ORCHIDEE_SITI/SITES_dataset_Century_new.xlsx'

C_input_exp = pd.read_excel(loc_exp)
site_names_all = C_input_exp['ID.Site'].unique()[2:len(C_input_exp)]
site_names_all = map(str, site_names_all)
N_sites_all=len(site_names_all)

site_T0_array_all=np.array(['CHNO3_Min', 'COL_T0', 'CREC3_Min', 'FEU_T0', 'JEU2_M0', 'LAJA2_Min', 'RHEU1_Min', 'RHEU2_T0','ARAZ_D0_N0', 'ULT_P0_B', 'BROAD_3_Nill', 'FOG_DwN0', 'TREV1_Min','AVRI_T1TR'])

################
#LOOP over sites

site_names = np.array(site_names_all)
N_sites = len(site_names)
SITE_litter_inc = np.zeros(N_sites)
SITE_SOC_steady_state = np.zeros(N_sites)
SITE_SOC_Young_steady_state = np.zeros(N_sites)
SITE_SOC_Old_steady_state = np.zeros(N_sites)
SITE_n_an = np.zeros(N_sites)
SITE_k1_k2_h_optim = np.zeros((N_sites,3))
SITE_Y_0_optim=np.zeros(N_sites)
SITE_O_0_optim=np.zeros(N_sites)
SITE_TOT_0_optim=np.zeros(N_sites)

site_T0_array = np.array(site_T0_array_all)

j=0
for site in site_names:
	##########################################
	print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        print "READING DATA OF SITE: ",site
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	site_df = C_input_exp[(C_input_exp['ID.Site'].values == [site])]
	site_T0_name = site_T0_array[j]
	site_T0= site_df[(site_df['ID.Treatment'].values == [site_T0_name])]

	#----------------------------------------------------------------------
	#  determine litter_inc for current site
	#----------------------------------------------------------------------
	# returns ABOVE, BELOW (kgC/m2/yr) and YEAR of measurements
	SAVE_FILE=-1
	ABOVE,BELOW,YEAR = AB_NanZeroRemover(site_T0,site_T0_name,SAVE_FILE,ROOTDIR)
	#---------------------------------------------------------------------------

	ABOVE_mean=np.mean(ABOVE)
	BELOW_mean=np.mean(BELOW)
	litter_inc = ABOVE_mean+BELOW_mean

	print 'litter inc',litter_inc

	SITE_litter_inc[j]=litter_inc


	if(optimization_param):
		#######################
        	########################
        	#Parameters optimization
        	########################
        	#===================================================
        	# SOC,VARIANCE, YEARS with nan and zero removed
        	#===================================================
        	sc,SOC_var,yy0=SOC_NanZeroRemover(site_T0,site_T0_name)
        	print sc
        	SOC_all_meas=np.asarray(site_T0['SOC']).astype(np.float64)/10. #kgC/m2 
	
		n_an = len(SOC_all_meas)
		SITE_n_an[j] = n_an
	
		#Standard dynamics over n_an years
	
		Y_0 = Young_SS(litter_inc,k1,r0)
		O_0 = Old_SS(litter_inc,k2,r0,h0)
	
		TOT_0 = Y_0+O_0
	
	        print'Steady state young',Y_0
	        print'Steady state Old',O_0
		print 'Total C steady state',TOT_0
	
		SITE_SOC_Young_steady_state[j]=Y_0
		SITE_SOC_Old_steady_state[j]=O_0
		SITE_SOC_steady_state[j]=TOT_0
			
		Y_C = Young_C(n_an,litter_inc,k1,rstraw,Y_0)
		O_C = Old_C(n_an,litter_inc,k1,k2,rstraw,hstraw,Y_0,O_0)
	
		TOT_C = Y_C+O_C
	
	
		print 'Young C ', Y_C
		print 'Old C ', O_C
		print '============='
		print 'Total C', TOT_C
	
		NEW_ITER = 0	
		########################
		#Parameters optimization
	
		opt_param=minimize(J_param, k1_k2_h_priors, method='trust-constr', bounds=bnds_param,constraints=cons, options={'disp':True})
	        k1_k2_h_optim=opt_param.x
		print "Optimum solution k1, k2, h:",k1_k2_h_optim
	
		SITE_k1_k2_h_optim[j,0] = k1_k2_h_optim[0]
		SITE_k1_k2_h_optim[j,1] = k1_k2_h_optim[1]
		SITE_k1_k2_h_optim[j,2] = k1_k2_h_optim[2]
	
		#Optimized initial conditions
		SITE_Y_0_optim[j] = Y_0site
		SITE_O_0_optim[j] = O_0site
		SITE_TOT_0_optim[j] = Y_0site+O_0site

	
	j+=1

if(optimization_param):
	print '*******************'
	print 'priors:',k1_k2_h_priors 
	for i in range(N_sites):
		print site_names[i],' optimized k2 and h: ',SITE_k1_k2_h_optim[i]
		
c=0
for site in site_names:
        print '**************'
        print site
        print '**************'
	litter_in=SITE_litter_inc[c]
	print 'LITTERIN'
	print litter_in

	Tot_SOC_SS_site = SITE_SOC_steady_state[c]
	Young_SS_site = SITE_SOC_Young_steady_state[c]
	Old_SS_site = SITE_SOC_Old_steady_state[c]

	if(optimization_param):
		Tot_SOC_SS_optim = SITE_TOT_0_optim[c]
		Young_SS_optim=SITE_Y_0_optim[c]
        	Old_SS_optim=SITE_O_0_optim[c]
		n_an =SITE_n_an[c]
		k1_opt = SITE_k1_k2_h_optim[c,0]
		k2_opt = SITE_k1_k2_h_optim[c,1]
		h_opt = SITE_k1_k2_h_optim[c,2]

		#Optimized dynamics over n_an years

		Opt_Y_C = Young_C(n_an,litter_in,k1_opt,rstraw,Young_SS_optim)
		Opt_O_C = Old_C(n_an,litter_in,k1_opt,k2_opt,rstraw,h_opt,Young_SS_optim,Old_SS_optim)
		Opt_C = Opt_Y_C+Opt_O_C

		#Standard dynamics over n_an years

		Opt_Y_st = Young_C(n_an,litter_in,k1,rstraw,Young_SS_site)
		Opt_O_st = Old_C(n_an,litter_in,k1,k2,rstraw,hstraw,Young_SS_site,Old_SS_site)
		Opt_st = Opt_Y_st+Opt_O_st

		#Measured SOC
		site_df = C_input_exp[(C_input_exp['ID.Site'].values == [site])]
		site_T0_name = site_T0_array[c]
		site_T0= site_df[(site_df['ID.Treatment'].values == [site_T0_name])]

		sc,SOC_var,yy0=SOC_NanZeroRemover(site_T0,site_T0_name) 

		
		#Select only available data
		SOC_all_meas=np.asarray(site_T0['SOC']).astype(np.float64)/10. #kgC/m2
		ss0=np.where(np.isnan(SOC_all_meas),0,SOC_all_meas)      # replace NaN with 0
		Opt_C_select=Opt_C[ss0>0]
		Opt_st_select = Opt_st[ss0>0]

		print 'Years: '
		print yy0
		
		print 'Standard SOC: '
		print Opt_st

		print 'Optimized SOC: '
		print Opt_C
		#print Opt_C_select
		
		print 'Measured SOC: '
		print sc
		#Save output
		#prova
		#Years, Standard, Optimized, Measured SOC
		SOCXX = np.stack((yy0,Opt_st_select,Opt_C_select,sc))
		np.save(out1,SOCXX)

	####################
        #4p1000 optimization
        ####################
	if(optimization_4p1000):
		tstart = time.time()
		n_an_4p1000=30.
                NEW_ITER = 0
	
		k1_opt=XX_optimized_param[c,0]
		k2_opt=XX_optimized_param[c,1]
		h_opt=XX_optimized_param[c,2]
		
		print 'optimized k1,k2,h:', XX_optimized_param[c,:]
		
		#Steady state with optimized param
		Y_SS_opt = Young_SS(litter_in,k1_opt,r0)
                O_SS_opt = Old_SS(litter_in,k2_opt,r0,h_opt)
                Tot_SS = Y_SS_opt+O_SS_opt


                #Set the target
                target = Tot_SS*0.004

                #Priors
                #in_opt = np.array(litter_inc[0])*(1+0.004)
                in_opt = np.array(litter_in)*(1+0.004)

                #Bounds
                bnds=[(0,10)]

                #Optimization algorithm
                opt_mean=minimize(J_new, in_opt, method='SLSQP',bounds=bnds, options={'disp':True})
                litter_opt = opt_mean.x

                ############################
                #Check the target is reached

                END=predict_c.shape[0]-1
                C_fin=predict_c[END]
                C_init=predict_c[0]
                SUMO=C_fin-C_init
                print 'Cfin',C_fin
                print 'C_init',C_init
                Target_reached=(C_fin-C_init)/(C_init*n_an_4p1000)
                if (Target_reached < 0.005 and Target_reached > 0.003):
                        print "Target reached successfully"
                        print "Target reached :", Target_reached
                else:
                        print "Target not reached", Target_reached

                #Some prints
                print 'initial SOC (kgC/m2):',Tot_SS
                print 'initial litter (kgC/m2/yr):',litter_in
                print "litter opt (kgC/m2/yr):",litter_opt

	c+=1


if(optimization_param):
	out3=open("k1_k2_h_optim6.txt","wb")
	np.save(out3,SITE_k1_k2_h_optim)
	out3.close()
	out1.close()


###################
#End 
if(optimization_4p1000):
	tend = time.time()
	tot_time=tend-tstart
	print " "
	print " Ok, done. Total time: ",tot_time
