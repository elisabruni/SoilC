############################################################
#IBCM (Andren and Katterer 1997)
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
###################
#Data

ROOTDIR='/Users/ebruni/Desktop/DOTTORATO/'
loc_exp = ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/SITES_dataset_multimodel.xlsx'
OUTPUT_files=ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OUTPUTS_4p1000_v8/ICBM/'

#Save data
if(optimization_param):
        out1=open("SOC_ICBM_optim.txt","wb")

if(optimization_4p1000):
        out_mo_pools=open(OUTPUT_files+"SOC_ICBM_pools.txt","wb")
        out_mo=open(OUTPUT_files+"SOC_ICBM.txt","wb")
        out_lit=open(OUTPUT_files+"Litter_income_ICBM.txt","wb")
	out_lit_tot = open(OUTPUT_files+"Litter_tot_ICBM.txt","wb")
        out_priors = open(OUTPUT_files+"priors_and_opt_in_ICBM.txt","wb")
		

##################
#Functions
##################
#########################################################
def AB_NanZeroRemover(site_T0,site_T0_name):
#########################################################
        #######################################################################
        # remove NaNs and ZEROs from ABOVE and BELOW
        #######################################################################
        yy=np.asarray(site_T0['Year']).astype(np.int16)             # array of years
        aa=np.asarray(site_T0['ABOVE']).astype(np.float64)/10.  # array of C_above (kgC/m2/yr)
        bb=np.asarray(site_T0['BELOW']).astype(np.float64)/10.  # array of C_below (kgC/m2/yr)
        aa0=np.where(np.isnan(aa),999,aa)  # replace NaN with 999 in "above"
        YEAR=yy[aa0<999]                     # select years where above>0
        abo=aa[aa0<999]                      # select ABOVE>0
        bel=bb[aa0<999]                      # select corresponding BELOW
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
	t = np.arange(1,n_an+1)
	Y_t = litterin/(k1*r) + (Y_SS - litterin/(k1*r))*np.exp(-k1*r*t)
	return Y_t

def Old_SS(litterin,k2,r,h):
        #litterin_0 = litterin[0]
	litterin_0 = litterin
        O_0 = h*litterin_0/(k2*r)
        return O_0

def Old_C(n_an,litterin,k1,k2,r,h,Y_SS,O_SS):
        t = np.arange(1,n_an+1)
	O_t = h*litterin/(k2*r) + (O_SS - h*litterin/(k2*r)-h*(k1*r*Y_SS-litterin)/(r*(k2-k1)))*np.exp(-k2*r*t)+(h*(k1*r*Y_SS-litterin)/(r*(k2-k1)))*np.exp(-k1*r*t)
	return O_t


############################################################
# J_new ######## OBJECTIVE FUNCTION
############################################################
def J_new(in_new):

        global NEW_ITER, n_an_4p1000, predict_c
        global k1_opt,k2_opt,r_opt,rstraw,hstraw,Y_SS_opt,O_SS_opt,r_ss
	global litter_in, h_pred

	#PROVA
	#Calculate proportion of in_new that is straw / OM
	prop_straw = np.minimum(litter_in/in_new,1.)
	h_pred = hstraw*prop_straw+h_fym*(1.-prop_straw)
	##

	predict_Young_C = Young_C(n_an_4p1000,in_new,k1_opt,r_opt,Y_SS_opt)
	predict_Old_C = Old_C(n_an_4p1000,in_new,k1_opt,k2_opt,r_opt,h_pred,Y_SS_opt,O_SS_opt)
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
def J_param(k1_new_k2_new_r_new):
	
        global NEW_ITER, predict_c, sc, SOC_var, SOC_all_meas
        global rstraw, hstraw, r0
	global Y_0site, O_0site

	k1_new=k1_new_k2_new_r_new[0]
	k2_new=k1_new_k2_new_r_new[1]
	r_new=k1_new_k2_new_r_new[2]
		
	litter_in = SITE_litter_inc[j]
	n_an = SITE_n_an[j]

	Y_0site = Young_SS(litter_in,k1_new,r0)
	O_0site = Old_SS(litter_in,k2_new,r0,hstraw)

	
        predict_Young_C = Young_C(n_an,litter_in,k1_new,r_new,Y_0site)
        predict_Old_C = Old_C(n_an,litter_in,k1_new,k2_new,r_new,hstraw,Y_0site,O_0site)
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
		print "r_new=",r_new

	        #print 'Delta', Delta
        	#print 'OBJ FUN', J_param

        	print 'sc',sc
        	print 'predict_c',predict_c


        param_opt=k1_new_k2_new_r_new

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
#r0 = 1.
#Straw (without N)
hstraw = 0.125
h_fym=0.31
h_ss=0.47

#rstraw_noN = 1.22
#rstraw_N = 1.

#r_fym = 1.1
#r_ss = 0.97

#Priors
#k1_k2_h_priors = np.array([k1,k2,hstraw])
#bnds_param = [(0.,20.),(0.,10.),(0.,1.)] 


if(optimization_4p1000):
        ##################
        #Import optimized params
        ##################
        n=0 # contatore
        inp=open(ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM8/optim_param_ICBM.txt','rb')
        while inp.read(1):
                inp.seek(-1,1)
                XX_optimized_param=np.load(inp)
                n+=1
	#Number of parameters optimized
	P_opt = len(XX_optimized_param[1])

        ##################
        #Import temp and moist response function
        #################
        #r_day forward
        n=0
        inp=open(ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM8/temp_moist_re_day_av_fw_ICBM_4p1000.txt','rb')
        while inp.read(1):
                inp.seek(-1,1)
                XX_temp_moist_resp_func_fw=np.load(inp)
                n+=1
        #################
        #re_day steady state
        n=0
        inp=open(ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM8/temp_moist_re_day_av_ss_ICBM_4p1000.txt','rb')
        while inp.read(1):
                inp.seek(-1,1)
                XX_temp_moist_resp_func_ss=np.load(inp)
                n+=1


if(optimization_param):
        ##################
        #Import temp and moist response function
        #################
	#r_day forward
	n=0
	inp=open(ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM8/temp_moist_re_day_av_fw_ICBM.txt','rb')
        while inp.read(1):
                inp.seek(-1,1)
                XX_temp_moist_resp_func_fw=np.load(inp)
                n+=1

        #################
	#re_day steady state
        n=0
        inp=open(ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM8/temp_moist_re_day_av_ss_ICBM.txt','rb')
        while inp.read(1):
                inp.seek(-1,1)
                XX_temp_moist_resp_func_ss=np.load(inp)
                n+=1


##############
C_input_exp = pd.read_excel(loc_exp)
site_names_all = C_input_exp['ID.Site'].unique()[2:len(C_input_exp)]
site_names_all = map(str, site_names_all)
N_sites_all=len(site_names_all)

site_T0_array_all=np.array(['CHNO3_Min', 'COL_T0', 'CREC3_Min', 'FEU_T0', 'LAJA2_Min', 'RHEU1_Min', 'RHEU2_T0','ARAZ_D0_N0', 'ULT_P0_B', 'BROAD_3_Nill',  'TREV1_Min','AVRI_T1TR','BOLO_T0', 'GRAB_CP','MUNCHE_CP','RITZ_CP'])

N_use_control = np.array([1,0,1,0,1,1,0,0,0,0,1,0,0,1,0,0])

################
#LOOP over sites

site_names = np.array(site_names_all)
N_sites = len(site_names)
SITE_litter_inc = np.zeros(N_sites)
SITE_ABOVE = []
SITE_BELOW = []
SITE_ABOVE_mean = np.zeros(N_sites)
SITE_BELOW_mean = np.zeros(N_sites)
SITE_ERR2_ABOVE_mean = np.zeros(N_sites)
SITE_ERR2_BELOW_mean = np.zeros(N_sites)
SITE_cov_mean = []
SITE_SOC_steady_state = np.zeros(N_sites)
SITE_SOC_Young_steady_state = np.zeros(N_sites)
SITE_SOC_Old_steady_state = np.zeros(N_sites)
SITE_n_an = np.zeros(N_sites)
SITE_k1_k2_r_optim = np.zeros((N_sites,3))
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
	ABOVE,BELOW,YEAR = AB_NanZeroRemover(site_T0,site_T0_name)
	#---------------------------------------------------------------------------

        SITE_ABOVE.append(ABOVE)
        SITE_BELOW.append(BELOW)

	ABOVE_mean=np.mean(ABOVE)
	BELOW_mean=np.mean(BELOW)
        SITE_ABOVE_mean[j]=ABOVE_mean
        SITE_BELOW_mean[j]=BELOW_mean

	litter_inc = ABOVE_mean+BELOW_mean
	print 'litter inc',litter_inc
	SITE_litter_inc[j]=litter_inc

        SITE_ERR2_ABOVE_mean[j]=np.std(ABOVE)**2/len(ABOVE)
        SITE_ERR2_BELOW_mean[j]=np.std(BELOW)**2/len(BELOW)
        if (SITE_ERR2_ABOVE_mean[j] == 0):
                SITE_ERR2_ABOVE_mean[j]=0.05 # to be checked
                SITE_ERR2_BELOW_mean[j]=0.05
	
	cov_AB_mean = np.cov(ABOVE,BELOW)/np.sqrt(len(ABOVE))
	SITE_cov_mean.append(cov_AB_mean)

	#if(N_use_control[j]==0):
	#	rstraw = rstraw_noN
	#else:
	#	rstraw = rstraw_N

	#r parameter calculated from temp and moist response function
	r0 = XX_temp_moist_resp_func_ss[j]
	rstraw = XX_temp_moist_resp_func_fw[j]

	#Priors
	#TEST
        #k1=0.032622
	k1_k2_r_priors = np.array([k1,k2,rstraw])
	bnds_param = [(0.,20.),(0.,10.),(-2.,10.)]

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
		#TEST (munche  doesn't work with constraints, but anyways they are respected)
		if(j==14):	
			opt_param=minimize(J_param, k1_k2_r_priors, method='trust-constr', bounds=bnds_param, options={'disp':True})
		else:
			opt_param=minimize(J_param, k1_k2_r_priors, method='trust-constr', bounds=bnds_param,constraints=cons, options={'disp':True,'maxiter':5000})

	        k1_k2_r_optim=opt_param.x
		print "Optimum solution k1, k2, r:",k1_k2_r_optim
	
		SITE_k1_k2_r_optim[j,0] = k1_k2_r_optim[0]
		SITE_k1_k2_r_optim[j,1] = k1_k2_r_optim[1]
		SITE_k1_k2_r_optim[j,2] = k1_k2_r_optim[2]
	
		#Optimized initial conditions

        	Y_0site = Young_SS(litter_inc,k1_k2_r_optim[0],r0)
        	O_0site = Old_SS(litter_inc,k1_k2_r_optim[1],r0,hstraw)
		SITE_Y_0_optim[j] = Y_0site
		SITE_O_0_optim[j] = O_0site
		SITE_TOT_0_optim[j] = Y_0site+O_0site

	
	j+=1

if(optimization_param):
	print '*******************'
	print 'priors:',k1_k2_r_priors 
	for i in range(N_sites):
		print site_names[i],' optimized k2 and h: ',SITE_k1_k2_r_optim[i]

optim_4p1000_sites = np.zeros(N_sites)		
c=0
for site in site_names:
        print '**************'
        print site
        print '**************'
	ABOVE=SITE_ABOVE[c]
	BELOW=SITE_BELOW[c]
	ABOVE_mean=SITE_ABOVE_mean[c]
	BELOW_mean = SITE_BELOW_mean[c]
        err2_above = SITE_ERR2_ABOVE_mean[c]
        err2_below = SITE_ERR2_BELOW_mean[c]
	#Approx CHECK if better
	litter_in_err = np.sqrt(err2_above+err2_below)

	litter_in=SITE_litter_inc[c]
	print 'LITTERIN'
	print litter_in

	Tot_SOC_SS_site = SITE_SOC_steady_state[c]
	Young_SS_site = SITE_SOC_Young_steady_state[c]
	Old_SS_site = SITE_SOC_Old_steady_state[c]

        #r parameter calculated from temp and moist response function
        r0 = XX_temp_moist_resp_func_ss[c]
        rstraw = XX_temp_moist_resp_func_fw[c]

	if(optimization_param):
		Tot_SOC_SS_optim = SITE_TOT_0_optim[c]
		Young_SS_optim=SITE_Y_0_optim[c]
        	Old_SS_optim=SITE_O_0_optim[c]
		n_an =SITE_n_an[c]
		k1_opt = SITE_k1_k2_r_optim[c,0]
		k2_opt = SITE_k1_k2_r_optim[c,1]
		r_opt = SITE_k1_k2_r_optim[c,2]

		#Optimized dynamics over n_an years

		Opt_Y_C = Young_C(n_an,litter_in,k1_opt,r_opt,Young_SS_optim)
		Opt_O_C = Old_C(n_an,litter_in,k1_opt,k2_opt,r_opt,hstraw,Young_SS_optim,Old_SS_optim)
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
		#Save output in tC/ha
		#Years, Standard, Optimized, Measured SOC, Error
		SOCXX = np.stack((yy0,Opt_st_select*10.,Opt_C_select*10.,sc*10.,SOC_var*100.)) #tC/ha
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
		r_opt=XX_optimized_param[c,2]
		
		print 'optimized k1,k2,r:', XX_optimized_param[c,:]
		
		#Steady state with optimized param
		Y_SS_opt = Young_SS(litter_in,k1_opt,r0)
                O_SS_opt = Old_SS(litter_in,k2_opt,r0,hstraw)
                Tot_SS = Y_SS_opt+O_SS_opt

		#Dynamics with optimized param
                Y_c_opt = Young_C(n_an_4p1000,litter_in,k1_opt,r_opt,Y_SS_opt)
                O_c_opt = Old_C(n_an_4p1000,litter_in,k1_opt,k2_opt,r_opt,hstraw,Y_SS_opt,O_SS_opt)
                opt_c = Y_c_opt+O_c_opt


		print 'Difference between equilibrium and 1st year exp (Dyn[0]):'
		print 'spinup:',Tot_SS
		print 'Dyn[0]:',opt_c[0]

                #Set the target
                target = opt_c[0]*0.004

                #Priors
                #in_opt = np.array(litter_inc[0])*(1+0.004)
                in_opt = np.array(litter_in)*(1+0.004)

                #Bounds
                bnds=[(0,10)]

                #Optimization algorithm
                opt_mean=minimize(J_new, in_opt, method='SLSQP',bounds=bnds, options={'disp':True})
                litter_opt = opt_mean.x

		optim_4p1000_sites[c]=litter_opt

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
		print "***********************"
		print "4p1000 OPTIMIZATION"
		print "***********************"
                print 'initial SOC (kgC/m2):',opt_c[0]
                print 'initial litter (kgC/m2/yr):',litter_in
                print "litter opt (kgC/m2/yr):",litter_opt
		print "Litter increase (%):", (litter_opt-litter_in)*100./litter_in
		print "***********************"

		############################################################
		#CREATE output file with PREDICTED CARBON over 30 years (standard and 4x100)
		#                       + EXPERIMENT SOC filled with Nan for 30 years (T0 + other treatments)
		###########################################################
		predict_c_standard_Y=Young_C(n_an_4p1000,litter_in,k1_opt,r_opt,Y_SS_opt)
		predict_c_standard_O=Old_C(n_an_4p1000,litter_in,k1_opt,k2_opt,r_opt,hstraw,Y_SS_opt,O_SS_opt)
		predict_c_standard_pools=np.stack((predict_c_standard_Y,predict_c_standard_O),axis=1)
		predict_c_standard = predict_c_standard_Y+predict_c_standard_O

		#4p1000
		predict_c_opt_Y=Young_C(n_an_4p1000,litter_opt,k1_opt,r_opt,Y_SS_opt)
                predict_c_opt_O=Old_C(n_an_4p1000,litter_opt,k1_opt,k2_opt,r_opt,h_pred,Y_SS_opt,O_SS_opt)
		predict_c_opt_pools = np.stack((predict_c_opt_Y,predict_c_opt_O),axis=1)
		predict_c_opt=predict_c_opt_Y+predict_c_opt_O

		year_out = np.arange(1,31)

		SOC_pools_out=np.stack((predict_c_standard_pools*10.,predict_c_opt_pools*10.))#tC/ha
		SOC_out=np.stack((year_out,predict_c_standard*10.,predict_c_opt*10.)) #tC/ha

		np.save(out_mo_pools,SOC_pools_out)
		np.save(out_mo,SOC_out)

		############################
		#UNCERTAINTIES
		###########################
		Uncert_Q = True
		if(Uncert_Q):

			MC_length=2 #set the number of Monte Carlo simulations

			#optimize for n variations of in_opt generatated randomly around the above/below covariance
			opt_parameters_MC=np.zeros(MC_length)

			#prior above_below
			AB_BE_array=np.array([ABOVE_mean,BELOW_mean])
			ab_be_init_est=AB_BE_array*(1+0.004)

			cov_AB_mean=np.cov(ABOVE,BELOW)/np.sqrt(len(ABOVE)) #covariance between ABOVE mean and BELOW mean if no Nans
			print 'cov',cov_AB_mean
			if (np.all(cov_AB_mean)==0): #if covariance is 0, take the mean amongst all sites' covariances to generate cov_AB_mean
				cov_AB_mean=np.mean(SITE_cov_mean)
			#print cov_AB_mean

			in_rand_param_MC = np.zeros(MC_length)
			sample_shape=0
			while sample_shape<MC_length: #generate random multinormal until sample_shape=MC_length
				print ' '
				print '************************'
				print 'UNCERTAINTIES ITERATION:'
				print sample_shape+1
				in_rand_gen=np.random.multivariate_normal(ab_be_init_est,cov_AB_mean)
				if all(i>0 for i in in_rand_gen): #test if all elements of random array are positive
					#Only positive arrays
					in_rand=in_rand_gen
					in_rand_param=np.sum(in_rand)  #array to optimize
					#Save all priors in an array
					in_rand_param_MC[sample_shape]=in_rand_param #add new generated sample to array on rand_in samples
					print '****************'
					print "LITTER IN generated randomly from ab_be_init_est:", in_rand
					print 'total litter in (kgC/m2/yr):',in_rand_param
					print '****************'
					#Minimize J_new for the generated array
					opt=minimize(J_new, in_rand_param, method='SLSQP',bounds=bnds, options={'disp':True})
					opt_param=opt.x
					opt_parameters_MC[sample_shape]=opt_param
					
					print "OPTIMUM LITTER IN (kgC/m2/yr)"
					print opt_param
					print '****************'
					print "Litter in increase (%)"
					print (opt_param-in_rand_param)*100./in_rand_param
					print '****************'
					sample_shape+=1


			#matrix of the optimum parameters
			#print "opt parameters:", opt_parameters_MC

			#save priors (litterin generated random) and out (optimized)
			out_priors_and_opt = np.stack((in_rand_param_MC*10.,opt_parameters_MC*10.)) #save as tC/ha/yr
			np.save(out_priors,out_priors_and_opt)

			#STANDARD ERROR CALCULATION for the optimized litter inputs to reach 4x1000

			litter_opt_err=np.std(opt_parameters_MC[:])/np.sqrt(MC_length)

			#print "Uncertainties:",litter_opt_err


			#Error litter = SE per litter in e in opt
			#print litter_in*10.,litter_in_err*10.,litter_opt*10.,litter_opt_err*10.
			save_lit = np.stack((litter_in*10.,litter_in_err*10.,litter_opt[0]*10.,litter_opt_err*10.))#save as tC/ha/yr
			np.save(out_lit,save_lit)		

			#in ICBM save_lit = litter_tot_save
			litter_tot_save = save_lit
			np.save(out_lit_tot,litter_tot_save)

	c+=1


if(optimization_param):
	out3=open("optim_param_ICBM.txt","wb")
	np.save(out3,SITE_k1_k2_r_optim)
	out3.close()
	out1.close()


###################
#End 
if(optimization_4p1000):
	out_lit_tot.close()
	out_lit.close()
	out_priors.close()
	out_mo_pools.close()
        out_mo.close()
        for i in range(N_sites):
                print site_names_all[i],' C inputs to 4p1000: ', optim_4p1000_sites[i]
	tend = time.time()
	tot_time=tend-tstart
	print " "
	print " Ok, done. Total time: ",tot_time
