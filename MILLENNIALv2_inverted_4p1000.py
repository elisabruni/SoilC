#Millennial_4p1000_inverted.py v2
###########################################################
#Millennial Abramof et al
#see https://github.com/PNNL-TES/millenial/blob/master/R/decomp.R
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
np.set_printoptions(threshold=sys.maxsize)
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


#DIRECTORIES
ROOTDIR='/Users/ebruni/Desktop/DOTTORATO/'
loc_out = ROOTDIR+"SCRIPT_MODELLI/MULTIMODEL/OUTPUTS_4p1000_v5/MILLENNIALv2/"
loc_exp = ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/SITES_dataset_multimodel.xlsx'
loc_optim_param = ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM5/'
loc_spinup = ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM5/'


#OUTPUT FILES
out_mo_pools=open(loc_out+"SOC_MILLENNIALv2_pools.txt","wb")
out_mo=open(loc_out+"SOC_MILLENNIALv2.txt","wb")
out_lit=open(loc_out+"Litter_income_MILLENNIALv2.txt","wb")
out_lit_tot=open(loc_out+"Litter_tot_MILLENNIALv2.txt","wb")
out_priors = open(loc_out+"priors_and_opt_in_MILLENNIALv2.txt","wb")

#########################################################
def DailyToYearly(filein,typeConv):
#########################################################
		dim_file = len(filein.shape)
                c=0
                x=0
                #Set number of years
		if(dim_file==1):
                	var_year=np.zeros(np.int(len(filein)/365))
                	lent = np.int(len(var_year))
                	for j in range(0,lent):
                	        if(typeConv==0.): #for precipitation and PET (cumulative)
                	                var_year[c]=np.sum(filein[x:365+x])
                	        else: #for Temperature (mean)
                	                var_year[c]=np.sum(filein[x:365+x])/365
                	        x=x+365
                	        c+=1
		else:
			var_year=np.zeros((np.int(filein.shape[0]/365),filein.shape[1]))
			lent = np.int(len(var_year[0]))
			for j in range(0,lent):
                                if(typeConv==0.): #for precipitation and PET (cumulative)
                                        var_year[c,:]=np.sum(filein[x:365+x,:],axis=0)
                                else: #for Temperature (mean)
                                        var_year[c,:]=np.sum(filein[x:365+x,:],axis=0)/365
                                x=x+365
                                c+=1
                return var_year

#########################################################
def AB_NanZeroRemover(site_T0,site_T0_name,iout,ROOTDIR):
#########################################################
        if ( iout > 0 ):
                out1=open(ROOTDIR+"AB.data","wb")
        #######################################################################
        # remove NaNs and ZEROs from ABOVE and BELOW
        #######################################################################
        yy=np.asarray(site_T0['Year']).astype(np.int16)             # array of years
        aa=np.asarray(site_T0['ABOVE']).astype(np.float64)*100/365  # array of C_above (gC/m2/day)
        bb=np.asarray(site_T0['BELOW']).astype(np.float64)*100/365  # array of C_below (gC/m2/day)
        aa0=np.where(np.isnan(aa),999,aa)  # replace NaN with 999 in "above"
        YEAR=yy[aa0<999]                     # select years where above>0
        abo=aa[aa0<999]                      # select ABOVE>0
        bel=bb[aa0<999]                      # select corresponding BELOW
        if (iout > 0):
                XX=np.stack((YEAR,abo,bel),axis=0)
                np.save(out1,XX)
        print site_T0_name,': AB_NanZeroRemover --> selected ',len(YEAR),' out of ',len(yy)
        return abo,bel,YEAR

#############################################
def SOC_NanZeroRemover(site_T0,site_T0_name):
#############################################
        BigRelativeError=0.15 # mean percentage variance amongst all sites
        # put year, soc, variance into numpy arrays
        yy0=np.asarray(site_T0['Year']).astype(np.int16)            # array of years
        ss0=np.asarray(site_T0['SOC']).astype(np.float64)*100           # array of SOCs (pass to gC/m2)
        vv0=np.asarray(site_T0['SOC variance']).astype(np.float64)  # array of SOC variances
        ss0=np.where(np.isnan(ss0),0,ss0)      # replace NaN with 0
        sc=ss0[ss0>0]                            # cut away all 0s, sc now corresponds to real measurements
        ########################################
        #No need to convert sc and vv0 in rothc
        #######################################
        YEAR=yy0[ss0>0]                          # select the years corresponding to sc
        vv0=vv0*10000
        std2=np.std(sc)**2                      # square standard deviation of the measurements (use when no error provided - ??)
        if (std2 == 0):
                std2=(BigRelativeError*sc)**2    # <-- check  # if std == 0 use BigRelativeError in spite of std
        vv0=np.where(np.isnan(vv0),std2,vv0)   # Replace NaN in variance array with std2
        vv0=np.where(vv0==0,std2,vv0)           # Replace 0 in variance with std2
        var=vv0[ss0>0]                           # Restrict variance corresponding to the selected SOCs data
        print site_T0_name,': SOC_NanZeroRemover (cleanup of SOC data) --> selected ',len(YEAR), ' years out of ',len(yy0)
        return sc,var,YEAR,ss0
    
############################################################
# J_new ######## OBJECTIVE FUNCTION
############################################################
def J_new(in_new):
        global NEW_ITER, target
        global n_an_4p1000,n_an
        global SOC_data_init, predict_c, SOC_stable
        global clay,temp_t,water_t 
        global state, t_fwd

        
	predict_c = odeint(decomp,state,t_fwd,args=(temp_t,water_t,in_new))

        J_new = abs(target*n_an_4p1000 - np.sum(predict_c[predict_c.shape[0]-1]-predict_c[0]))

        print 'OBJ FUN', J_new
        #print 'predict_c',predict_c

        NEW_ITER+=1
        if ( NEW_ITER % 100 == 0 ):
                print "NEW_ITER ",NEW_ITER," in_new=",in_new," J_new=",J_new
        param_opt=in_new
        return J_new


########################
#Decomposition function
#########################

def decomp(state, t_now, forc_st, forc_sw, forc_npp):

        global param_pi, param_pa, kaff_pl, alpha_pl
        global eact_pl,rate_pa, rate_break, rate_leach
        global kaff_des, param_p1, param_p2, kaff_lb, alpha_lb
	global eact_lb,rate_bd, rate_ma, cue_ref
	global cue_t, tae_ref, matpot, lambda_p
	global porosity, kamin, param_pb, param_c1
	global clay, BD, pH

        #print 'DECOMP'

        POM = state[0]
        LMWC = state[1]
        AGG = state[2]
        MIC = state[3]
        MAOM = state[4]


        t_now=np.int(t_now)


    	#Equation 11
    	#Mayes 2012, SSAJ
    	kaff_lm = np.exp(-param_p1 * 7. - param_p2) * kaff_des

    	#Equation 12
    	#Georgiou in review
    	param_qmax = param_bulkd * param_c1 * param_claysilt

	# Hydrological properties
    
    	#Equation 5
    	scalar_wd = (forc_sw[t_now]/ porosity)**0.5
    
   	#Equation 4
    	scalar_wb = np.exp(lambda_p * -matpot) * (kamin + (1 - kamin) * ((porosity - forc_sw[t_now]) / porosity)**0.5) * scalar_wd

	# Decomposition
    
	gas_const = 8.31446
	
	#Equation 3
	vmax_pl = alpha_pl * np.exp(-eact_pl / (gas_const * (forc_st[t_now] + 273.15)))
	
	#Equation 2
	# POM -> LMWC
	if(POM>0 and MIC>0):
		f_PO_LM = vmax_pl * scalar_wd * POM * MIC / (kaff_pl + MIC)
	else:
		f_PO_LM=0

	#Equation 6
	# POM -> AGG
	if(POM>0):
		f_PO_AG = rate_pa * scalar_wd * POM
	else:
		f_PO_AG=0
	
	#Equation 7
	# AGG -> MAOM
	if(AGG>0):
		f_AG_break = rate_break * scalar_wd * AGG
	else:
		f_AG_break=0
	
	#Equation 9
	# LMWC -> out of system leaching
	if(LMWC>0):
		f_LM_leach = rate_leach * scalar_wd * LMWC
	else:
		f_LM_leach=0

	#Equation 10
	# LMWC -> MAOM
	if(LMWC>0 and MAOM>0):
		f_LM_MA = scalar_wd * kaff_lm * LMWC * (1 - MAOM / param_qmax)
	else:
		f_LM_MA=0
	
	#Equation 13
	# MAOM -> LMWC
	if(MAOM>0):
		f_MA_LM = kaff_des * MAOM / param_qmax
	else:
		f_MA_LM=0
	
	#Equation 15
	vmax_lb = alpha_lb * np.exp(-eact_lb / (gas_const * (forc_st[t_now] + 273.15)))
	
	#Equation 14
	# LMWC -> MIC
	if(LMWC>0 and MIC>0):
		f_LM_MB = vmax_lb * scalar_wb * MIC * LMWC / (kaff_lb + LMWC)
	else:
		f_LM_MB=0

	#Equation 16
	# MIC -> MAOM/LMWC
	if(MIC>0):
		f_MB_turn = rate_bd * MIC**2.0
	else:
		f_MB_turn=0

	#Equation 18
	# MAOM -> AGG
	if(MAOM>0):
		f_MA_AG = rate_ma * scalar_wd * MAOM
	else:
		f_MA_AG=0
	
	#Equation 22
	# microbial growth flux, but is not used in mass balance
	
	#Equation 21
	# MIC -> atmosphere
	if(MIC>0 and LMWC>0):
		f_MB_atm = f_LM_MB * (1 - (cue_ref- cue_t * (forc_st[t_now] - tae_ref) ) )
	else:
		f_MB_atm=0

	# Update state variables
    
	#Equation 1
	dPOM = forc_npp * param_pi + f_AG_break * param_pa - f_PO_AG - f_PO_LM
	
	#Equation 8
	dLMWC = forc_npp * (1. - param_pi) - f_LM_leach + f_PO_LM - f_LM_MA - f_LM_MB + f_MB_turn * (1. - param_pb) + f_MA_LM
	
	#Equation 17
	dAGG = f_MA_AG + f_PO_AG - f_AG_break
	
	#Equation 19
	dMAOM = f_LM_MA - f_MA_LM + f_MB_turn * param_pb - f_MA_AG + f_AG_break * (1. - param_pa)
	
	#Equation 20
	dMIC = f_LM_MB - f_MB_turn - f_MB_atm

        return [dPOM, dLMWC, dAGG, dMIC, dMAOM]

############################################################

############################################################
#Initialization
tstart = time.time()
NEW_ITER = 0
CHI2_PRINT_FREQUENCY=2

#########################
######################################################
#For each site: Set SITE name and experiment duration
#####################################################


C_input_exp = pd.read_excel(loc_exp)
site_names_all = C_input_exp['ID.Site'].unique()[2:len(C_input_exp)]
site_names_all = map(str, site_names_all)

N_sites_all=len(site_names_all)


#Control plots
site_T0_array_all=np.array(['CHNO3_Min', 'COL_T0', 'CREC3_Min', 'FEU_T0', 'JEU2_M0', 'LAJA2_Min', 
                            'RHEU1_Min', 'RHEU2_T0','ARAZ_D0_N0', 'ULT_P0_B', 'BROAD_3_Nill', 
                            'TREV1_Min','AVRI_T1TR','BOLO_T0', 'GRAB_CP','MUNCHE_CP','RITZ_CP'])

#Stationary solution array for each experiment T0
N_sites=len(site_names_all)

print N_sites

#Read numerical spinup

SITES_spinup=pd.read_csv(loc_spinup+'list_SOC_spinup_opti_MILv2.csv')
SITES_spinup.drop(SITES_spinup.columns[[0,1]], axis = 1, inplace = True)

#Import optimized parameters
my_cols = ['param','value']
SITES_optim_param=pd.read_csv(loc_optim_param+'optim_param_MILLENNIALv2.csv', names=my_cols, engine="python")

SITE_eact_lb=SITES_optim_param['value'][(SITES_optim_param['param']=='eact_lb')]
SITE_eact_pl=SITES_optim_param['value'][(SITES_optim_param['param']=='eact_pl')]
SITE_kaff_pl=SITES_optim_param['value'][(SITES_optim_param['param']=='kaff_pl')]


SOC_exp_array=[] #interpolated SOC dynamics experiments
SOC_clean_exp_array=[]
SOC_clean_exp_variance=[]
SOC_clean_year=[]
SOC_clean_where=[]

SITE_SOC_data_init=np.zeros(N_sites)
SITE_year0=np.zeros(N_sites)
SITE_date_init=np.zeros(N_sites)
SITE_date_end=np.zeros(N_sites)
SITE_exper_len = np.zeros(N_sites)
SITE_clay=np.zeros(N_sites)
SITE_silt=np.zeros(N_sites)
SITE_litterinc = np.zeros(N_sites)
SITE_litterinc_err2 = np.zeros(N_sites)
SITE_water_in=[]
SITE_temp_in=[]
SITE_water_t=[]
SITE_temp_t=[]
SITE_ABOVE_mean=np.zeros(N_sites)
SITE_BELOW_mean=np.zeros(N_sites)
SITE_ERR2_ABOVE_mean=np.zeros(N_sites)
SITE_ERR2_BELOW_mean=np.zeros(N_sites)
SITE_cov_mean= []
SITE_TREATMENTS=[]
SITE_BD=np.zeros(N_sites)
SITE_pH=np.zeros(N_sites)
SITE_n_an = np.zeros(N_sites)
                              
                              
                              
j=0
for site in site_names_all:
        ##########################################
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        print "READING DATA OF SITE: ",site
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        site_df = C_input_exp[(C_input_exp['ID.Site'].values == [site])]
        year_0 = np.min(site_df['Year'])
        year_end = np.max(site_df['Year'])
        year_30=year_0+30
        missing_years_to30 = np.int(year_30-year_end-1)

        site_T0_name = site_T0_array_all[j]
        site_treatments = site_df['ID.Treatment'].unique()[0:len(site_df)]
        site_treatments = map(str, site_treatments)
                              
        #INTERPOLATE each treatment from year0 to year30, fill with Nan
        TREATMENTS=np.zeros((30,len(site_treatments)))
        count=0
        for i in site_treatments:
                site_T= site_df[(site_df['ID.Treatment'].values == [i])]
                SOC_dyn_T=site_T['SOC'] #(tC/ha)
                if(missing_years_to30>0): #if experiment has less than 30y, fill missing years with Nans
                        empty_ar = np.empty(missing_years_to30)
                        empty_ar[:]=np.NaN
                        SOC_dyn_T_30=np.append(SOC_dyn_T,empty_ar)
                else: #cut experiment to 30th year
                        SOC_dyn_T_30 = np.array(SOC_dyn_T[0:30])
                TREATMENTS[:,count]=SOC_dyn_T_30
                count+=1
        SITE_TREATMENTS.append(TREATMENTS)
                              
        #GET initial years for 4p1000
        site_T0= site_df[(site_df['ID.Treatment'].values == [site_T0_name])]
        SITE_year0[j] = np.min(site_T0['Year'])
        date_init = np.str(1980)
        date_end = np.str(2010)
        exper_len = np.int(date_end) - np.int(date_init) + 1
        SITE_exper_len[j] = exper_len
                              
        clay = np.mean(site_T0['Clay'])*100. #%
	silt = np.mean(site_T0['Silt'])*100. #%
	if(np.isnan(silt)):
		silt=(100.-clay)/2.
        bd=np.mean(site_T0['Bulk density'])*1000. #kgsoil/m3
        ph=np.mean(site_T0['pH'])

        SITE_BD[j]=bd
        SITE_pH[j]=ph
        SITE_clay[j]=clay
	SITE_silt[j]=silt
        SITE_date_init[j]=date_init
        SITE_date_end[j]=date_end

        #soil_temp_ss = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"temp_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
        #soil_hum_ss  = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"hum_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
        soil_temp    = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"temp_"+site+"_"+date_init+"_"+date_end+".txt"
        soil_hum     = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"hum_"+site+"_"+date_init+"_"+date_end+".txt"
                              
        #.....................

        #with open(soil_temp_ss) as fileID:
        #        # C
        #        temp_in = np.array(map(float,fileID))
        #        SITE_temp_in.append(temp_in)

        with open(soil_temp) as fileID:
                # C
                temp_t = np.array(map(float,fileID))
                SITE_temp_t.append(temp_t)

        #with open(soil_hum_ss) as fileID:
        #        #conversion kg(H2O)/m2(soil) to m3(H2O)/m3(soil)
        #        water_in = np.array(map(float,fileID))/100.
        #        SITE_water_in.append(water_in)

        with open(soil_hum) as fileID:
                #conversion kg(H2O)/m2(soil) to m3(H2O)/m3(soil)
                water_t = np.array(map(float,fileID))/100.
                SITE_water_t.append(water_t)

        #----------------------------------------------------------------------
        #  determine litter_inc for current site
        #----------------------------------------------------------------------
        # returns ABOVE, BELOW (gC/m2/day) and YEAR of measurements
        SAVE_FILE=-1
        ABOVE,BELOW,YEAR = AB_NanZeroRemover(site_T0,site_T0_name,SAVE_FILE,ROOTDIR)

        #---------------------------------------------------------------------------

        ABOVE_mean=np.mean(ABOVE)
        BELOW_mean=np.mean(BELOW)

        SITE_ABOVE_mean[j]=ABOVE_mean
        SITE_BELOW_mean[j]=BELOW_mean

        SITE_ERR2_ABOVE_mean[j]=np.std(ABOVE)**2/len(ABOVE)
        SITE_ERR2_BELOW_mean[j]=np.std(BELOW)**2/len(BELOW)
        if (SITE_ERR2_ABOVE_mean[j] == 0):
                SITE_ERR2_ABOVE_mean[j]=0.05 # to be checked
                SITE_ERR2_BELOW_mean[j]=0.05

	cov_AB_mean = np.cov(ABOVE,BELOW)/np.sqrt(len(ABOVE))
        SITE_cov_mean.append(cov_AB_mean)

        #Litter calculated as the sum of ABOVE and BELOW
        LITTER_mean = ABOVE_mean + BELOW_mean
        LITTER_var = SITE_ERR2_ABOVE_mean[j]+SITE_ERR2_BELOW_mean[j]

        SITE_litterinc[j]      = LITTER_mean
        SITE_litterinc_err2[j] = LITTER_var

        litter_inc = np.array([ABOVE_mean,BELOW_mean])    # litter C inputs parameters (gC/m2/day)
        Err2_litter_inc = np.array([SITE_ERR2_ABOVE_mean,SITE_ERR2_BELOW_mean])

        #===================================================
        # SOC,VARIANCE, YEARS with nan and zero removed
        #===================================================
        sc,var,yy0,ss0=SOC_NanZeroRemover(site_T0,site_T0_name)
        print sc
        SOC_clean_exp_array.append(sc)
        SOC_clean_exp_variance.append(var)
        SOC_clean_year.append(yy0)
        SOC_clean_where.append(ss0)

        SOC_all_meas = np.asarray(site_T0['SOC']).astype(np.float64)
        n_an = np.int(len(SOC_all_meas))
        SITE_n_an[j]=n_an
        SITE_SOC_data_init[j]=sc[0]
        print 'Initial SOC data: ',SITE_SOC_data_init[j]

        #===================================================

        j+=1

#>>>>>>>> END_OF_INITIALIZATION <<<<<<<<<<<<<<<<<<<<<<<<<<<<
tstart = time.time()
#>>>>>>>>>>>>>>>>>>>>>>>>

                              
######################################
#Set bounds and constraints for the optimization
#########
bnds=[(0,10)]
                              
#############################################
#PARAMETERS
#############################################
param_pi=0.66
param_pa=0.33
alpha_pl=2.5e12
rate_pa=0.02
rate_break=0.019
rate_leach=0.0015
kaff_des=1
param_p1=0.186
param_p2=0.216
kaff_lb=290.
alpha_lb=2.6e12
rate_bd=0.0036
rate_ma=0.02
cue_ref=0.6
cue_t=0.012
tae_ref=15.
matpot=15.
lambda_p=2.1e-4 
porosity = 0.6
kamin=0.2
param_pb=0.5
param_c1=0.86

########################                              
                              
litterin_sites = np.zeros((N_sites,4))
SOC_out_all = []
optim_4p1000_sites = np.zeros(N_sites)
SITE_initial_cond = np.zeros((N_sites,5))
for j in range(N_sites):
        Current_Site_Index=j
        site       = site_names_all[j]
        print '**********'
        print 'SITE',site
        print '**********'
        YEARS      = SOC_clean_year[j]
        SOC_data   = SOC_clean_exp_array[j]
        ss0 = SOC_clean_where[j]
                              
        n_an = np.int(SITE_n_an[j])
        n_days=np.int(n_an*365)
        n_an_4p1000 = 30
        n_days_4p1000 = n_an_4p1000*365
                           
        SOC_data_init = SITE_SOC_data_init[j] #first year SOC
        print 'initial SOC',SOC_data_init
        print 'check year number'
        SOC_var    = SOC_clean_exp_variance[j]
        clay       = SITE_clay[j]
	silt       = SITE_silt[j]
	print '###########'
	print 'SILT',silt
	print '###########'
	param_claysilt = clay+silt
        param_bulkd      = SITE_BD[j]
        pH      = SITE_pH[j]   
                              
        temp_t     = SITE_temp_t[j]
        water_t    = SITE_water_t[j]

        err_above = SITE_ERR2_ABOVE_mean[j]
        err_below = SITE_ERR2_BELOW_mean[j]
        ABOVE_mean = SITE_ABOVE_mean[j]
        BELOW_mean = SITE_BELOW_mean[j]
        ab_ratio = ABOVE_mean/BELOW_mean


        #LITTER INCOME AT SITE
        #above-below array to calculate uncertainties
        AB_BE_array=np.array([ABOVE_mean,BELOW_mean])
	cov_mean = SITE_cov_mean[j]

        #litter income at site (obs gC/m2/day)
        litter_mean = SITE_litterinc[j]
        Err2_litter_mean = SITE_litterinc_err2[j]

        #litter_inc in above and below
        litter_inc=np.array([ABOVE_mean,BELOW_mean])

        #litter income prior (1D)
        in_opt = litter_mean*(1+0.004)

        #Spinup 
        POM_ss=SITES_spinup.loc[j]['POM']
        LMWC_ss=SITES_spinup.loc[j]['LMWC']
        AGG_ss=SITES_spinup.loc[j]['AGG']
        MIC_ss=SITES_spinup.loc[j]['MIC']
        MAOM_ss=SITES_spinup.loc[j]['MAOM'] 
                              
        state=np.array([POM_ss,LMWC_ss,AGG_ss,MIC_ss,MAOM_ss])    
                              
        #Optimized parameters
        eact_lb = SITE_eact_lb.iloc[j]
	eact_pl = SITE_eact_pl.iloc[j]
        kaff_pl = SITE_kaff_pl.iloc[j]
        
                              
        #4p1000 analysis
        print '>>>>>>>>>>>>>>>>'
        print '4p1000 ',site
        target = np.sum(state)*0.004


        t_fwd=np.arange(0,n_days_4p1000-365) #Togli un anno perche dopo aggiungi spinup
        print ' '

	#fwd_output = odeint(decomp,state,t_fwd,args=(temp_t,water_t,litter_mean))
	#print fwd_output

        opt_mean=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,options={'disp':True})
        litter_opt = opt_mean.x
        print "SLSQP: Optimum solution:", litter_opt
        total_opt_in=np.sum(opt_mean.x)

        print 'Initial litter (gC/m2/day):'
        print litter_mean
        print '4p1000 litter (gC/m2/day):'
        print total_opt_in


        #optimized litter pools and total litter (save)
        in_opt_save = np.append(litter_opt*365./100.,total_opt_in*365./100.) #tC/ha/year
	#np.save(out2,in_opt_save)

        optim_4p1000_sites[j]=litter_opt

        #calculate percentage increase/decrease of inputs
        input_change=(total_opt_in-litter_mean)*100/litter_mean
        print "% change of litter inputs:",input_change

        #Check 4p1000 target is reached
        END=predict_c.shape[0]-1
        print 'END',END
        C_fin=np.sum(predict_c[END])
        C_init=np.sum(predict_c[0])
        SUMO=C_fin-C_init
        Target_reached=(C_fin-C_init)/(C_init*n_an_4p1000)
        if (Target_reached < 0.005 and Target_reached > 0.003):
                print "Target reached successfully"
                print "Target reached :", Target_reached
        else:
                print "Target not reached",Target_reached
                print 'C init', C_init
                print predict_c[0]
                print 'C_fin', C_fin
		
      

	############################################################
        #CREATE output file with PREDICTED CARBON over 30 years (standard and 4x100)
        #                       + EXPERIMENT SOC filled with Nan for 30 years (T0 + other treatments)
        ###########################################################	

	predict_c_standard_pools_day = odeint(decomp,state,t_fwd,args=(temp_t,water_t,litter_mean))
	predict_c_standard_pools=DailyToYearly(predict_c_standard_pools_day,1)
	predict_c_standard = np.sum(predict_c_standard_pools,axis=1)

	SOC_model_standard_pools = np.concatenate(([state],predict_c_standard_pools))
	SITE_SOC_model=np.sum(state)
        SOC_model_standard = np.concatenate((np.array([SITE_SOC_model]),predict_c_standard))
	
        #4p1000
        predict_c_opt_pools_day=odeint(decomp,state,t_fwd,args=(temp_t,water_t,litter_opt))
	predict_c_opt_pools=DailyToYearly(predict_c_opt_pools_day,1)
        predict_c_opt = np.sum(predict_c_opt_pools,axis=1)
        
	SOC_model_opt_pools =np.concatenate(([state],predict_c_opt_pools))
        SOC_model_opt = np.concatenate((np.array([SITE_SOC_model]),predict_c_opt))
	
	
	year_out = np.arange(1,31)

	SOC_pools_out=np.stack((SOC_model_standard_pools/100.,SOC_model_opt_pools/100.))#tC/ha
	print len(year_out)
	print len(SOC_model_standard/100.)
	print len(SOC_model_opt/100.)
	SOC_out=np.stack((year_out,SOC_model_standard/100.,SOC_model_opt/100.)) #tC/ha
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

                cov_AB_mean=cov_mean #covariance between ABOVE mean and BELOW mean if no Nans
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
                                print 'total litter in (gC/m2/day):',in_rand_param
                                print '****************'
                                #Minimize J_new for the generated array
                                opt=minimize(J_new, in_rand_param, method='SLSQP',bounds=bnds, options={'disp':True})
                                opt_param=opt.x
                                opt_parameters_MC[sample_shape]=opt_param

                                print "OPTIMUM LITTER IN (gC/m2/day)"
                                print opt_param
                                print '****************'
                                print "Litter in increase (%)"
                                print (opt_param-in_rand_param)*100./in_rand_param
                                print '****************'
                                sample_shape+=1
		
		out_priors_and_opt = np.stack((in_rand_param_MC*365./100.,opt_parameters_MC*365./100.)) #save as tC/ha/yr
		np.save(out_priors,out_priors_and_opt)

		#STANDARD ERROR CALCULATION for the optimized litter inputs to reach 4x1000
		litter_opt_err=np.std(opt_parameters_MC[:])/np.sqrt(MC_length)
		
		#Error litter = SE per litter in e in opt
		save_lit = np.stack((litter_mean*365./100.,Err2_litter_mean*365./100.,litter_opt[0]*365./100.,litter_opt_err*365./100.))#save as tC/ha/yr
                np.save(out_lit,save_lit)      

                #in MILLENNIAL save_lit = litter_tot_save
                litter_tot_save = save_lit
                np.save(out_lit_tot,litter_tot_save)


#Print optimized parameters                              
for i in range(N_sites):
        print site_names_all[i],' C inputs to 4p1000: ', optim_4p1000_sites[i]                              

out_lit_tot.close()
out_lit.close()
out_priors.close()
out_mo_pools.close()
out_mo.close()
#out2.close()

tend = time.time()
tot_time=tend-tstart
print " "
print " Ok, done. Total time: ",tot_time                              


