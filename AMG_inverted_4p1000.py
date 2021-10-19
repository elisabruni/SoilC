############################################################
#AMGv1 and AMGv2
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


#############
#Type of simulation
AMGv2=False

optimization_param=False
optimization_4p1000=True

###################
#Data

ROOTDIR='/Users/ebruni/Desktop/DOTTORATO/'
loc_exp = ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/SITES_dataset_multimodel.xlsx'
OUTPUT_files=ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OUTPUTS_4p1000_v5/AMG/'

############
#Save data
if(optimization_param):
        out1=open("SOC_AMG_optim.txt","wb")
if(optimization_4p1000):
        out_mo_pools=open(OUTPUT_files+"SOC_AMG_pools.txt","wb")
        out_mo=open(OUTPUT_files+"SOC_AMG.txt","wb")
	out_lit=open(OUTPUT_files+"Litter_income_AMG.txt","wb")
        out_priors = open(OUTPUT_files+"priors_and_opt_in_AMG.txt","wb")
	out_lit_tot = open(OUTPUT_files+"Litter_tot_AMG.txt","wb")

#########################################################
def DailyToYearly(filein,typeConv):
#########################################################
                c=0
                x=0
		#Set number of years
                var_year=np.zeros(np.int(len(filein)/365))
                lent = np.int(len(var_year))
                for j in range(0,lent):
                        if(typeConv==0.): #for precipitation and PET (cumulative)
                                var_year[c]=np.sum(filein[x:365+x])
                        else: #for Temperature (mean)
                                var_year[c]=np.sum(filein[x:365+x])/365
                        x=x+365
                        c+=1
                return var_year

#########################################################
def AB_NanZeroRemover(site_T0,site_T0_name):
#########################################################
        #######################################################################
        # remove NaNs and ZEROs from ABOVE and BELOW
        #######################################################################
        yy=np.asarray(site_T0['Year']).astype(np.int16)             # array of years
        aa=np.asarray(site_T0['ABOVE']).astype(np.float64)  # array of C_above (tC/ha/yr)
        bb=np.asarray(site_T0['BELOW']).astype(np.float64)  # array of C_below (tC/ha/yr)
        aa0=np.where(np.isnan(aa),999,aa)  # replace NaN with 999 in "above"
        YEAR=yy[aa0<999]                     # select years where above>0
        abo=aa[aa0<999]                      # select ABOVE>0
        bel=bb[aa0<999]                      # select corresponding BELOW
        print site_T0_name,': AB_NanZeroRemover --> selected ',len(YEAR),' out of ',len(yy)
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
        ########################################
        #No need to convert sc and vv0 in rothc
        #######################################
        YEAR=yy0[ss0>0]                          # select the years corresponding to sc
        std2=np.std(sc)**2                      # square standard deviation of the measurements (use when no error provided - ??)
        if (std2 == 0):
                std2=(BigRelativeError*sc)**2    # <-- check  # if std == 0 use BigRelativeError in spite of std
        vv0=np.where(np.isnan(vv0),std2,vv0)   # Replace NaN in variance array with std2
        vv0=np.where(vv0==0,std2,vv0)           # Replace 0 in variance with std2
        var=vv0[ss0>0]                           # Restrict variance corresponding to the selected SOCs data
        print site_T0_name,': SOC_NanZeroRemover (cleanup of SOC data) --> selected ',len(YEAR), ' years out of ',len(yy0)
        return sc,var,YEAR,ss0

################################
################################
#ENVIRONMENTAL FUNCTIONS
################################
################################
#############################################
def control_temp_func(temp_mean):
#############################################
	#temp_mean = mean annual temperature

	aT = 25. 
	cT = 0.120 #1/Kelvin 
	TRef = 15. #Celcius degrees

	bT = (aT - 1)*np.exp(cT*TRef)

	if(temp_mean>0):
		temp_func = aT/(1+bT*np.exp(-cT*temp_mean))
	else:
		temp_func = 0.

 	return temp_func


		
#############################################
def control_moist_func(water_mean,pot_evapot):
#############################################
	#water_mean = cumulative annual water inputs (precipitation and irrigation water)
	#pot_evapot = potential evapotranspiration (mm)

	aH = 3.0/100
	bH = 5.247 #1/m

	moist_func = 1/(1+ aH*np.exp(-bH*(water_mean-pot_evapot)/1000))

	return moist_func


#############################################
def control_clay_func(clay):
#############################################
	#clay (gclay/kgsoil)
	global AMGv2
	if(AMGv2):
		aM =2.519/1000 #AMGv2
	else:
		aM = 2.720/1000 #g/kg2
	
	clay_func = np.exp(-aM*clay)

	return clay_func

#############################################
def control_carbonate_func(carbonate):
#############################################	
	#carbonate (CaCO3) content (gCaCo3/kgsoil)
	global AMGv2
	if(AMGv2):
		cM = 1.50/1000 #AMGv2
	else:
		cM = 1.67/1000  #g/kg
	
	carbonate_func = 1/(1+cM*carbonate)
	return carbonate_func

#In AMGv2
#############################################
def control_pH_func(pH):
#############################################
	apH = 0.112 
	bpH = 8.5
	pH_func = np.exp(-apH*(pH-bpH)**2)
	return pH_func

#In AMGv2
#############################################
def control_CN_ratio_func(CN_ratio):
#############################################
	aCN = 0.060 
	bCN = 11.
	CN_ratio_func = 0.8*np.exp(-aCN*(CN_ratio-bCN)**2) + 0.2
	return CN_ratio_func

#############################################
#Mineralization rate function
#############################################
#k0 is the potential mineralization rate set at 0.165 and 0.290 for AMGv1 and AMGv2
#To be optimized
#############################################
if(AMGv2):
	def k_rate_func(k0, temp_mean, water_mean, pot_evapot, clay, carbonate, pH, CN_ratio):
		k_rate = k0*control_temp_func(temp_mean)*control_moist_func(water_mean,pot_evapot)*control_clay_func(clay)*control_carbonate_func(carbonate)*control_pH_func(pH)*control_CN_ratio_func(CN_ratio)
		return k_rate
else:
	def k_rate_func(k0, temp_mean, water_mean, pot_evapot, clay, carbonate):
		k_rate = k0*control_temp_func(temp_mean)*control_moist_func(water_mean,pot_evapot)*control_clay_func(clay)*control_carbonate_func(carbonate)
		return k_rate

#############################################
def humcoef_func(CN_crop):
#############################################
	ecoef = 0.69
	fcoef = 11.2
	humcoef_AG = 1- (ecoef*CN_crop)/(fcoef+CN_crop)
	humcoef_BG = 0.39
	#add other humcoef for OM
	humcoef = np.array([humcoef_AG,humcoef_BG])
	return humcoef

#############################################
#If humcoefficients already defined
def forward_ACT(n_an,init,litterin,temp_mean,water_mean,pot_evapot,k0,history_land,humcoef):
#############################################
	#litterin 2x array (AG,BG)
	#humcoef 2x array (AG,BG) humus coefficient
	global AMGv2

	#number of years = n_an because you don't have initialization
	ACT_cpools_tyear = np.zeros(n_an)
	#Set first year of active pool
	ACT_cpools_tyear[0]=init*(1-history_land)
	dt=1
	for x in range(1,n_an):
		temp_t = temp_mean[x]
		water_t = water_mean[x]
		pot_evapot_t = pot_evapot[x]
		ACT_current= ACT_cpools_tyear[x-1]
		if(AMGv2):
			kk_ma = k_rate_func(k0,temp_t, water_t, pot_evapot_t, clay, carbonate, pH, CN_ratio)
		else:
			kk_ma = k_rate_func(k0,temp_t, water_t, pot_evapot_t, clay, carbonate)

		ACT_next = ACT_current + (np.sum(litterin*humcoef) - kk_ma*ACT_current)*dt
		ACT_cpools_tyear[x]=ACT_next

	return ACT_cpools_tyear	

#############################################
def SOC_stable_func(init,history_land):
#############################################

	SOC_stable = init*history_land
	
	return SOC_stable

#############################################
def SOC_tot_func(SOC_stable, SOC_act):
#############################################
	SOC_tot = SOC_stable + SOC_act
	return SOC_tot

############################################################
# param_opt ######## OBJECTIVE FUNCTION
############################################################
def Jparam_new(k0_new):
    global param_predict_c_forward, SOC_data, SOC_var, ss0
    global NEW_ITER, clay, SOC_data_init,litter_inc, temp_in,water_in, PET_in,c_vegetated,t_param, n_an
    global clay,temp_t,water_t,PET_t,ssp_select_optim

    Delta=np.zeros(len(YEARS)-1)
    j=Current_Site_Index

    #Stable SOC
    stable_SOC = SOC_stable_func(SOC_data_init,arable_scoef)
    #Active SOC
    param_predict_c_forward =forward_ACT(n_an,SOC_data_init,litter_inc,temp_t,water_t,PET_t,k0_new,arable_scoef,humcoef)

    #Total SOC
    tot_SOC  = SOC_tot_func(stable_SOC,param_predict_c_forward)


    ssp_select_optim=tot_SOC[ss0>0]

    #print 'SOC_data',SOC_data
    #print 'model',ssp_select_optim

    #per l'ottimizzazione togli il primo anno, visto che e preso direttam dai dati
    Delta[0:]=ssp_select_optim[1:]-SOC_data[1:]
    Delta2w=Delta**2/SOC_var[1:]
    Jparam_new = np.sum(Delta2w)

    if ( NEW_ITER % CHI2_PRINT_FREQUENCY == 0 ):
        print "==> NEW_ITER ",NEW_ITER," param_new=",k0_new," Jparam_new=",Jparam_new
        print "param_predict_c_forward",param_predict_c_forward
        print "SOC_data",SOC_data
        print "Delta",Delta
        print "Delta2w",Delta2w
        print "Error %",np.sqrt(SOC_var)/SOC_data

    NEW_ITER+=1
    return Jparam_new


############################################################
# J_new ######## OBJECTIVE FUNCTION
############################################################
def J_new(in_new):
        global NEW_ITER
        global n_an_4p1000,n_an,optimized_k
        global SOC_data_init, predict_c, SOC_stable
        global clay,temp_t,water_t,PET_t,arable_scoef
	global humcoef_OM,litter_inc

	#Calculate proportion of straw / OM
	#in_new_tot = np.sum(in_new)
	#litter_inc_tot = np.sum(litter_inc)
	prop_straw = litter_inc/in_new
	humcoef_pred = humcoef*prop_straw+humcoef_OM*(1.-prop_straw)
	humcoef_pred = np.maximum(humcoef_pred,humcoef)
	print humcoef_pred


	#Active
	predict_c_fwd_4p1000=forward_ACT(n_an_4p1000,SOC_data_init,in_new,temp_t,water_t,PET_t,optimized_k,arable_scoef,humcoef_pred)
	#Stable
	stable_SOC_4p1000 = SOC_stable_func(SOC_data_init,arable_scoef)
	#Total
	predict_c = predict_c_fwd_4p1000+stable_SOC_4p1000
	
        J_new = abs(target*n_an_4p1000 - np.sum(predict_c[predict_c.shape[0]-1]-predict_c[0]))

        #print 'OBJ FUN', J_new
        #print 'predict_c',predict_c

        NEW_ITER+=1
        if ( NEW_ITER % 100 == 0 ):
                print "NEW_ITER ",NEW_ITER," in_new=",in_new," J_new=",J_new
        param_opt=in_new
        return J_new

############################################################
#Initialization
NEW_ITER = 0
CHI2_PRINT_FREQUENCY=2

#############################################
#PARAMETERS
#############################################
#Fraction of initial SOC in the stable part
#############################################
arable_scoef=0.65
grassland_scoef=0.4


k0_prior=0.165 #need to be optimized



######################################################
#For each site: Set SITE name and experiment duration
#####################################################

C_input_exp = pd.read_excel(loc_exp)
site_names_all = C_input_exp['ID.Site'].unique()[2:len(C_input_exp)]
site_names_all = map(str, site_names_all)

N_sites_all=len(site_names_all)

#Control plots
site_T0_array_all=np.array(['CHNO3_Min', 'COL_T0', 'CREC3_Min', 'FEU_T0', 'JEU2_M0', 'LAJA2_Min', 'RHEU1_Min', 'RHEU2_T0','ARAZ_D0_N0', 'ULT_P0_B', 'BROAD_3_Nill', 'TREV1_Min','AVRI_T1TR','BOLO_T0', 'GRAB_CP','MUNCHE_CP','RITZ_CP'])


#Stationary solution array for each experiment T0
N_sites=len(site_names_all)

print N_sites

if(optimization_4p1000):
	##################
	#Import optimized k
	##################
	n=0 # contatore
	inp=open(ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM5/optim_param_AMG.txt','rb')
	while inp.read(1):
		inp.seek(-1,1)
		XX_optimized_k=np.load(inp)
		n+=1

SOC_exp_array=[] #interpolated SOC dynamics experiments
SOC_clean_exp_array=[]
SOC_clean_exp_variance=[]
SOC_clean_year=[]
SOC_clean_where=[]

SITE_SOC_data_init=np.zeros(N_sites)
SITE_year0=np.zeros(N_sites)
SITE_date_init=np.zeros(N_sites)
SITE_date_end=np.zeros(N_sites)
SITE_date_init_ss=np.zeros(N_sites)
SITE_date_end_ss=np.zeros(N_sites)
SITE_exper_len = np.zeros(N_sites)
SITE_clay=np.zeros(N_sites)
SITE_carbonate=np.zeros(N_sites)
SITE_litterinc = np.zeros(N_sites)
SITE_litterinc_err2 = np.zeros(N_sites) #sum
SITE_Err2_litter_arr= np.zeros((N_sites,2)) #array
SITE_water_in=[]
SITE_temp_in=[]
SITE_PET_in=[]
SITE_water_t=[]
SITE_temp_t=[]
SITE_PET_t=[]
SITE_ABOVE_mean=np.zeros(N_sites)
SITE_BELOW_mean=np.zeros(N_sites)
SITE_ERR2_ABOVE_mean=np.zeros(N_sites)
SITE_ERR2_BELOW_mean=np.zeros(N_sites)
SITE_ABOVE=[]
SITE_BELOW=[]
SITE_cov_mean = []
SITE_TREATMENTS=[]
SITE_humcoef=np.zeros((N_sites,2))
SITE_humcoef_OM = np.zeros(N_sites)
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
	#GET initial years for ss and forward
        site_T0= site_df[(site_df['ID.Treatment'].values == [site_T0_name])]
        SITE_year0[j] = np.min(site_T0['Year'])
        if(optimization_4p1000):
                date_init = np.str(1980)
                date_end = np.str(2010)
        else:
                date_init = np.str(np.int(year_0))
                date_end = np.str(np.int(year_end))

	date_init_ss = np.str(np.int(year_0 - 30))
        date_end_ss = np.str(np.int(year_0 - 1))
        exper_len = np.int(date_end) - np.int(date_init) + 1
        SITE_exper_len[j] = exper_len
        clay = np.mean(site_T0['Clay'])
	carbonate = np.mean(site_T0['CaCO3'])
	humcoefAG=np.mean(site_T0['humcoefAG'])
	humcoefBG=np.mean(site_T0['humcoefBG'])
	humcoef = np.array([humcoefAG,humcoefBG])
	
	humcoefOM=np.mean(site_T0['humcoefOM'])
	
	SITE_humcoef[j,:]= humcoef
	SITE_humcoef_OM[j]=humcoefOM
	SITE_carbonate[j]=carbonate
        SITE_clay[j]=clay
        SITE_date_init[j]=date_init
        SITE_date_end[j]=date_end
        SITE_date_init_ss[j]=date_init_ss
        SITE_date_end_ss[j]=date_end_ss

        soil_temp_ss = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"temp_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
        soil_hum_ss  = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"water_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
        soil_PET_ss  = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"PETcorr_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
        soil_temp    = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"temp_"+site+"_"+date_init+"_"+date_end+".txt"
        soil_hum     = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"water_"+site+"_"+date_init+"_"+date_end+".txt"
        soil_PET     = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"PETcorr_"+site+"_"+date_init+"_"+date_end+".txt"
        #.....................

	with open(soil_temp_ss) as fileID:
                # C to K
                temp_in = np.array(map(float,fileID))
                temp_in_monthly=DailyToYearly(temp_in,1)
                SITE_temp_in.append(temp_in_monthly)

        with open(soil_temp) as fileID:
                # C to K
                temp_t = np.array(map(float,fileID))
                temp_t_monthly=DailyToYearly(temp_t,1)
                SITE_temp_t.append(temp_t_monthly)

        with open(soil_hum_ss) as fileID:
                #rain+snow(mm/day)
                water_in = np.array(map(float,fileID))
                water_in_monthly=DailyToYearly(water_in,0)
                SITE_water_in.append(water_in_monthly)

        with open(soil_hum) as fileID:
                #rain+snow(mm/day)
                water_t = np.array(map(float,fileID))
                water_t_monthly=DailyToYearly(water_t,0)
                SITE_water_t.append(water_t_monthly)

        with open(soil_PET_ss) as fileID:
                #potential evapot Pennman-Monteith (mm/s)
                PET_in = np.array(map(float,fileID))
                PET_in_monthly=DailyToYearly(PET_in,0)*3600.*24.
                SITE_PET_in.append(PET_in_monthly)

        with open(soil_PET) as fileID:
                #potential evapot Pennman-Monteith (mm/s)
                PET_t = np.array(map(float,fileID))
		PET_t_monthly=DailyToYearly(PET_t,0)*3600.*24.
                SITE_PET_t.append(PET_t_monthly)
        #.....................

	#----------------------------------------------------------------------
        #  determine litter_inc for current site
        #----------------------------------------------------------------------
        # returns ABOVE, BELOW (tC/ha/month) and YEAR of measurements
        ABOVE,BELOW,YEAR = AB_NanZeroRemover(site_T0,site_T0_name)
        #---------------------------------------------------------------------------

	SITE_ABOVE.append(ABOVE)
	SITE_BELOW.append(BELOW)

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

	litter_inc = np.array([ABOVE_mean,BELOW_mean])    # litter C inputs parameters (tC/ha/yr)
        Err2_litter_inc_arr = np.array([SITE_ERR2_ABOVE_mean[j],SITE_ERR2_BELOW_mean[j]])
	SITE_Err2_litter_arr[j]=Err2_litter_inc_arr

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


################################################
# OPTIMIZATION CONSTRAINTS
#ordine AG,BG
################################################
def constr1(x):
        global a_constr,ab_ratio
        return x[0]-a_constr*x[1]*ab_ratio

def constr2(x):
        global b_constr,ab_ratio
        return b_constr*x[1]*ab_ratio-x[0]


tstart = time.time()

######################################
#Set bounds and constraints for the optimization
#########
bnds=[(0,10),(0,10)]

a_constr=1
b_constr=1
ab_ratio_var=(1-a_constr)*100
con1={'type':'ineq','fun':constr1}
con2={'type':'ineq','fun':constr2}
cons=[con1,con2]

#MAIN

litterin_sites = np.zeros((N_sites,4))
SOC_out_all = []
optim_param_sites =np.zeros(N_sites)
optim_4p1000_sites = np.zeros((N_sites,2))
for j in range(N_sites):
        Current_Site_Index=j
        site       = site_names_all[j]
        YEARS      = SOC_clean_year[j]
        SOC_data   = SOC_clean_exp_array[j]
        ss0 = SOC_clean_where[j]
        n_an = np.int(SITE_n_an[j])
        n_an_4p1000 = 30
        n_an_ss = 30
        SOC_data_init = SITE_SOC_data_init[j] #first year SOC
	print 'initial SOC',SOC_data_init
        print 'check year number'
        SOC_var    = SOC_clean_exp_variance[j]
        clay       = SITE_clay[j]
	carbonate  = SITE_carbonate[j]
	humcoef = SITE_humcoef[j,:]
	humcoef_OM = SITE_humcoef_OM[j]
        temp_in    = SITE_temp_in[j]
        #print 'TEMP IN', temp_in
        water_in   = SITE_water_in[j]
        #print 'WATER IN', water_in
        PET_in     = SITE_PET_in[j]
        temp_t     = SITE_temp_t[j]
        water_t    = SITE_water_t[j]
        PET_t     = SITE_PET_t[j]
	ABOVE = SITE_ABOVE[j]
	BELOW = SITE_BELOW[j]
        err2_above = SITE_ERR2_ABOVE_mean[j]
        err2_below = SITE_ERR2_BELOW_mean[j]
        ABOVE_mean = SITE_ABOVE_mean[j]
        BELOW_mean = SITE_BELOW_mean[j]
	#Approx CHECK if better
        litter_in_err = np.sqrt(err2_above+err2_below)

        ab_ratio = ABOVE_mean/BELOW_mean
	if(optimization_4p1000):
		optimized_k = XX_optimized_k[j]

        #LITTER INCOME AT SITE
        #above-below array to calculate uncertainties
        AB_BE_array=np.array([ABOVE_mean,BELOW_mean])

        #litter income at site (obs tC/ha/yr)
        litter_mean = SITE_litterinc[j]
        #litter_inc in above and below
        litter_inc=np.array([ABOVE_mean,BELOW_mean])

	#Error
        Err2_litter_mean = SITE_litterinc_err2[j] #sum
	Err2_litter_inc = SITE_Err2_litter_arr[j] #array

	#to be saved
        litter_inc_save=np.append(litter_inc,litter_mean)
        litter_inc_err_save=np.append(np.sqrt(Err2_litter_inc),np.sqrt(Err2_litter_mean))
	
	#print 'Err2_litter_inc'
	#print Err2_litter_inc

        #litter income prior (2D)
        in_opt = litter_inc*(1+0.004)

        print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        print '>>>>>>>>>>>>>> Analysis for SITE ',site
        print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        #--------------------------------------------------------------#

	#Calculate SOC stocks dynamics
	#SOC first year (equivalent to spinup) -> SOC_data_init
	print 'Stationary solution before 4x1000: '
	print SOC_data_init
        if(not optimization_4p1000):
		#Standard simulation (k0 non optim)
        	SOC_stable = SOC_stable_func(SOC_data_init,arable_scoef)
  		SOC_act = forward_ACT(n_an,SOC_data_init,litter_inc,temp_t,water_t,PET_t,k0_prior,arable_scoef,humcoef)
         	SOC_tot = SOC_tot_func(SOC_stable,SOC_act)
		print 'Site SOC dynamics before 4x1000: '
		print SOC_tot
		#print 'Active'
		#print SOC_act
		#print 'Stable'
		#print SOC_stable
	
		ssp_select=SOC_tot[ss0>0]
		print YEARS	
		print 'Modelled SOC (tC/ha) where data avail'
		print ssp_select
		print 'Observed (tC/ha)'
		print SOC_data
		print 'sso'
		print ss0
	
        #--------------------------------------------------------------#

        print 'LITTER:',litter_inc
        print 'CLAY:',clay
	print 'CACO3:',carbonate
        print 'TEMP:',np.mean(temp_t)
        print 'WATER:',np.mean(water_t)

	#OPTIMIZATION PARAMETERS
        if(optimization_param):
                print 'OPTIMIZATION PARAM'
                bnds_param=bnds=[(0.03,1.)]
                opt_param_mean=minimize(Jparam_new, k0_prior, method='SLSQP', bounds=bnds_param, options={'disp':True})
                print "SLSQP: Optimum solution param_temp:",opt_param_mean.x
                param_min=opt_param_mean.x
                CHI2=Jparam_new(opt_param_mean.x)
                print 'SLSQP: CHI2 ',CHI2

		SOC_model_std = ssp_select	
                SOC_model_optim = ssp_select_optim
                SOC_error=SOC_var
                print 'SOC_error',SOC_error
                XX=np.stack((YEARS,SOC_model_std,SOC_model_optim,SOC_data,SOC_error))
                np.save(out1,XX)

                print "SOC model optim", SOC_model_optim
                print "SOC EXPER", SOC_data
                print "Optimized SOC"
                print " YEAR    SOCmodel     SOCdata"
                for k in range(len(SOC_model_optim)):
                        print '{0:6d} {1:7.1f} {2:7.1f}'.format(int(YEARS[k]), SOC_model_optim[k],SOC_data[k])
                optim_param_sites[j]=param_min


        if(optimization_4p1000):
		tstart = time.time()
                print '>>>>>>>>>>>>>>>>'
                print '4p1000 ',site
                target = SOC_data_init*0.004
                n_an_4p1000 = 30

                #SITE_SOC_dyn = np.concatenate((np.array([SOC_data_init]),fwd))
                #print 'SOC dynamics before opt', SITE_SOC_dyn

                print ' '
                opt_mean=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,constraints=cons,options={'disp':True})
                litter_opt = opt_mean.x
                print "SLSQP: Optimum solution:", litter_opt
                total_opt_in=np.sum(opt_mean.x)
		
		print 'Initial litter (tC/ha/yr):'	
		print litter_mean
		print '4p1000 litter (tC/ha/yr):'
		print total_opt_in


                #optimized litter pools and total litter (save)
                litter_opt_save = np.append(litter_opt,total_opt_in)
	
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
                        print "Target not reached"
                        print 'C init', C_init
                        print predict_c[0]
                        print 'C_fin', C_fin
		

                ############################################################
                #CREATE output file with PREDICTED CARBON over 30 years (standard and 4x100)
                #                       + EXPERIMENT SOC filled with Nan for 30 years (T0 + other treatments)
                ###########################################################
                #WORKING : metti condition su EOM, per scegliere humcoef

                predict_c_standard_stable=SOC_stable_func(SOC_data_init,arable_scoef)
                predict_c_standard_active=forward_ACT(n_an_4p1000,SOC_data_init,litter_inc,temp_t,water_t,PET_t,optimized_k,arable_scoef,humcoef)
		predict_c_standard_pools=np.stack((predict_c_standard_active,np.repeat(predict_c_standard_stable,n_an_4p1000)),axis=1)
		predict_c_standard =SOC_tot_func(predict_c_standard_stable,predict_c_standard_active)

                #4p1000
                predict_c_opt_stable=SOC_stable_func(SOC_data_init,arable_scoef)
                predict_c_opt_active=forward_ACT(n_an_4p1000,SOC_data_init,litter_opt,temp_t,water_t,PET_t,optimized_k,arable_scoef,humcoef)
                predict_c_opt_pools = np.stack((predict_c_opt_active,np.repeat(predict_c_opt_stable,n_an_4p1000)),axis=1)
                predict_c_opt=SOC_tot_func(predict_c_opt_stable,predict_c_opt_active)

                year_out = np.arange(1,31)

                SOC_pools_out=np.stack((predict_c_standard_pools,predict_c_opt_pools))#tC/ha
                SOC_out=np.stack((year_out,predict_c_standard,predict_c_opt)) #tC/ha

                np.save(out_mo_pools,SOC_pools_out)
                np.save(out_mo,SOC_out)

		############################
                #UNCERTAINTIES
                ###########################
                Uncert_Q = True
                if(Uncert_Q):

                        MC_length=2 #set the number of Monte Carlo simulations

                        #optimize for n variations of in_opt generatated randomly around the above/below covariance
                        opt_parameters_MC=np.zeros((MC_length,len(litter_inc)))

                        #prior above_below
                        ab_be_init_est=AB_BE_array*(1+0.004)

                        cov_AB_mean=np.cov(ABOVE,BELOW)/np.sqrt(len(ABOVE)) #covariance between ABOVE mean and BELOW mean if no Nans
                        #print 'cov',cov_AB_mean

                        if (np.all(cov_AB_mean)==0): #if covariance is 0, take the mean amongst all sites' covariances to generate cov_AB_mean
                                cov_AB_mean=np.mean(SITE_cov_mean)
                        #print cov_AB_mean

                        in_rand_param_MC = np.zeros((MC_length,len(litter_inc)))
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
                                        #in_rand_param=np.array(in_rand)  #array to optimize
                                        #Save all priors in an array
                                        in_rand_param_MC[sample_shape]=in_rand #add new generated sample to array on rand_in samples
                                        print '****************'
                                        print "LITTER IN generated randomly from ab_be_init_est:", in_rand
                                        print 'litter in (tC/ha/yr):',in_rand
					print 'Total litter in (tC/ha/yr):', np.sum(in_rand)
                                        print '****************'
                                        #Minimize J_new for the generated array
                                        opt=minimize(J_new, in_rand, method='SLSQP',bounds=bnds, options={'disp':True})
                                        opt_param=opt.x
					opt_param_tot = np.sum(opt_param)
                                        opt_parameters_MC[sample_shape]=opt_param

                                        print "OPTIMUM LITTER IN (tC/ha/yr)"
                                        print opt_param_tot
                                        print '****************'
                                        print "Litter in increase (%)"
                                        print (opt_param_tot-np.sum(in_rand))*100./np.sum(in_rand)
                                        print '****************'
					sample_shape+=1


                        #matrix of the optimum parameters
                        #print "opt parameters:", opt_parameters_MC

                        #save priors (litterin generated random) and out (optimized)

			print 'in_rand_param_MC,opt_parameters_MC'
			print in_rand_param_MC,opt_parameters_MC
			
			
                        out_priors_and_opt = np.stack((in_rand_param_MC,opt_parameters_MC)) #save as tC/ha/yr
                        np.save(out_priors,out_priors_and_opt)

                        #STANDARD ERROR CALCULATION for the optimized litter inputs to reach 4x1000

                        litter_opt_AB=np.std(opt_parameters_MC[:,0])/np.sqrt(MC_length)
			litter_opt_BE=np.std(opt_parameters_MC[:,1])/np.sqrt(MC_length)
			litter_opt_err = np.array([litter_opt_AB,litter_opt_BE])
			litter_opt_err_sum = np.sum(litter_opt_err)

			litter_opt_err_save = np.append(litter_opt_err,litter_opt_err_sum) #tC/ha/yr
                        #print "Uncertainties:",litter_opt_err


                        #Error litter = SE per litter in e in opt
                        save_lit = np.stack((litter_inc_save,litter_inc_err_save,litter_opt_save,litter_opt_err_save))#save as tC/ha/yr
                        np.save(out_lit,save_lit)
		
			litter_tot_save = np.stack((litter_mean, np.sqrt(Err2_litter_mean),total_opt_in,litter_opt_err_sum)) #tC/ha/yr
			np.save(out_lit_tot,litter_tot_save)


if(optimization_param):
        for i in range(N_sites):
                print site_names_all[i],' optimized param: ', optim_param_sites[i]
        out3=open("optim_param_AMG.txt","wb")
        np.save(out3,optim_param_sites)
        out1.close()
        out3.close()

	out1.close()
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
