#RothC-26.3 (https://www.rothamsted.ac.uk/sites/default/files/RothC_guide_WIN.pdf)
#Discrete to continuous form from Parshotam (1996) (https://doi.org/10.1016/0304-3800(95)00065-8)
#Matricial form
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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random

#################
#Type of simulation
##################
optimization_4p1000=True
optimization_param=False
#######################

#Save data

if(optimization_param):
        out1=open("SOC_RothC_optim_7.txt","wb")
if(optimization_4p1000):
        out2=open("SOC_RothC_4p1000_7.txt","wb")

plot_CM=False
if(plot_CM):
	pp0 = PdfPages("control_moisture_sites_test.pdf")

plot_CT=False
if(plot_CT):
        pp1 = PdfPages("control_temp_sites_test.pdf")

##################
#Functions
##################
##############################
#Convert daily data to monthly
##############################
def DailyToMonthly(filein,typeConv):
		c=0
		x=0
		var_month=np.zeros(np.int(len(filein)/365.*12.))
		month_days = [31,28,31,30,31,30,31,31,30,31,30,31]
		lent = np.int(len(var_month))
		for j in range(0,lent):
			for i in month_days:
				if(typeConv==0.): #for precipitation and PET (cumulative)
					var_month[c]=np.sum(filein[x:x+i])
				else: #for Temperature (mean)
					var_month[c]=np.sum(filein[x:x+i])/np.float(i)
				x=x+i
				c+=1
			if c>len(var_month)-12:
				break
		return var_month
		

#########################################################
def AB_NanZeroRemover(site_T0,site_T0_name,iout,ROOTDIR):
#########################################################
        if ( iout > 0 ):
                out1=open(ROOTDIR+"AB.data","wb")
        #######################################################################
        # remove NaNs and ZEROs from ABOVE and BELOW
        #######################################################################
        yy=np.asarray(site_T0['Year']).astype(np.int16)             # array of years
        aa=np.asarray(site_T0['ABOVE']).astype(np.float64)/12.  # array of C_above (tC/ha/month)
        bb=np.asarray(site_T0['BELOW']).astype(np.float64)/12.  # array of C_below (tC/ha/month)
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
def control_temp_func(temp_mean,param_temp):
#############################################
	#temp mean = monthly average temperature (C)

	#a = 47.91/(1.+np.exp(106.06/(temp_mean+18.27)))
	a = 47.91/(1.+np.exp(106.06/(temp_mean+param_temp)))
	#print 'contro temp',a
	return a

#############################################
def control_moist_func(clay, rain, mean_pot_evapot,c_param):
#############################################
	#clay %, rain (mm), mean_pot_evapot (mm)

	global n_an,n_month,accTSMD, maxTSMD

	#Maximum Topsoil Moisture Deficit (TSMD)

	if(c_param <1.):	#For vegetated soil
		#For a soil depth 0-23cm
		#If soil depth is for ex  X(cm) -> maxTSMD=maxTSMD*X/23
		maxTSMD = -(20.+1.3*clay*100-0.01*(clay*100)**2) 

	else:			#For bare soil divide by 1.8
		maxTSMD = -(20.+1.3*clay*100-0.01*(clay*100)**2)/1.8

	#Accumulated TSMD

	#Attenzione formula originale: invece di mean_pot_evapot c'e (open_pan_evapot*0.75). Siccome uso mean_pot_evapot (from Muller o calculated), open_pan_evapot = mean_pot_evapot/0.75
	
	number_of_months = n_month
	

	pot_ev_month = mean_pot_evapot
	rain_month = rain
        if(pot_ev_month < rain_month and accTSMD==0. ):
                accTSMD = 0.
        else:
                accTSMD = accTSMD+(rain_month-pot_ev_month)
                accTSMD = np.maximum(maxTSMD, accTSMD)
                accTSMD = np.minimum(accTSMD,0.)

        #Moisture rate modifier (b)
        if(np.abs(accTSMD)<0.444*np.abs(maxTSMD)):
                b=1.
        else:
                b = 0.2+(1.-0.2)*(np.abs(maxTSMD)-np.abs(accTSMD))/(np.abs(maxTSMD)-0.444*np.abs(maxTSMD))	
		#b = np.minimum(1.,b)

	#print 'control moist',b
	#print 'rain',rain 
	#print 'mean_pot_evapot',mean_pot_evapot
	#print 'max',maxTSMD
	#print 'acc',accTSMD
	#print 'denom',(maxTSMD-0.444*maxTSMD)
	
	return b,np.abs(accTSMD)


############################################################
#BIO_HUM partition
############################################################

def BIO_HUM_partition(clay):
	#Partitioning of C between CO2 and BIO+HUM
	X_ratio = 1.67*(1.85+1.6*np.exp(-0.0786*clay*100.))
	BIO_HUM_frac = 1./(X_ratio+1.)
	return BIO_HUM_frac

############################################################
#IOM partition
############################################################

def IOM_partition(SOC):
	#tC/ha
	IOM_frac = 0.049*(SOC**(1.139))
	return IOM_frac

############################################################
#A_matrix
############################################################
def A_matrix(clay):
	
	alpha = BIO_frac_plant*BIO_HUM_partition(clay)
	beta = HUM_frac_plant*BIO_HUM_partition(clay)	

	A_matrix = np.zeros((4,4))
	np.fill_diagonal(A_matrix, -1)	

	A_matrix[2,0]=alpha
	A_matrix[2,1]=alpha
	A_matrix[2,2]=(alpha-1.)
	A_matrix[2,3]=alpha

	A_matrix[3,0]=beta
        A_matrix[3,1]=beta
        A_matrix[3,2]=beta
        A_matrix[3,3]=(beta-1.)

	A_out = A_matrix
	
	return A_out

############################################################
#K_matrix
############################################################

def K_matrix(clay,temp_mean,rain, mean_pot_evapot,c_param,t_param,param_temp_in):
	global moist_rate_b, accTSMD_abs, temp_rate_a
	K_matrix=np.zeros((4,4))
	moist_rate_b, accTSMD_abs = control_moist_func(clay, rain, mean_pot_evapot,c_param)
	temp_rate_a = control_temp_func(temp_mean,param_temp_in)
	K_matrix[0,0]=kDPM*temp_rate_a*moist_rate_b*c_param*t_param
	K_matrix[1,1]=kRPM*temp_rate_a*moist_rate_b*c_param*t_param
	K_matrix[2,2]=kBIO*temp_rate_a*moist_rate_b*c_param*t_param
	K_matrix[3,3]=kHUM*temp_rate_a*moist_rate_b*c_param*t_param
	K_out = K_matrix
	return K_out

############################################################
#Spinup
############################################################
def spinup(SOC_data_init,litterin,clay,temp,rain, pot_evapot,c_param,t_param,param_temp_in):
    #global ABOVE_mean,BELOW_mean, err_above, err_below
    global accTSMD
    global b_rate_array, accTSMD_abs_array
    global a_rate_array, temp_mean_array

    matrix_in_mean=np.append(litterin,[0.,0.])

    #Initialize accTSMD
    accTSMD =0. #check where to place it
    b_rate_array=np.zeros(n_month_ss)
    accTSMD_abs_array=np.zeros(n_month_ss)
    a_rate_array=np.zeros(n_month_ss)
    temp_mean_array=np.zeros(n_month_ss)
    for ts in range(0,n_month_ss):
        temp_mean = temp[ts]
        rain_mean = rain[ts]
        pot_evapot_mean = pot_evapot[ts]
	#Calculate mean K_matrix
        if (ts == 0):
            kk_ma_mean=K_matrix(clay,temp_mean,rain_mean, pot_evapot_mean,c_param,t_param,param_temp_in)
        else:
            kk_ma_mean+=K_matrix(clay,temp_mean,rain_mean, pot_evapot_mean,c_param,t_param,param_temp_in)
	b_rate_array[ts]=moist_rate_b
	accTSMD_abs_array[ts]=accTSMD_abs
	a_rate_array[ts]=temp_rate_a
	temp_mean_array[ts]=temp_mean

    kk_ma_mean=kk_ma_mean/n_month_ss
    a_ma_mean=A_matrix(clay)
    ss_spinup=-np.linalg.solve(np.dot(a_ma_mean,kk_ma_mean),matrix_in_mean)
    IOM_pool = IOM_partition(SOC_data_init)
    ss_spinup_IOM = np.append(ss_spinup,IOM_pool)

    #print 'IN'
    #print matrix_in_mean
    #print 'A'
    #print a_ma_mean
    #print 'K'
    #print kk_ma_mean


    return ss_spinup_IOM


############################################################
#Forward
############################################################
def forward_ACTIVE(SOC_data_init,n_an,init,litterin,clay,temp,water,mean_pot_evapot,c_param,t_param,param_temp_in):
    global temp_mean, rain_mean,accTSMD # needed!!
    global dt,n_month,one_month # not really needed

    matrix_cpools_tmean = np.zeros((n_an-1,5))

    length=n_month
    matrix_cpools_t=np.zeros((length+1,4))
    matrix_in = np.zeros(4)
    matrix_cpools_t[0]=init[0:4]

    for i in range(0,len(litterin)):
        matrix_in[i]=litterin[i]

    a_ma=A_matrix(clay)

    #Initialize accTSMD 
    accTSMD=0.
    #Annual average
    for x in range(0,n_an-1):
        matrix_cpools_ymean = np.zeros(4)
	#Monthly Euler
        for ts in range(x*one_month,(one_month*x)+one_month):
            temp_mean = temp[ts]
            rain_mean = water[ts]
	    evapot_mean = mean_pot_evapot[ts]
            matrix_current= matrix_cpools_t[ts]
            kk_ma = K_matrix(clay,temp_mean,rain_mean,evapot_mean,c_param,t_param,param_temp_in)
            matrix_next = matrix_current + matrix_in + np.dot(a_ma,np.dot(kk_ma,matrix_current))*dt
            matrix_cpools_t[ts+1]=matrix_next
            matrix_cpools_ymean += matrix_next

	#Yearly average C pools		
        matrix_cpools_ymean = matrix_cpools_ymean/one_month

	IOM_pool = IOM_partition(SOC_data_init)
        matrix_cpools_tmean[x] = np.append(matrix_cpools_ymean,IOM_pool)

    return matrix_cpools_tmean

############################################################
# param_opt ######## OBJECTIVE FUNCTION
############################################################
def Jparam_new(param_new_t):
    global ssp, param_predict_c_forward, SOC_data, SOC_var, ss0
    global NEW_ITER, clay, SOC_data_init,litter_inc, temp_in,water_in, PET_in,c_vegetated,t_param, n_an
    global clay,temp_t,water_t,PET_t,ssp_select

    Delta=np.zeros(len(YEARS))
    j=Current_Site_Index

    #Check they change
    #litter income at site (obs tC/ha/month)
    #litter_mean = SITE_litterinc[j]
    #litter_inc in DPM and RPM
    #litter_inc=np.array([litter_mean*frac_array[0],litter_mean*frac_array[1]])

    param_predict_c = spinup(SOC_data_init,litter_inc,clay,temp_in,water_in, PET_in,c_vegetated,t_param,param_new_t)
    ssp = np.sum(param_predict_c)

    param_predict_c_forward = forward_ACTIVE(SOC_data_init,n_an,param_predict_c,litter_inc,clay,temp_t,water_t,PET_t,c_vegetated,t_param,param_new_t)
    ssp_fwd = np.sum(param_predict_c_forward,axis=1)

    ssp_all = np.concatenate((np.array([ssp]),ssp_fwd))

    #selct model values where data avail
    ssp_select=ssp_all[ss0>0]

    #print 'ss0',ss0
    #print 'SOC_data',SOC_data
    #print 'model all', ssp_all
    #print 'model',ssp_select

    Delta[0]=ssp-SOC_data[0]
    Delta[1:]=ssp_select[1:]-SOC_data[1:]
    Delta2w=Delta**2/SOC_var
    Jparam_new = np.sum(Delta2w)

    if ( NEW_ITER % CHI2_PRINT_FREQUENCY == 0 ):
        print "==> NEW_ITER ",NEW_ITER," param_new=",param_new_t," Jparam_new=",Jparam_new
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
	global n_an_4p1000, c_vegetated, DPM_frac_plant, RPM_frac_plant
        global spinup_c_4p1000, predict_c
	global clay,temp_t,water_t,PET_t,c_vegetated,t_param

	litter_new = np.zeros(2)
	litter_new[0] = DPM_frac_plant*in_new
	litter_new[1] = RPM_frac_plant*in_new 
     
	predict_c_pools  = forward_ACTIVE(SOC_data_init,n_an_4p1000,spinup_c_4p1000,litter_new,clay,temp_t,water_t,PET_t,c_vegetated,t_param,optimized_t_param)
	predict_c = np.sum(predict_c_pools,axis=1)

        J_new = abs(target*n_an_4p1000 - np.sum(predict_c[predict_c.shape[0]-1]-predict_c[0]))

        #print 'OBJ FUN', J_new
        #print 'predict_c',predict_c
	#print 'predict_c_pools', predict_c_pools

        NEW_ITER+=1
        if ( NEW_ITER % 100 == 0 ):
                print "NEW_ITER ",NEW_ITER," in_new=",in_new," J_new=",J_new
        param_opt=in_new
        return J_new

#######
#4p1000

################################################
# OPTIMIZATION CONSTRAINTS
#ordine AM,BM,AS,BS
################################################
#def constr1(x):
#        global a_constr,ab_ratio
#        return x[0]+x[2]-a_constr*(x[1]+x[3])*ab_ratio
#
#def constr2(x):
#        global b_constr,ab_ratio
#        return b_constr*(x[1]+x[3])*ab_ratio-x[0]-x[2]


tstart = time.time()

######################################
#Set bounds and constraints for the optimization
#########
bnds=bnds=[(0,10)]

#a_constr=1
#b_constr=1
#ab_ratio_var=(1-a_constr)*100
#con1={'type':'ineq','fun':constr1}
#con2={'type':'ineq','fun':constr2}
#cons=[con1,con2]

############
#Parameters

CHI2_PRINT_FREQUENCY=50
prior_param_temp = 18.27

#Decomposition rate

kDPM=10.
kRPM=0.3
kBIO=0.66
kHUM=0.02


#Since k is vased on yearly decomposition rates
t_param=1./12. #year->month

#Soil cover
c_vegetated = 0.6
c_bare = 1.

#Pools partitioning for plant residues
#For agricultural crops and improved grassland
DPM_frac_plant = 0.59
RPM_frac_plant = 1-DPM_frac_plant

DPM_RPM_frac = np.array([DPM_frac_plant,RPM_frac_plant])

BIO_frac_plant = 0.46
HUM_frac_plant = 1.-BIO_frac_plant

#Initialization
NEW_ITER = 0
one_month = 12
dt=1.

if(optimization_4p1000):
        ##################
        #Import optimized params
        ##################
        n=0 # contatore
        inp=open('temp_param_optim_7.txt','rb')
        while inp.read(1):
                inp.seek(-1,1)
                XX_optimized_param=np.load(inp)
                n+=1

######################################################
#For each site: Set SITE name and experiment duration
#####################################################
ROOTDIR='/Users/ebruni/Desktop/DOTTORATO/'
loc_exp = ROOTDIR+'SIMULAZIONI_ORCHIDEE_SITI/SITES_dataset_Century_new.xlsx'

C_input_exp = pd.read_excel(loc_exp)
site_names_all = C_input_exp['ID.Site'].unique()[2:len(C_input_exp)]
site_names_all = map(str, site_names_all)

N_sites_all=len(site_names_all)

#Control plots
site_T0_array_all=np.array(['CHNO3_Min', 'COL_T0', 'CREC3_Min', 'FEU_T0', 'JEU2_M0', 'LAJA2_Min', 'RHEU1_Min', 'RHEU2_T0','ARAZ_D0_N0', 'ULT_P0_B', 'BROAD_3_Nill', 'FOG_DwN0', 'TREV1_Min','AVRI_T1TR'])


#Stationary solution array for each experiment T0
N_sites=len(site_names_all)

print N_sites

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
SITE_litterinc = np.zeros(N_sites)
SITE_litterinc_err2 = np.zeros(N_sites)
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
SITE_TREATMENTS=[]
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
                SOC_dyn_T=site_T['SOC']*100 #(gC/m2)
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
		temp_in_monthly=DailyToMonthly(temp_in,1)
		SITE_temp_in.append(temp_in_monthly)

        with open(soil_temp) as fileID:
                # C to K
                temp_t = np.array(map(float,fileID))
		temp_t_monthly=DailyToMonthly(temp_t,1)
                SITE_temp_t.append(temp_t_monthly)

	with open(soil_hum_ss) as fileID:
                #rain+snow(mm/day)
                water_in = np.array(map(float,fileID))
		water_in_monthly=DailyToMonthly(water_in,0)
		SITE_water_in.append(water_in_monthly)

        with open(soil_hum) as fileID:
                #rain+snow(mm/day)
		water_t = np.array(map(float,fileID))
		water_t_monthly=DailyToMonthly(water_t,0)
                SITE_water_t.append(water_t_monthly)

        with open(soil_PET_ss) as fileID:
                #potential evapot Pennman-Monteith (mm/s)
                PET_in = np.array(map(float,fileID))
		PET_in_monthly=DailyToMonthly(PET_in,0)*3600.*24.
                SITE_PET_in.append(PET_in_monthly)

        with open(soil_PET) as fileID:
                #potential evapot Pennman-Monteith (mm/s)
                PET_t = np.array(map(float,fileID))
		PET_t_monthly=DailyToMonthly(PET_t,0)*3600.*24.
                SITE_PET_t.append(PET_t_monthly)

        #.....................



        #----------------------------------------------------------------------
        #  determine litter_inc for current site
        #----------------------------------------------------------------------
        # returns ABOVE, BELOW (tC/ha/month) and YEAR of measurements
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

	#Litter calculated as the sum of ABOVE and BELOW
	LITTER_mean = ABOVE_mean + BELOW_mean
	LITTER_var = SITE_ERR2_ABOVE_mean[j]+SITE_ERR2_BELOW_mean[j]

        SITE_litterinc[j]      = LITTER_mean
        SITE_litterinc_err2[j] = LITTER_var


        #mean litter C inputs (tC/ha/month)
        DPM_m = LITTER_mean*DPM_frac_plant
        RPM_m = LITTER_mean*RPM_frac_plant
	
	DPM_var = LITTER_var*DPM_frac_plant*DPM_frac_plant
	RPM_var = LITTER_var*RPM_frac_plant*RPM_frac_plant

        litter_inc = np.array([DPM_m,RPM_m])    # litter C inputs parameters (tC/ha/month)
	Err2_litter_inc = np.array([DPM_var,RPM_var])



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



#MAIN

litterin_sites = np.zeros((N_sites,4))
SOC_out_all = []
optim_param_sites =np.zeros(N_sites)
optim_4p1000_sites = np.zeros((N_sites))
for j in range(N_sites):
        Current_Site_Index=j
        site       = site_names_all[j]
        YEARS      = SOC_clean_year[j]
        SOC_data   = SOC_clean_exp_array[j]
	ss0 = SOC_clean_where[j]
	n_an = np.int(SITE_n_an[j])
	n_an_4p1000 = 30
	n_an_ss = 30
	n_month_ss = 30*12
	if(optimization_param):	
		print 'NUMBER OF YEARS', n_an
		n_month = np.int(n_an*12)
	else:
		print 'NUMBER OF YEARS', n_an_4p1000
		n_month = np.int(n_an_4p1000*12)

	SOC_data_init = SITE_SOC_data_init[j]
	print 'initial SOC',SOC_data_init
        SOC_var    = SOC_clean_exp_variance[j]
        clay       = SITE_clay[j]
        temp_in    = SITE_temp_in[j]
	#print 'TEMP IN', temp_in
        water_in   = SITE_water_in[j]
	#print 'WATER IN', water_in
	PET_in	   = SITE_PET_in[j]
        temp_t     = SITE_temp_t[j]
        water_t    = SITE_water_t[j]
	PET_t     = SITE_PET_t[j]
        err_above = SITE_ERR2_ABOVE_mean[j]
        err_below = SITE_ERR2_BELOW_mean[j]
        ABOVE_mean = SITE_ABOVE_mean[j]
        BELOW_mean = SITE_BELOW_mean[j]
        ab_ratio = ABOVE_mean/BELOW_mean
        frac_array=DPM_RPM_frac
        print 'DPM:RPM fractions',frac_array
        if(optimization_4p1000):
                optimized_t_param = XX_optimized_param[j]

	#LITTER INCOME AT SITE
        #above-below array to calculate uncertainties
        AB_BE_array=np.array([ABOVE_mean,BELOW_mean])

        #litter income at site (obs tC/ha/month)
        litter_mean = SITE_litterinc[j]
        Err2_litter_mean = SITE_litterinc_err2[j]

	#litter_inc in DPM and RPM
	litter_inc=np.array([litter_mean*frac_array[0],litter_mean*frac_array[1]])

        #to be saved
        #litter_inc_save=np.append(litter_inc,LITTER_mean)
        #litter_inc_err_save=np.append(np.sqrt(Err2_litter_inc),np.sqrt(tot_err2_litter_inc))

        #litter income prior (1D) ->splitted in J_new
        in_opt = litter_mean*(1+0.004)

	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        print '>>>>>>>>>>>>>> Analysis for SITE ',site
        print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        #--------------------------------------------------------------#
	if(optimization_param):
        	spinup_c=spinup(SOC_data_init,litter_inc,clay,temp_in,water_in, PET_in,c_vegetated,t_param,prior_param_temp)
		#Standard forward (non optimzed)
        	fwd=forward_ACTIVE(SOC_data_init,n_an,spinup_c,litter_inc,clay,temp_t,water_t,PET_t,c_vegetated,t_param,prior_param_temp)
		SITE_SOC_model=np.sum(spinup_c) #tC/ha
		print 'Stationary solution before optimiz: '
		print SITE_SOC_model
		#print 'Site SOC dynamics before optimiz: '
		#print fwd

		Dynamics_all=np.concatenate((np.array([spinup_c]),fwd))
		print '  DPM         RPM         BIO         HUM         IOM'
		print Dynamics_all
		print "--"
		Dynamics_select_exper = Dynamics_all[ss0>0]
		Dynamics_select_exper_sum = np.sum(Dynamics_select_exper,axis=1)
		print Dynamics_select_exper_sum

	if(optimization_4p1000):
                spinup_c_4p1000=spinup(SOC_data_init,litter_inc,clay,temp_in,water_in, PET_in,c_vegetated,t_param,optimized_t_param)
		SITE_SOC_model=np.sum(spinup_c_4p1000) #tC/ha
                print 'Stationary solution before 4p1000: '
                print SITE_SOC_model
                #Standard forward (non optimzed)
                #fwd_4p1000=forward_ACTIVE(SOC_data_init,n_an_4p1000,spinup_c,litter_inc,clay,temp_t,water_t,PET_t,c_vegetated,t_param,optimized_t_param)
        #--------------------------------------------------------------#
	
	if(plot_CM):
	#Plot control moisture
		plt.figure(figsize=(15,5))
		plt.scatter(accTSMD_abs_array,b_rate_array)
		plt.vlines(np.abs(maxTSMD),0.,1.,colors='k',linestyles='dashed')
		plt.vlines(0.444*np.abs(maxTSMD),0.,1.,colors='r',linestyles='dashed')
		plt.xlabel('accumulated absolute topsoil moisture deficit (mm/month)')
		plt.ylabel('moisture rate modifying factor (b)')
		plt.title(site)
		plt.xlim((0,60))
		plt.savefig(pp0,format='pdf')
   		#plt.show()

	if(plot_CT):
	        #Plot control moisture
                plt.figure(figsize=(15,5))
                plt.scatter(temp_mean_array,a_rate_array)
                plt.xlabel('Monthly mean temperature (C)')
                plt.ylabel('temperature rate modifying factor (a)')
                plt.title(site)
                plt.xlim((-10,30))
                plt.savefig(pp1,format='pdf')
                plt.show()
	

	print 'LITTER:',litter_inc
	print 'CLAY:',clay
	print 'TEMP:',np.mean(temp_t)
	print 'WATER:',np.mean(water_t)

	#Years, Standard, Measured SOC
        #SOCXX = np.stack((YEARS,Dynamics_select_exper_sum,SOC_data))
        #np.save(out,SOCXX)

	#OPTIMIZATION PARAMETERS
	if(optimization_param):
		print 'OPTIMIZATION PARAM'
		bnds_param=bnds=[(10.,40.)]
		opt_param_mean=minimize(Jparam_new, prior_param_temp, method='SLSQP', bounds=bnds_param, options={'disp':True})
		print "SLSQP: Optimum solution param_temp:",opt_param_mean.x
		param_min=opt_param_mean.x
		CHI2=Jparam_new(opt_param_mean.x)
		print 'SLSQP: CHI2 ',CHI2
		ssp_o = np.array([ssp])
		ssp_o_fwd = np.sum(param_predict_c_forward,axis=1)

		#SOC_model = np.concatenate((ssp_o,ssp_o_fwd))
		SOC_model_optim = ssp_select
		SOC_error=SOC_var
		print 'SOC_error',SOC_error
		#YEARS, SOC_model_standard, SOC_model_optim, SOC_data
		XX=np.stack((YEARS,Dynamics_select_exper_sum,SOC_model_optim,SOC_data,SOC_error))
		np.save(out1,XX)

		print "SOC model optim", SOC_model_optim
		print "len SOC model optim:",len(SOC_model_optim)
		print "len SOC data",len(SOC_data)
		print "SOC EXPER", SOC_data
		print "Optimized SOC"
		print " YEAR      SOCmodelOptim     	  SOCdata"
		for k in range(len(SOC_model_optim)-1):
			print '{0:6d} {1:7.1f} {2:7.1f}'.format(int(YEARS[k]), SOC_model_optim[k],SOC_data[k])
		optim_param_sites[j]=param_min

	#-------------------------
	#-------------------------
	#4p1000
	if(optimization_4p1000):
		target = SITE_SOC_model*0.004
		
		#SITE_SOC_dyn = np.concatenate((np.array([SITE_SOC_model]),fwd))
		#print 'SOC dynamics before opt', SITE_SOC_dyn

		print ' '
		opt_mean=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,options={'disp':True})
		litter_opt = opt_mean.x
		print "SLSQP: Optimum solution:", litter_opt

		#total_opt_in=np.sum(opt_mean.x)

		#optimized litter pools and total litter (save)
		#in_opt_save = np.append(litter_opt,total_opt_in)

                print 'Initial litter (tC/ha/month):'
                print litter_mean
                print '4p1000 litter (tC/ha/month):'
                print litter_opt

		#calculate percentage increase/decrease of inputs
		input_change=(litter_opt-litter_mean)/litter_mean
		print "% change of litter inputs:",input_change


		#Check 4p1000 target is reached
		END=predict_c.shape[0]-1
		print 'END',END
		C_fin=np.sum(predict_c[END])
		#C_init=np.sum(predict_c[0])
		C_init=SITE_SOC_model
		SUMO=C_fin-C_init
		Target_reached=(C_fin-C_init)/(C_init*n_an_4p1000)
		if (Target_reached < 0.005 and Target_reached > 0.003):
			print "Target reached successfully"
			print "Target reached :", Target_reached
		else:
			print "Target not reached"
			print 'C init', C_init
			#print predict_c[0]
			print SITE_SOC_model
			print 'C_fin', C_fin

		optim_4p1000_sites[j]=litter_opt

if(optimization_param):
        for i in range(N_sites):
                print site_names_all[i],' optimized param: ', optim_param_sites[i]
	out3=open("temp_param_optim_7.txt","wb")
        np.save(out3,optim_param_sites)
        out3.close()
        out1.close()

if(optimization_4p1000):
        for i in range(N_sites):
                print site_names_all[i],' C inputs to 4p1000: ', optim_4p1000_sites[i]
	out4=open("Cinput_4p1000.txt","wb")
        np.save(out4,optim_4p1000_sites)
        out4.close()
        out2.close()

if(plot_CM):
	pp0.close()
if(plot_CT):
        pp1.close()

