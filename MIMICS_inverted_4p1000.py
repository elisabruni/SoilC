#MIMICS inverted 4p1000
#Wieder et al. 2015
#see https://github.com/wwieder/MIMICS

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
loc_out = ROOTDIR+"SCRIPT_MODELLI/MULTIMODEL/OUTPUTS_4p1000_v5/MIMICS/"
loc_exp = ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/SITES_dataset_multimodel.xlsx'
loc_optim_param = ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM5/'
loc_spinup = ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM5/'


#OUTPUT FILES
out_mo_pools=open(loc_out+"SOC_MIMICS_pools.txt","wb")
out_mo=open(loc_out+"SOC_MIMICS.txt","wb")
out_lit=open(loc_out+"Litter_income_MIMICS.txt","wb")
out_lit_tot=open(loc_out+"Litter_tot_MIMICS.txt","wb")
out_priors = open(loc_out+"priors_and_opt_in_MIMICS.txt","wb")

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
        aa=np.asarray(site_T0['ABOVE']).astype(np.float64)*100./365.  # array of C_above (gC/m2/day)
        bb=np.asarray(site_T0['BELOW']).astype(np.float64)*100./365.  # array of C_below (gC/m2/day)
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
        global predict_c
        global fCLAY,temp_t
        global state, t_fwd


        predict_c = odeint(XEQ_fw,state,t_fwd,args=(in_new,))
	#print 'PREDICT 0'
	#print predict_c[0]
	#print 'PREDICT END'
	#print predict_c[predict_c.shape[0]-1]
	#print 'somma della differenza END - 0'
	#print np.sum(predict_c[predict_c.shape[0]-1]-predict_c[0])

	#print 'TARGET',target

        J_new = abs(target*n_an_4p1000 - np.sum(predict_c[predict_c.shape[0]-1]-predict_c[0]))

        print 'OBJ FUN', J_new
        #print 'predict_c',predict_c

        NEW_ITER+=1
        if ( NEW_ITER % 100 == 0 ):
                print "NEW_ITER ",NEW_ITER," in_new=",in_new," J_new=",J_new
        param_opt=in_new
        return J_new


#----------------forward-----------------------

def XEQ_fw(y, t, forc_npp): #y=st_state , t=time_step

        global depth,Vslope,Vint,Kslope,Kint
        global CUE,k,a,cMAX,cMIN,cSLOPE
        global KO,FI
        global aV,aK,fMET
        global LITmin,MICtrn,SOMmin,DEsorb,OXIDAT
        global fCLAY, pH
	global temp_t
	global Tao_MOD1
	#print t
        t=np.int(t)

        LIT_1=y[0]
        LIT_2=y[1]
        MIC_1=y[2]
        MIC_2=y[3]
        SOM_1=y[4]
        SOM_2=y[5]
        SOM_3=y[6]
 
	#Forward
    	TSOI=temp_t[t]
	#print 'TEMPERATURE:',TSOI, ' t:',t
   
	###################
        #EST_LIT_in  = forc_npp                  #gC/m2/day (Knapp et al. Science 2001)
        #EST_LIT     = EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/day
        #I        = np.zeros(2)              #Litter inputs to MET/STR
        #I[0]     = (EST_LIT / depth) * fMET      #partitioned to layers (mgC/cm3/day)
        #I[1]     = (EST_LIT / depth) * (1-fMET )
        #Tao_MOD1 = np.sqrt(forc_npp*365./100.)  #basicaily standardize against NWT (forc_npp*365=ANPP (gC/m2/y))
        tao      = np.array([5.2e-4*np.exp(0.3*fMET ), 2.4e-4*np.exp(0.1*fMET )])
        tao      = tao * Tao_MOD1 *24. #daily
        fPHYS    = np.array([0.3 * np.exp(1.3*fCLAY), 0.2 * np.exp(0.8*fCLAY)])  #fraction to SOMp
        fCHEM    = np.array([0.1 * np.exp(-3.*fMET )  , 0.3 * np.exp(-3.*fMET )])   #fraction to SOMc
        fAVAI    = 1.- (fPHYS + fCHEM)
        #desorb   = 9e-4 * np.exp(-3.*(np.sqrt(fCLAY))) #if modified by MIC!
        #desorb   = 3e-4 * np.exp(-4.*(np.sqrt(fCLAY))) #if stand alone rate constant
        desorb   = 1.5e-5 * np.exp(-1.5*(fCLAY))*24. #daily      #CHANGED FOR GLOBAL RUN!!!
        pSCALAR  = a * np.exp(-k*(np.sqrt(fCLAY)))  #Scalar for texture effects on SOMp
        MOD1     = np.array([10, 2, 10, 3, 3, 2])
        MOD2     = np.array([8, 2 ,4 * pSCALAR, 2, 4, 6 * pSCALAR])

        Vmax     = np.exp(TSOI * Vslope + Vint) * aV *24.
        Km       = np.exp(Kslope * TSOI + Kint) * aK
        VMAX     = Vmax * MOD1
        KM       = Km / MOD2
        #????????????????????

        #Flows to and from MIC_1
        #EQA1
        LITmin[0] = MIC_1 * VMAX[0] * LIT_1 / (KM[0] + LIT_1)   #MIC_1 decomp of MET lit
        #EQA2
        LITmin[1] = MIC_1 * VMAX[1] * LIT_2 / (KM[1] + LIT_2)   #MIC_1 decomp of STRUC lit
        #EQA4
        MICtrn[0] = MIC_1 * tao[0]  * fPHYS[0]                  #MIC_1 turnover to PHYSICAL SOM
        MICtrn[1] = MIC_1 * tao[0]  * fCHEM[0]                  #MIC_1 turnover to CHEMICAL SOM
        MICtrn[2] = MIC_1 * tao[0]  * fAVAI[0]                  #MIC_1 turnover to AVAILABLE SOM
        #EQA3
        SOMmin[0] = MIC_1 * VMAX[2] * SOM_3 / (KM[2] + SOM_3)   #decomp of SOMa by MIC_1

        #Flows to and from MIC_2
        #EQA5
        LITmin[2] = MIC_2 * VMAX[3] * LIT_1 / (KM[3] + LIT_1)   #decomp of MET litter
        #EQA6
        LITmin[3] = MIC_2 * VMAX[4] * LIT_2 / (KM[4] + LIT_2)   #decomp of SRUCTURAL litter

        #EQA8
        MICtrn[3] = MIC_2 * tao[1]  * fPHYS[1]                  #MIC_2 turnover to PHYSICAL  SOM
        MICtrn[4] = MIC_2 * tao[1]  * fCHEM[1]                  #MIC_2 turnover to CHEMICAL  SOM
        MICtrn[5] = MIC_2 * tao[1]  * fAVAI[1]                  #MIC_2 turnover to AVAILABLE SOM

        #EQA7
        SOMmin[1] = MIC_2 * VMAX[5] * SOM_3 / (KM[5] + SOM_3)   #decomp of SOMa by MIC_2	

 
        #EQA9
        DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)          #desorbtion of PHYS to AVAIL (function of fCLAY)

        #EQA10
        OXIDAT    = ((MIC_2 * VMAX[4] * SOM_2 / (KO[1]*KM[4] + SOM_2)) +
                           (MIC_1 * VMAX[1] * SOM_2 / (KO[0]*KM[1] + SOM_2)))  #oxidation of C to A

        #can make fluxes from CHEM a function of microbial biomass size?

        dLIT_1 = forc_npp[0]*(1-FI[0]) - LITmin[0] - LITmin[2]
        dMIC_1 = CUE[0]*(LITmin[0]+ SOMmin[0]) + CUE[1]*(LITmin[1]) - np.sum(MICtrn[0:3])
        dSOM_1 = forc_npp[0]*FI[0] + MICtrn[0] + MICtrn[3]- DEsorb

        dLIT_2 = forc_npp[1] * (1.-FI[1]) - LITmin[1] - LITmin[3]
        dMIC_2 = CUE[2]*(LITmin[2]+ SOMmin[1]) + CUE[3]*(LITmin[3]) - np.sum(MICtrn[3:6])
        dSOM_2 = forc_npp[1]*FI[1] + MICtrn[1] + MICtrn[4] - OXIDAT

        dSOM_3  = MICtrn[2] + MICtrn[5] + DEsorb + OXIDAT - SOMmin[0] - SOMmin[1]

        diff_eq = [dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3]
    
    	fluxes=[LITmin,KM, LIT_1,LIT_2,MICtrn,MIC_2,tao,fAVAI,fCHEM,fPHYS,VMAX,SOMmin,SOM_3,DEsorb,OXIDAT]
	for flux in fluxes:
		if(np.any(np.isnan(flux))):
			print "Il y a des NA"
			exit()

	return diff_eq		


############################################################
#Initialization
tstart = time.time()
NEW_ITER = 0
CHI2_PRINT_FREQUENCY=2

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

SITES_spinup=pd.read_csv(loc_spinup+'list_SOC_spinup_opti_MIMICS.csv')
SITES_spinup.drop(SITES_spinup.columns[[0,1]], axis = 1, inplace = True)

#Import optimized parameters
my_cols = ['param','value']
SITES_optim_param=pd.read_csv(loc_optim_param+'optim_param_MIMICS.csv', names=my_cols, engine="python")

SITE_aV=SITES_optim_param['value'][(SITES_optim_param['param']=='1')] #aV
SITE_aK=SITES_optim_param['value'][(SITES_optim_param['param']=='2')] #aK
SITE_fMET=SITES_optim_param['value'][(SITES_optim_param['param']=='3')] #fMET

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
#SITE_temp_t=[]
SITE_TSOI_fw=[]
SITE_ABOVE_mean=np.zeros(N_sites)
SITE_BELOW_mean=np.zeros(N_sites)
SITE_ERR2_ABOVE_mean=np.zeros(N_sites)
SITE_ERR2_BELOW_mean=np.zeros(N_sites)
SITE_cov_mean= []
SITE_TREATMENTS=[]
#SITE_BD=np.zeros(N_sites)
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

        clay = np.mean(site_T0['Clay'])
        #bd=np.mean(site_T0['Bulk density'])*1000. #kgsoil/m3
        ph=np.mean(site_T0['pH'])

        #SITE_BD[j]=bd
        SITE_pH[j]=ph
        SITE_clay[j]=clay

        SITE_date_init[j]=date_init
        SITE_date_end[j]=date_end

        soil_temp    = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"temp_"+site+"_"+date_init+"_"+date_end+".txt"

        with open(soil_temp) as fileID:
                # C
                temp_t = np.array(map(float,fileID))
                #SITE_temp_t.append(temp_t)

	TSOI_mean_fw = DailyToYearly(temp_t,1)

	#Repeate mean annual temperature 365 times each year
	append_temp = []
	for y_t in TSOI_mean_fw:
		rep = [y_t]*365
		append_temp.append(rep)
	TSOI_fw = np.concatenate(append_temp)
	SITE_TSOI_fw.append(TSOI_fw)

        #.....................
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
        Err2_litter_inc_AB = np.array([SITE_ERR2_ABOVE_mean,SITE_ERR2_BELOW_mean])

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
bnds=[(0,10),(0,10)]

################################################
#Add constraint to Met:struc ratio (not to be lower or higher than 50% old ratio
def constr1(x):
        global ms_constr, ms_ratio
        return x[0]/x[1]-ms_constr*ms_ratio
def constr2(x):
        global ms_constr, ms_ratio
        return -(x[0]/x[1])+(1+ms_constr)*ms_ratio

ms_constr=0.5
con1={'type':'ineq','fun':constr1}
con2={'type':'ineq','fun':constr2}
cons=[con1,con2]

#WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
#Parameters
#WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

depth=30. #cm
#-----------------caclulate parameters---------------------------
#Calculate Vmax & (using parameters from German 2012, as in Wieder et al. 2013 Nature Climate Change)
Vslope   = np.array([0.063]*6) #daily (6 Dimensional array)
Vint    = 5.47 #daily
#aV       = 8e-6

Kslope   = np.zeros(6)
Kslope[0]= 0.017 #META LIT to MIC_1
Kslope[1]= 0.027 #STRU LIT to MIC_1 
Kslope[2]= 0.017 #AVAI SOM to MIC_1 
Kslope[3]= 0.017 #META LIT to MIC_2
Kslope[4]= 0.027 #STRU LIT to MIC_2
Kslope[5]= 0.017 #AVAI SOM to MIC_2
Kint     = 3.19
#aK       = 10

CUE        =np.array([0.55, 0.25, 0.75, 0.35])  #for LITm and LITs entering MICr and MICK, respectively

#------NEW Parameters--------------
k        = 2.0    #2.0			#REDUCED FROM 3 TO 1, REDUCES TEXTURE EFFECTS ON SOMa decay
a        = 2.0    #2.2			#increased from 4.0 to 4.5

cMAX     = 1.4                    #ORIG 1.4 Maximum CHEM SOM scalar w/   0% Clay 
cMIN     = 1.2                    #ORIG 1.4 Minimum CHEM SOM scalar w/ 100% Clay 
cSLOPE   = cMIN - cMAX            #Slope of linear function of cSCALAR for CHEM SOM  

#------------!!MODIFIERS AS IN MIMICS2_b!!---------------

KO       = np.array([4,4])      #scalar modifies Km of Oxidat	
FI       = np.array([0.05, 0.05])

LITmin  = np.zeros(4)
MICtrn  = np.zeros(6)
SOMmin  = np.zeros(2)
DEsorb  = np.zeros(1)
OXIDAT  = np.zeros(1)

#------------#------------#------------

litterin_sites = np.zeros((N_sites,4))
SOC_out_all = []
optim_4p1000_sites = np.zeros((N_sites,2))
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
        fCLAY       = SITE_clay[j]

        pH      = SITE_pH[j]
        temp_t     = SITE_TSOI_fw[j]

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


        #Spinup
        LIT_1_ss=SITES_spinup.loc[j]['LIT_1']
        LIT_2_ss=SITES_spinup.loc[j]['LIT_2']
        MIC_1_ss=SITES_spinup.loc[j]['MIC_1']
        MIC_2_ss=SITES_spinup.loc[j]['MIC_2']
        SOM_1_ss=SITES_spinup.loc[j]['SOM_1']
	SOM_2_ss=SITES_spinup.loc[j]['SOM_2']
	SOM_3_ss=SITES_spinup.loc[j]['SOM_3']

        state=np.array([LIT_1_ss,LIT_2_ss,MIC_1_ss,MIC_2_ss,SOM_1_ss,SOM_2_ss,SOM_3_ss])

        #Optimized parameters
        aK = SITE_aK.iloc[j]
        aV = SITE_aV.iloc[j]
        fMET = SITE_fMET.iloc[j]

	#Transform litter inputs in mgC/cm3/day and divide in met:str
        EST_LIT_in  = litter_mean                  #gC/m2/day (Knapp et al. Science 2001)
        EST_LIT     = EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/day
        I        = np.zeros(2)              #Litter inputs to MET/STR
        I[0]     = (EST_LIT / depth) * fMET      #partitioned to layers (mgC/cm3/day)
        I[1]     = (EST_LIT / depth) * (1-fMET )

	ms_ratio = I[0]/I[1] #metabolic:structural ratio

	Tao_MOD1 = np.sqrt(litter_mean*365./100.)  #basicaily standardize against NWT (forc_npp*365=ANPP (gC/m2/y))


        #Error of the metabolic and struc fractions
        Err2_litter_inc_M = Err2_litter_mean*fMET
        Err2_litter_inc_S = Err2_litter_mean*(1-fMET)
        Err2_litter_inc = np.array([Err2_litter_inc_M,Err2_litter_inc_S])

        #to be saved
        litter_inc_save=np.append(I*365.*depth/10.,litter_mean*365./100.) #tC/ha/yr
        litter_inc_err_save=np.append(np.sqrt(Err2_litter_inc),np.sqrt(Err2_litter_mean))

        #litter income prior (2D)
        in_opt = I*(1+0.004)

        #4p1000 analysis
        print '>>>>>>>>>>>>>>>>'
        print '4p1000 ',site
        target = np.sum(state)*0.004

        t_fwd=np.arange(0,n_days_4p1000-365) #Togli un anno perche dopo aggiungi spinup
        print ' '


        #opt_mean=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,constraints=cons,options={'disp':True})
	opt_mean=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,options={'disp':True})
        litter_opt = opt_mean.x
        print "SLSQP: Optimum solution:", litter_opt
        total_opt_in=np.sum(opt_mean.x)*depth*1e4/1e3 #gC/m2/day

        print 'Initial litter (gC/m2/day):'
        print litter_mean
        print '4p1000 litter (gC/m2/day):'
        print total_opt_in


        #optimized litter pools and total litter (save)
        litter_opt_save = np.append(litter_opt*365.*depth/10.,total_opt_in*365./100.) #tC/ha/year
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

        predict_c_standard_pools_day = odeint(XEQ_fw,state,t_fwd,args=(I,))
        predict_c_standard_pools=DailyToYearly(predict_c_standard_pools_day,1)
        predict_c_standard = np.sum(predict_c_standard_pools,axis=1)

        SOC_model_standard_pools = np.concatenate(([state],predict_c_standard_pools))
        SITE_SOC_model=np.sum(state)
        SOC_model_standard = np.concatenate((np.array([SITE_SOC_model]),predict_c_standard))

        #4p1000
        predict_c_opt_pools_day=odeint(XEQ_fw,state,t_fwd,args=(litter_opt,))
        predict_c_opt_pools=DailyToYearly(predict_c_opt_pools_day,1)
        predict_c_opt = np.sum(predict_c_opt_pools,axis=1)

        SOC_model_opt_pools =np.concatenate(([state],predict_c_opt_pools))
        SOC_model_opt = np.concatenate((np.array([SITE_SOC_model]),predict_c_opt))


        year_out = np.arange(1,31)

        SOC_pools_out=np.stack((SOC_model_standard_pools*depth/10.,SOC_model_opt_pools*depth/10.))#tC/ha
        #print len(year_out)
        #print len(SOC_model_standard/100.)
        #print len(SOC_model_opt/100.)
        SOC_out=np.stack((year_out,SOC_model_standard*depth/10.,SOC_model_opt*depth/10.)) #tC/ha
        np.save(out_mo_pools,SOC_pools_out)
        np.save(out_mo,SOC_out)


        ############################
        #UNCERTAINTIES
        ###########################
        Uncert_Q = True
        if(Uncert_Q):

                MC_length=2 #set the number of Monte Carlo simulations

                #optimize for n variations of in_opt generatated randomly around the above/below covariance
                opt_parameters_MC=np.zeros((MC_length,len(I)))

                #prior above_below
                AB_BE_array=np.array([ABOVE_mean,BELOW_mean])
                ab_be_init_est=AB_BE_array*(1+0.004)

                cov_AB_mean=cov_mean #covariance between ABOVE mean and BELOW mean if no Nans
                print 'cov',cov_AB_mean
                if (np.all(cov_AB_mean)==0): #if covariance is 0, take the mean amongst all sites' covariances to generate cov_AB_mean
                        cov_AB_mean=np.mean(SITE_cov_mean)
                #print cov_AB_mean

                in_rand_param_MC = np.zeros((MC_length,len(I)))
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
				in_rand_sum=np.sum(in_rand)	

        			#Transform litter inputs in mgC/cm3/day and divide in met:str
       				EST_LIT_in  = in_rand_sum                  #gC/m2/day (Knapp et al. Science 2001)
        			EST_LIT     = EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/day
        			I_rand        = np.zeros(2)              #Litter inputs to MET/STR
        			I_rand[0]     = (EST_LIT / depth) * fMET      #partitioned to layers (mgC/cm3/day)
        			I_rand[1]     = (EST_LIT / depth) * (1-fMET )

				in_rand_param=I_rand #array to optimize

                                #Save all priors in an array
				print '***********'
				print 'in_rand_param'
				print in_rand_param
				print 'in_rand_param_MC[sample_shape]'
				print in_rand_param_MC[sample_shape]
				print '***********'
                                in_rand_param_MC[sample_shape]=in_rand_param #add new generated sample to array on rand_in samples
                                print '****************'
                                print "LITTER IN generated randomly from ab_be_init_est:", in_rand
                                print 'total litter in (gC/m2/day):',in_rand_sum
                                print '****************'
                                #Minimize J_new for the generated array
                                #opt=minimize(J_new, in_rand_param, method='SLSQP',bounds=bnds, constraints=cons,options={'disp':True})
				opt=minimize(J_new, in_rand_param, method='SLSQP',bounds=bnds,options={'disp':True})
                                opt_param=opt.x
                                opt_parameters_MC[sample_shape]=opt_param

                                print "OPTIMUM LITTER IN (gC/m2/day)"
                                print opt_param*depth*1e4/1e3
                                print '****************'
                                print "Litter in increase (%)"
                                print (opt_param-in_rand_param)*100./in_rand_param
                                print '****************'
                                sample_shape+=1
                out_priors_and_opt = np.stack((in_rand_param_MC*365.*depth/10.,opt_parameters_MC*365.*depth/10.)) #save as tC/ha/yr
                np.save(out_priors,out_priors_and_opt)

                #STANDARD ERROR CALCULATION for the optimized litter inputs to reach 4x1000
                #litter_opt_err=np.std(opt_parameters_MC[:])/np.sqrt(MC_length)

                litter_opt_M=np.std(opt_parameters_MC[:,0])/np.sqrt(MC_length)
                litter_opt_S=np.std(opt_parameters_MC[:,1])/np.sqrt(MC_length)
                litter_opt_err = np.array([litter_opt_M,litter_opt_S])
                litter_opt_err_sum = np.sum(litter_opt_err)

		litter_opt_err_save = np.append(litter_opt_err*365.*depth/10.,litter_opt_err_sum*365.*depth/10.) #tC/ha/yr

                #Error litter = SE per litter in e in opt
		save_lit = np.stack((litter_inc_save,litter_inc_err_save,litter_opt_save,litter_opt_err_save))#save as tC/ha/yr
                np.save(out_lit,save_lit)

                litter_tot_save =  np.stack((litter_mean*365./100., np.sqrt(Err2_litter_mean)*365./100.,total_opt_in*365./100.,litter_opt_err_sum*365.*depth/10.)) #tC/ha/yr
                np.save(out_lit_tot,litter_tot_save)


#Print optimized parameters
for i in range(N_sites):
        print site_names_all[i],' C inputs to 4p1000 (tC/ha/yr): ', optim_4p1000_sites[i]*365.*depth/10.                           

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

