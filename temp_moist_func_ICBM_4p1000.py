import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

########
#Param
########
#--------
#Moisture
#--------
gamma=0.85 #From graph in Karlsson et al 2011
rmin=0.55 #From graph in Karlsson et al 2011
alpha=0.5 #From Fortin et al. 2011
#----------
#Temperature
#----------
Tmin=-3.8 #From Karlsson et al 2011
#----------

###################
#Moisture function
###################
def rwat2(clay, water_in_m3):

        # mcw = wilting point
        # mcfc = field capacity

	#From katterer et al 2006
	mcfc=0.2699+0.3247*clay
        mcw=0.0284+0.3790*clay

	#water_in_m3 = np.maximum(water_in_m3,alpha*mcw)

        if(water_in_m3<alpha*mcw):
                re_wat=0.18*(water_in_m3/(alpha*mcw))

        elif((water_in_m3>alpha*mcw or water_in_m3==alpha*mcw) and (water_in_m3<gamma*mcfc or water_in_m3==gamma*mcfc)):
                re_wat=0.18+(1.-0.18)*((water_in_m3-alpha*mcw)/(gamma*mcfc-alpha*mcw))

        elif(water_in_m3>gamma*mcfc):
                re_wat=1+(1-rmin)*((water_in_m3-gamma*mcfc)/(gamma*mcfc-mcfc))

        return re_wat


def func_temp(Tsoil,Tmin):
        temp_func=np.maximum(0.,((Tsoil-Tmin)**2)/((30-Tmin)**2))
        return temp_func


######################################################
#For each site: Set SITE name and experiment duration
#####################################################
ROOTDIR='/Users/ebruni/Desktop/DOTTORATO/'
loc_exp = ROOTDIR+'SCRIPT_MODELLI/MULTIMODEL/SITES_dataset_multimodel.xlsx'

OUTPUT_files=ROOTDIR+'SCRIPT_MODELLI/IBCM/TEMP_MOIST_RESP_FUNC'

C_input_exp = pd.read_excel(loc_exp)
site_names = C_input_exp['ID.Site'].unique()[2:len(C_input_exp)]
site_names = map(str, site_names)

N_sites=len(site_names)

site_T0_array=np.array(['CHNO3_Min', 'COL_T0', 'CREC3_Min', 'FEU_T0', 'LAJA2_Min', 'RHEU1_Min', 'RHEU2_T0','ARAZ_D0_N0', 'ULT_P0_B', 'BROAD_3_Nill', 'TREV1_Min','AVRI_T1TR','BOLO_T0', 'GRAB_CP','MUNCHE_CP','RITZ_CP'])


SITE_av_re_day_fw=np.zeros(N_sites)
SITE_av_re_day_ss=np.zeros(N_sites)

j=0
for site in site_names:

        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        print "READING DATA OF SITE: ",site
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        site_df = C_input_exp[(C_input_exp['ID.Site'].values == [site])]
        year_0 = np.min(site_df['Year'])
        year_end = np.max(site_df['Year'])
        year_30=year_0+30

        site_T0_name = site_T0_array[j]

        #GET initial years for ss and forward
        site_T0= site_df[(site_df['ID.Treatment'].values == [site_T0_name])]
        date_init = np.str(1980)
        date_end = np.str(2010)
        date_init_ss = np.str(np.int(year_0 - 30))
        date_end_ss = np.str(np.int(year_0 - 1))
        exper_len = np.int(date_end) - np.int(date_init) + 1

        soil_temp_ss = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"temp_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
        soil_hum_ss  = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"hum_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
        soil_temp    = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"temp_"+site+"_"+date_init+"_"+date_end+".txt"
        soil_hum     = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"hum_"+site+"_"+date_init+"_"+date_end+".txt"

        with open(soil_temp_ss) as fileID:
                # C to K
                temp_in = np.array(map(float,fileID))

        with open(soil_temp) as fileID:
                # C to K
                temp_t = np.array(map(float,fileID))

        with open(soil_hum_ss) as fileID:
                #  conversion kg(H2O)/m2(soil) to m3(H2O)/m3(soil)
                water_in_m3 = np.array(map(float,fileID))/100

        with open(soil_hum) as fileID:
                #  conversion kg(H2O)/m2(soil) to m3(H2O)/m3(soil)
                water_t_m3 = np.array(map(float,fileID))/100

        clay = np.mean(site_T0['Clay'])


	#Moisture response funciton parameter
	#cum_re_wat_ss=0.
	#for moist_day_ss in water_in_m3:	#steady state
	#	re_wat_ss = rwat2(clay, moist_day_ss)
	#	cum_re_wat_ss+=re_wat_ss

	#Temperaure response funciton parameter
	#cum_re_temp_ss=0.
	#for temp_day_ss in temp_in:		#steady state
	#	re_temp_ss = func_temp(temp_day_ss,Tmin)
	#	cum_re_temp_ss+=re_temp_ss
		
	re_day_ss = np.zeros(len(water_in_m3))
        for i in range(0,len(water_in_m3)):
                #Steady state
                moist_day_ss=water_in_m3[i]
                temp_day_ss=temp_in[i]

                re_wat_ss = rwat2(clay, moist_day_ss)
                re_temp_ss = func_temp(temp_day_ss,Tmin)

                re_dayS = re_wat_ss*re_temp_ss
                re_day_ss[i]=re_dayS


	re_day_fw = np.zeros(len(water_t_m3))
	for i in range(0,len(water_t_m3)):
		#Forward
		moist_day_fw=water_t_m3[i]
		temp_day_fw=temp_t[i]

		re_wat_fw = rwat2(clay, moist_day_fw)
		re_temp_fw = func_temp(temp_day_fw,Tmin)

		re_day = re_wat_fw*re_temp_fw
		re_day_fw[i]=re_day

	#Average re_day
        av_re_day_ss = np.mean(re_day_ss)
	av_re_day_fw = np.mean(re_day_fw)

	print 'Daily average re_day over entire experim'
	print av_re_day_fw	
	#cum_re_day_fw = np.sum(re_day_fw)
	#print 'Cumulative re_day over entire experim'
	#print cum_re_day_fw
	#print 'Annual cumulative average (sum of re_day over experiment divided by number of years)'
	#ann_cum_av_re_day_fw = cum_re_day_fw/exper_len 
	#print ann_cum_av_re_day_fw

	SITE_av_re_day_fw[j]=av_re_day_fw
        SITE_av_re_day_ss[j]=av_re_day_ss
	j+=1


df = pd.DataFrame({'site':site_names,'av_re_day_fw':SITE_av_re_day_fw,'av_re_day_ss':SITE_av_re_day_ss})

ULTU_re_day_ss = df.at[9,'av_re_day_ss']
ULTU_re_day_fw = df.at[9,'av_re_day_fw']
print ' '

df['norm_av_re_day_fw'] = df['av_re_day_fw']/ULTU_re_day_fw
df['norm_av_re_day_ss'] = df['av_re_day_ss']/ULTU_re_day_ss
df['norm_av_re_day_ss_fw'] = df['av_re_day_ss']/ULTU_re_day_fw #SS normalized against fw ULTUNA

print 'Average re_day normalized against Ultuna'
print df

norm_av_re_day_fw=np.array(df['norm_av_re_day_fw'])
norm_av_re_day_ss=np.array(df['norm_av_re_day_ss'])
norm_av_re_day_ss_fw=np.array(df['norm_av_re_day_ss_fw'])

out1=open("temp_moist_re_day_av_fw_ICBM_4p1000.txt","wb")
np.save(out1,norm_av_re_day_fw)
out1.close()
out2=open("temp_moist_re_day_av_ss_ICBM_4p1000.txt","wb")
np.save(out2,norm_av_re_day_ss_fw)
out2.close()



print '*************************'
