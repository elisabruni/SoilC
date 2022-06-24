#from https://github.com/wwieder/MIMICS/blob/master/MIMIC2_LIDET_GMDD_2015_public.R
setRepositories(addURLs =
                  c(CRANxtras = "https://ftp.igh.cnrs.fr/pub/CRAN/")) #set 1 5
library(rootSolve)
library(readxl)
library(deSolve)
library(xts)

optim_loc <- "OPTIM8"
nonoptim_loc <-"RUN_FWD8"

#OUTPUT FILES
location_dataframe<-paste0("/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/",optim_loc,"/SOC_MIMICS_optim.csv")
location_dataframe2<-paste0("/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/",optim_loc,"/optim_param_MIMICS.csv")
location_spinup_OPT<-paste0("/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/",optim_loc,"/list_SOC_spinup_opti_MIMICS.csv")
location_spinup_nonOPT<-paste0("/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/",nonoptim_loc,"/list_SOC_spinup_nonOPTI_MIMICS.csv")
location_plot<-paste0("/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/",optim_loc,"/GRAFICI_FIT/FIT_MODELLI/fit_MIMICS.pdf")

############################################################
# param_opt ######## OBJECTIVE FUNCTION
############################################################
Jparam_new<-function(param_new){
  Delta=rep(NA,exper_len)
  #j=Current_Site_Index
  
  #new param
  param_new_t<-param_new[1]
  param_new_k<-param_new[2]
  param_new_f<-param_new[3]
  #Aggiorna parametri fw
  Tpars_fw['aV']=param_new_t
  Tpars_fw['aK']=param_new_k
  Tpars_fw['fMET']=param_new_f
  #Aggiorna parametri spiup
  Tpars_ss['aV']=param_new_t
  Tpars_ss['aK']=param_new_k
  Tpars_ss['fMET']=param_new_f
  print(paste('New param',param_new))
  
  #spinup
  param_predict_c <- stode(y = Ty, time = SStime, 
                           func = XEQ_ss, parms = Tpars_ss, positive=TRUE)
  ssp = sum(param_predict_c$y)
  ssp_tC = ssp*depth/10.
  
  #print(paste("steady state optim",ssp))
  
  #Forcing forward
  #Model forcings
  state.model = param_predict_c$y
  
  #Forward
  model_out_pred <- Model.fwd(state.model,Tpars_fw,run.steps)
  
  #Annual simulated SOC
  n<-365
  yearly_model_out_pred<-aggregate(model_out_pred,list(rep(1:(nrow(model_out_pred)%/%n+1),each=n,len=nrow(model_out_pred))),mean)[-1]
  yearly_SOC_pred <- rowSums(yearly_model_out_pred[c(2:8)])*depth/10. #tC/ha
  
  #print("yearly_SOC_pred")
  #print(yearly_SOC_pred)
  
  #Add spinup
  #yearly_SOC_pred<-c(ssp_tC,yearly_SOC_pred)
  yearly_SOC_df_pred <- cbind(seq(1:exper_len),yearly_SOC_pred)
  
  #Annual observed SOC (and var)
  pred_yearly_SOC_simu_vs_obs <- cbind(yearly_SOC_df_pred,SOC_site_na)

  #selct model values where data avail
  pred_yearly_SOC_simu_vs_obs <- na.omit(pred_yearly_SOC_simu_vs_obs)
  
  #Togli NA da SOC_var
  nan_index <-attr(pred_yearly_SOC_simu_vs_obs,"na.action")[1]
  #If nan_index is not NULL
  if(!is.null(nan_index)){
    SOC_var_na<-SOC_var_na[-nan_index]
  }
  
  SOC_nonans<-pred_yearly_SOC_simu_vs_obs[,'SOC_site_na']
  count=1
  for(i in SOC_var_na){
    if(is.na(i)==TRUE | i==0.0){
      SOC<-as.numeric(SOC_nonans[count])
      #print(SOC)
      SOC_var_na[count]=rel_variance_df*SOC
    }  
    count=count+1
  }
  
  pred_yearly_SOC_simu_vs_obs=cbind(pred_yearly_SOC_simu_vs_obs,SOC_var_na)
  
  #Opt function
  Delta <- as.numeric(pred_yearly_SOC_simu_vs_obs[,'yearly_SOC_pred']) - as.numeric(pred_yearly_SOC_simu_vs_obs[,'SOC_site_na'])
  Delta2w=Delta^2/as.numeric(pred_yearly_SOC_simu_vs_obs[,'SOC_var_na'])
  Jparam_new_val = sum(Delta2w)
  
  print(paste("Jparam=",Jparam_new_val))
  
  return(Jparam_new_val)
}

#----------------analytical solutin using stode function-----------------------
XEQ_ss <- function(t, y, pars) {
  with (as.list(c(y, pars)),{
    
    #print("SPINUP")
    #print("pars")
    #print(pars)
    #print("t")
    #print(t)
    
    #Spinup
    TSOI<-forc_st
    
    #print("Vslope")
    #print(Vslope)
    #print("TSOI")
    #print(TSOI)
    
    #print(class(pars["aV"]))
    EST_LIT_in  <- forc_npp   		#gC/m2/day (Knapp et al. Science 2001)
    EST_LIT     <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/day
    I        <- array(NA, dim=2)              #Litter inputs to MET/STR 
    I[1]     <- (EST_LIT / depth) * pars["fMET"]      #partitioned to layers (mgC/cm3/day)
    I[2]     <- (EST_LIT / depth) * (1-pars["fMET"] )
    Tao_MOD1 <- sqrt(forc_npp*365/100)  #basicaily standardize against NWT (forc_npp*365=ANPP (gC/m2/y))
    tao      <- c(5.2e-4*exp(0.3*pars["fMET"] ), 2.4e-4*exp(0.1*pars["fMET"] ))	
    tao      <- tao * Tao_MOD1 *24 #daily
    fPHYS    <- c(0.3 * exp(1.3*fCLAY), 0.2 * exp(0.8*fCLAY)) 	#fraction to SOMp
    fCHEM    <- c(0.1 * exp(-3*pars["fMET"] )  , 0.3 * exp(-3*pars["fMET"] )  ) 	#fraction to SOMc
    fAVAI    <- 1- (fPHYS + fCHEM)
    #desorb   <- 9e-4 * exp(-3*(sqrt(fCLAY))) #if modified by MIC!
    #desorb   <- 3e-4 * exp(-4*(sqrt(fCLAY))) #if stand alone rate constant
    desorb   <- 1.5e-5 * exp(-1.5*(fCLAY))*24 #daily      #CHANGED FOR GLOBAL RUN!!! 
    pSCALAR  <- a * exp(-k*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
    MOD1     <- c(10, 2, 10, 3, 3, 2) 
    MOD2     <- c( 8, 2 ,4 * pSCALAR, 2, 4, 6 * pSCALAR) 	
    
    Vmax     <- exp(TSOI * Vslope + Vint) * pars["aV"] *24
    Km       <- exp(Kslope * TSOI + Kint) * pars["aK"]
    VMAX     <- Vmax * MOD1 
    #print("VMAX")
    #print(VMAX)
    KM       <- Km / MOD2
    #????????????????????
    
    #Flows to and from MIC_1
    #EQA1
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
    #EQA2
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
    #EQA4
    MICtrn[1] = MIC_1 * tao[1]  * fPHYS[1]                  #MIC_1 turnover to PHYSICAL SOM 
    MICtrn[2] = MIC_1 * tao[1]  * fCHEM[1]                  #MIC_1 turnover to CHEMICAL SOM  
    MICtrn[3] = MIC_1 * tao[1]  * fAVAI[1]                  #MIC_1 turnover to AVAILABLE SOM  
    #EQA3
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #decomp of SOMa by MIC_1
    
    #Flows to and from MIC_2
    #EQA5
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #decomp of MET litter
    #EQA6
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #decomp of SRUCTURAL litter
    
    #print('-----')
    #print("LITmin")
    #print(LITmin)
    #print("KM")
    #print(KM)
    #print("LIT_1")
    #print(LIT_1)
    #print("LIT_2")
    #print(LIT_2)
    #print('-----') 
    
    #EQA8
    MICtrn[4] = MIC_2 * tao[2]  * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM 
    MICtrn[5] = MIC_2 * tao[2]  * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM  
    MICtrn[6] = MIC_2 * tao[2]  * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM  
    
    #print('-----')
    #print("fPHYS")
    #print(fPHYS)
    #print("fCHEM")
    #print(fCHEM)
    #print("fAVAI")
    #print(fAVAI)
    
    #print("tao")
    #print(tao)
    #print("MIC_2")
    #print(MIC_2)
    #print("MICtrn")
    #print(MICtrn)
    #print('-----')
    
    #EQA7
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of SOMa by MIC_2
    
    #print('-----')
    #print("VMAX")
    #print(VMAX)
    #print("SOMmin")
    #print(SOMmin)
    #print("SOM_3")
    #print(SOM_3)
    #print('-----')
    
    #EQA9
    DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)		#desorbtion of PHYS to AVAIL (function of fCLAY)
    
    #print('-----')
    #print("DEsorb")
    #print(DEsorb)
    #print('-----')
    
    #EQA10
    OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +
                   (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))  #oxidation of C to A
    #print('-----')
    #print("OXIDAT")
    #print(OXIDAT)
    #print('-----')
    
    #can make fluxes from CHEM a function of microbial biomass size?
    
    dLIT_1 = I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    dMIC_1 = CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
    dSOM_1 = I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
    
    dLIT_2 = I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
    dMIC_2 = CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
    dSOM_2 = I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
    
    dSOM_3  = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
    
    diff_eq = c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3)
    names(diff_eq) = c("dLIT_1", "dLIT_2", "dMIC_1", "dMIC_2", "dSOM_1", "dSOM_2", "dSOM_3")
    #print(diff_eq)
    #print(length(diff_eq))
    #print(list(diff_eq))
    #print(length(list(diff_eq)))
    #pools<-c(LITmin,MICtrn,SOMmin)
    #print("pools")
    #print(pools)
    return(list(diff_eq))
  })
}

#----------------forward-----------------------

XEQ_fw <- function(t, y, pars) {
  with (as.list(c(y, pars)),{
    
    #print("FORWARD")
    #print("pars")
    #print(pars)
    
    #print("t")
    #print(t)
    
    #Forward
    TSOI<-forc_st[t]
    #print("TSOI")
    #print(TSOI)
    
    #print("length(TSOI)")
    #print(length(TSOI))
    #print("forc_st")
    #print(forc_st)
    
    #print("Vslope")
    #print(Vslope)
    #print("TSOI")
    #print(TSOI)
    #print("t")
    #print(t)
    #print(pars$aV)
    #print(pars$aK)
    #print(pars$fMET)
    #print(class(pars["aV"]))
    EST_LIT_in  <- forc_npp   		#gC/m2/day (Knapp et al. Science 2001)
    EST_LIT     <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/day
    I        <- array(NA, dim=2)              #Litter inputs to MET/STR 
    I[1]     <- (EST_LIT / depth) * pars$fMET     #partitioned to layers (mgC/cm3/day)
    I[2]     <- (EST_LIT / depth) * (1-pars$fMET )
    Tao_MOD1 <- sqrt(forc_npp*365/100)  #basicaily standardize against NWT (forc_npp*365=ANPP (gC/m2/y))
    tao      <- c(5.2e-4*exp(0.3*pars$fMET ), 2.4e-4*exp(0.1*pars$fMET ))	
    tao      <- tao * Tao_MOD1 *24 #daily
    fPHYS    <- c(0.3 * exp(1.3*fCLAY), 0.2 * exp(0.8*fCLAY)) 	#fraction to SOMp
    fCHEM    <- c(0.1 * exp(-3*pars$fMET )  , 0.3 * exp(-3*pars$fMET )  ) 	#fraction to SOMc
    fAVAI    <- 1- (fPHYS + fCHEM)
    #desorb   <- 9e-4 * exp(-3*(sqrt(fCLAY))) #if modified by MIC!
    #desorb   <- 3e-4 * exp(-4*(sqrt(fCLAY))) #if stand alone rate constant
    desorb   <- 1.5e-5 * exp(-1.5*(fCLAY))*24 #daily      #CHANGED FOR GLOBAL RUN!!! 
    pSCALAR  <- a * exp(-k*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
    MOD1     <- c(10, 2, 10, 3, 3, 2) 
    MOD2     <- c( 8, 2 ,4 * pSCALAR, 2, 4, 6 * pSCALAR) 	
    
    Vmax     <- exp(TSOI * Vslope + Vint) * pars$aV *24
    Km       <- exp(Kslope * TSOI + Kint) * pars$aK
    VMAX     <- Vmax * MOD1 
    KM       <- Km / MOD2
    
    #print(Vmax)
    ###################
    
    #Flows to and from MIC_1
    #EQA1
    LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + LIT_1)   #MIC_1 decomp of MET lit
    #EQA2
    LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + LIT_2)   #MIC_1 decomp of STRUC lit
    #EQA4
    MICtrn[1] = MIC_1 * tao[1]  * fPHYS[1]                  #MIC_1 turnover to PHYSICAL SOM 
    MICtrn[2] = MIC_1 * tao[1]  * fCHEM[1]                  #MIC_1 turnover to CHEMICAL SOM  
    MICtrn[3] = MIC_1 * tao[1]  * fAVAI[1]                  #MIC_1 turnover to AVAILABLE SOM  
    #EQA3
    SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + SOM_3)   #decomp of SOMa by MIC_1
    
    #Flows to and from MIC_2
    #EQA5
    LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + LIT_1)   #decomp of MET litter
    #EQA6
    LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + LIT_2)   #decomp of SRUCTURAL litter
    
    #print('-----')
    #print("LITmin")
    #print(LITmin)
    #print("KM")
    #print(KM)
    #print("LIT_1")
    #print(LIT_1)
    #print("LIT_2")
    #print(LIT_2)
    #print('-----') 
    
    
    #EQA8
    MICtrn[4] = MIC_2 * tao[2]  * fPHYS[2]                  #MIC_2 turnover to PHYSICAL  SOM 
    MICtrn[5] = MIC_2 * tao[2]  * fCHEM[2]                  #MIC_2 turnover to CHEMICAL  SOM  
    MICtrn[6] = MIC_2 * tao[2]  * fAVAI[2]                  #MIC_2 turnover to AVAILABLE SOM  
    
    #print('-----')
    #print("fPHYS")
    #print(fPHYS)
    #print("fCHEM")
    #print(fCHEM)
    #print("fAVAI")
    #print(fAVAI)
    
    #print("tao")
    #print(tao)
    #print("MIC_2")
    #print(MIC_2)
    #print("MICtrn")
    #print(MICtrn)
    #print('-----')
    
    
    #EQA7
    SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + SOM_3)   #decomp of SOMa by MIC_2
    
    #print('-----')
    #print("VMAX")
    #print(VMAX)
    #print("SOMmin")
    #print(SOMmin)
    #print("SOM_3")
    #print(SOM_3)
    #print('-----')
    
    #EQA9
    DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)		#desorbtion of PHYS to AVAIL (function of fCLAY)
    
    #print('-----')
    #print("DEsorb")
    #print(DEsorb)
    #print('-----')
    
    #EQA10
    OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + SOM_2)) +
                   (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + SOM_2)))  #oxidation of C to A
    #print('-----')
    #print("OXIDAT")
    #print(OXIDAT)
    #print('-----')
    
    #can make fluxes from CHEM a function of microbial biomass size?
    
    dLIT_1 = I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
    dMIC_1 = CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
    dSOM_1 = I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
    
    dLIT_2 = I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
    dMIC_2 = CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
    dSOM_2 = I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
    
    dSOM_3  = MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]
    
    diff_eq = c(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3)
    names(diff_eq) = c("dLIT_1", "dLIT_2", "dMIC_1", "dMIC_2", "dSOM_1", "dSOM_2", "dSOM_3")
    #print(diff_eq)
    #print(length(diff_eq))
    #print(list(diff_eq))
    #print(length(list(diff_eq)))
    #pools<-c(LITmin,MICtrn,SOMmin)
    #print("pools")
    #print(pools)
    
    fluxes<-c(LITmin,KM, LIT_1,LIT_2,MICtrn,MIC_2,tao,fAVAI,fCHEM,fPHYS,VMAX,SOMmin,SOM_3,DEsorb,OXIDAT)
    #print("c(LITmin,KM, LIT_1,LIT_2,MICtrn,MIC_2,tao,fAVAI,fCHEM,fPHYS,VMAX,SOMmin,SOM_3,DEsorb,OXIDAT)")
    #print(fluxes)
    if(anyNA(fluxes)){
      stop("Il y a des NA")
    } 
    return(list(diff_eq))
  })
}

#-----Solve forward-------
Model.fwd <- function (state, parameters, times=run.steps) {
  output <- ode(y = state, t=run.steps, func=XEQ_fw, parms = parameters, method="rk4") #solve ode, return output
  return(as.data.frame(cbind(time = output[run.steps.minus.one,"time"], 
                             LIT_1 = output[run.steps.minus.one,"LIT_1"],LIT_2 = output[run.steps.minus.one,"LIT_2"], 
                             MIC_1 = output[run.steps.minus.one,"MIC_1"], MIC_2 = output[run.steps.minus.one,"MIC_2"],
                             SOM_1 = output[run.steps.minus.one,"SOM_1"],SOM_2 = output[run.steps.minus.one,"SOM_2"], SOM_3 = output[run.steps.minus.one,"SOM_3"])))
}
#---------------------------------------------------------


#Define input directory
inputdir <- "/Users/ebruni/Desktop/DOTTORATO/"
loc_exp = paste0(inputdir,'SCRIPT_MODELLI/MULTIMODEL/SITES_dataset_multimodel.xlsx')

#Upload excel data
C_input_exp = read_excel(loc_exp, range="A1:AK1539")
C_input_exp = C_input_exp[-c(1,2),]
site_names_all = unique(C_input_exp$ID.Site)
site_names_all=site_names_all[!is.na(site_names_all)]
N_sites_all=length(site_names_all)

#Control plots

site_T0_array_all=c('CHNO3_Min', 'COL_T0', 'CREC3_Min', 'FEU_T0', 'LAJA2_Min', 
                    'RHEU1_Min', 'RHEU2_T0','ARAZ_D0_N0', 'ULT_P0_B', 'BROAD_3_Nill', 
                    'TREV1_Min','AVRI_T1TR','BOLO_T0', 'GRAB_CP','MUNCHE_CP','RITZ_CP')

SStime<-5000*365 #daily

list_SOC_0 <-rep(NA,length(site_names_all))

list_forc_npp <-rep(NA,length(site_names_all))
list_clay <-rep(NA,length(site_names_all))
list_pH <-rep(NA,length(site_names_all))
list_exper_len <-rep(NA,length(site_names_all))

list_date_init_ss <- rep(NA,length(site_names_all)) 
list_date_end_ss <- rep(NA,length(site_names_all))
list_date_init <- rep(NA,length(site_names_all))
list_date_end <-  rep(NA,length(site_names_all))

list_temp_ss <- list()
list_water_ss <- list()
list_temp_fw <- list()
list_water_fw <- list()

SOC_list <-list()
SOC_var_list<-list()

TSOI_fw <-list()

j<-1
for(site in site_names_all){
  #print(site)
  #j<-1
  #site='CHNO3'
  site_df <- subset(C_input_exp, C_input_exp$ID.Site == site)
  site_T0_name<-site_T0_array_all[j]
  site_T0 <- subset(site_df,site_df$ID.Treatment==site_T0_name)
  #print(site_T0)
  #Choose years of simulations
  #site='RHEU1'
  
  #date_init_ss <- 1980
  #date_end_ss <- 2010
  
  year_site <-subset(site_T0$Year, site_T0$ID.Site == site)
  exper_len <- length(year_site)
  year_0 <- min(year_site)
  year_end <- max(year_site)
  date_init_ss <- as.character(year_0-30)
  date_end_ss <- as.character(year_0-1)
  
  #print(date_init_ss)
  #print(date_end_ss)
  
  #Read data spinup
  soil_temp_ss <- paste0(inputdir,'SCRIPT_MODELLI/',site,'/',"temp_",site,"_",date_init_ss,"_",date_end_ss,".txt")
  forc_st_table <- read.table(soil_temp_ss) #degreeC
  
  above <-mean(as.numeric(subset(site_T0$ABOVE, site_T0$ID.Site == site)), na.rm=TRUE)
  below<-mean(as.numeric(subset(site_T0$BELOW, site_T0$ID.Site == site)), na.rm=TRUE)
  forc_npp <- (above+below)*100/365 #gC/m2/day
  pH <-mean(subset(site_T0$pH, site_T0$ID.Site == site))
  clay <-mean(as.numeric(subset(site_T0$Clay, site_T0$ID.Site == site)))
  
  #Forcing forward
  #Read data forward
  date_init <- year_0
  date_end <- year_end
  soil_temp <- paste0(inputdir,'SCRIPT_MODELLI/',site,'/',"temp_",site,"_",date_init,"_",date_end,".txt")
  forc_st_table_fw <- read.table(soil_temp) #degreeC
 
  ######################################################
  #Convert daily temperature to annual average temperature (Forward)
  ######################################################
  #Convert dataframe to list
  forc_st_table_fw_series <- as.list(forc_st_table_fw)$V1
  #Select series values for the length of experiment
  forc_st_table_fw_series <-forc_st_table_fw_series[1:(exper_len*365)]
  #Add dates to values
  forc_st_table_fw_xts <- xts(forc_st_table_fw_series,as.Date(1:length(forc_st_table_fw_series),origin=paste0(date_init,"-01-01")))
  #Calculate the annual mean
  forc_st_table_fw_yearly <- apply.yearly(forc_st_table_fw_xts,mean)
  #Convert xts to matrix
  forc_st_table_fw_yearly_ma <- as.matrix(forc_st_table_fw_yearly)
  #Convert matrix to list of list
  forc_st_table_fw_yearly_lol<-as.list(forc_st_table_fw_yearly_ma)
  #Convert list of list to list
  forc_st_table_fw_yearly_list <- do.call(what = "rbind", forc_st_table_fw_yearly_lol)[,1]
  
  TSOI_fw <- append(TSOI_fw,list(forc_st_table_fw_yearly_list))
  ############################################################
  
  SOC_site_na <-subset(site_T0$SOC, site_T0$ID.Site == site)
  SOC_var_na  <-subset(site_T0$`SOC variance`, site_T0$ID.Site == site)
  
  SOC_site <- na.omit(SOC_site_na)
  SOC_0 <- SOC_site[1]
  print(paste('Observed initial SOC (tC/ha):',SOC_0))
  
  SOC_list<-append(SOC_list,list(SOC_site_na))
  SOC_var_list<-append(SOC_var_list,list(SOC_var_na))
  
  #Salva dati
  list_SOC_0[j]<-SOC_0
  
  list_forc_npp[j] <-forc_npp
  list_clay[j] <- clay
  list_pH[j] <- pH
  list_exper_len[j] <-exper_len
  
  list_date_init_ss[j] <-date_init_ss
  list_date_end_ss[j] <-date_end_ss
  list_date_init[j] <-date_init
  list_date_end[j] <-date_end
  
  j<-j+1
  
}

#print(SOC_list)
#print(SOC_var_list)

#Calculate percentage variance among the dataset
jnew<-1
relative_variances<-rep(NA,length(site_names_all))
for(i in SOC_var_list){
  
  SOC<-SOC_list[jnew]
  SOC<-as.numeric(SOC[[1]])
  SOC<-replace(SOC, is.na(SOC), 0) #togli nan
  i<-replace(i, is.na(i), 0) #togli nan
  var_rel<-i/SOC
  mean_var_rel<-mean(var_rel,na.rm = TRUE)
  relative_variances[jnew]<-mean_var_rel
  jnew<-jnew+1
}

#print(relative_variances)

rel_variance_df<-mean(relative_variances[relative_variances>0])


##################
#import metab ratios (mean of AB, BE  optimized Century metabolic fractions)
#################
#inp_MET_frac=paste0(inputdir,'SCRIPT_MODELLI/MIMICS/MET_frac.csv')
#MET_frac<-read.table(inp_MET_frac,header=TRUE,sep=",",fill=TRUE,na.strings=c(""," ","NA"))$X0
#WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
fMET_prior <- 0.6916

depth=30. #cm
#-----------------caclulate parameters---------------------------
#Calculate Vmax & (using parameters from German 2012, as in Wieder et al. 2013 Nature Climate Change)
Vslope   <- array(0.063,dim=6) #daily
Vint     <- 5.47 #daily
aV       <- 8e-6

Kslope   <- array(NA,dim=6)
Kslope[1]<- 0.017 #META LIT to MIC_1
Kslope[2]<- 0.027 #STRU LIT to MIC_1 
Kslope[3]<- 0.017 #AVAI SOM to MIC_1 
Kslope[4]<- 0.017 #META LIT to MIC_2
Kslope[5]<- 0.027 #STRU LIT to MIC_2
Kslope[6]<- 0.017 #AVAI SOM to MIC_2
Kint     <- 3.19
aK       <- 10

CUE        <- c(0.55, 0.25, 0.75, 0.35)  #for LITm and LITs entering MICr and MICK, respectively

#------NEW Parameters--------------
k        <- 2.0    #2.0			#REDUCED FROM 3 TO 1, REDUCES TEXTURE EFFECTS ON SOMa decay
a        <- 2.0    #2.2			#increased from 4.0 to 4.5

cMAX     <- 1.4                    #ORIG 1.4 Maximum CHEM SOM scalar w/   0% Clay 
cMIN     <- 1.2                    #ORIG 1.4 Minimum CHEM SOM scalar w/ 100% Clay 
cSLOPE   <- cMIN - cMAX            #Slope of linear function of cSCALAR for CHEM SOM  

#------------!!MODIFIERS AS IN MIMICS2_b!!---------------

KO       <- c(4,4)      #scalar modifies Km of Oxidat	
FI       <- c(0.05, 0.05)

LITmin  <- rep(NA, dim=4)
MICtrn  <- rep(NA, dim=6)
SOMmin  <- rep(NA, dim=2)
DEsorb  <- rep(NA, dim=1)
OXIDAT  <- rep(NA, dim=1)

#------------#------------#------------

list_SOC_sdvers<-list()
list_SS <- rep(NA,length(site_names_all))
list_SOC_obs <-list()
list_spinup_nonOPT <- list()
j<-1
for(site in site_names_all){
  
  #site='CHNO3'
  #j=1
  #Initialize
  forc_npp <-list_forc_npp[j] #gC/m2/day
  
  fCLAY <-list_clay[j]
  pH <-list_pH[j]
  exper_len<-list_exper_len[j]
  fMET<-fMET_prior
  EST_LIT_in  <- forc_npp   		# ( g/m2/day, Knapp et al. Science 2001)
  EST_LIT     <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h 
  I        <- array(NA, dim=2)              #Litter inputs to MET/STR 
  I[1]     <- (EST_LIT / depth) * fMET      #partitioned to layers (mgC/cm3/h)
  I[2]     <- (EST_LIT / depth) * (1-fMET)
  
  #Esperiment length, temperature and moisture
  date_init_ss<-list_date_init_ss[j]
  date_end_ss<-list_date_end_ss[j]
  date_init<-list_date_init[j]
  date_end<-list_date_end[j]
  
  #TOGLI
  #Tao_MOD1 <- sqrt(forc_npp*365/100)  #basicaily standardize against NWT (forc_npp*365=ANPP (gC/m2/y))
  #tao      <- c(5.2e-4*exp(0.3*fMET), 2.4e-4*exp(0.1*fMET))	
  #tao      <- tao * Tao_MOD1*24 #daily
  #fPHYS    <- c(0.3 * exp(1.3*fCLAY), 0.2 * exp(0.8*fCLAY)) 	#fraction to SOMp
  #fCHEM    <- c(0.1 * exp(-3*fMET)  , 0.3 * exp(-3*fMET)  ) 	#fraction to SOMc
  #fAVAI    <- 1- (fPHYS + fCHEM)
  #desorb   <- 9e-4 * exp(-3*(sqrt(fCLAY))) #if modified by MIC!
  #desorb   <- 3e-4 * exp(-4*(sqrt(fCLAY))) #if stand alone rate constant
  #desorb   <- 1.5e-5 * exp(-1.5*(fCLAY))*24 #daily      #CHANGED FOR GLOBAL RUN!!! 
  #pSCALAR  <- a * exp(-k*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  #MOD1     <- c(10, 2, 10, 3, 3, 2) 
  #MOD2     <- c( 8, 2 ,4 * pSCALAR, 2, 4, 6 * pSCALAR) 	
  #--
  
  #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  #print(paste("exper_len",exper_len))
  
  #print(paste("len forc st fw", length(forc_st_fw[[1]])/365))
  
  SOC_0<-list_SOC_0[j]
  SOC_site_na<-SOC_list[[j]]
  SOC_var_na<-SOC_var_list[[j]]
  
  #run.steps <- seq(1,(exper_len-1)*365)
  run.steps <- seq(1,exper_len*365)
  run.steps.minus.one <- run.steps[1:length(run.steps)-1]
  
  #Forcing spinup
  forc_st_ss<-mean(as.numeric(TSOI_fw[j][[1]])) #Average temperature among different experiment's years
  
  #Forcing forward  
  #Repeat each year 365 times
  forc_st_fw<-rep(TSOI_fw[j][[1]], each=365)
  #Remove temperature of first year  (because already simulated by spinup)
  #forc_st_fw<-list(tail(forc_st_fw,-365))
  forc_st_fw<-list(forc_st_fw)
  #print(length(forc_st_fw[[1]]))
  
  #initialize pools
  LIT     <- I   # * 1e3
  MIC     <- I   # * 1e2
  SOM     <- rep(NA, 3) 
  SOM[1]  <- I[1]
  SOM[2]  <- I[2]
  SOM[3]  <- I[1] 
  
  #Spinup 
  Tpars_ss <- c( forc_npp=forc_npp,
                 Kslope=Kslope,Kint=Kint,Vslope=Vslope,Vint=Vint,aV=aV,aK=aK,#Added
                 forc_st=forc_st_ss,fMET=fMET) #mean annual temperature 
  #print(Tpars_ss)
  Ty    <- c( LIT_1 = LIT[1], LIT_2 = LIT[2], 
              MIC_1 = MIC[1], MIC_2 = MIC[2], 
              SOM_1 = SOM[1], SOM_2 = SOM[2], SOM_3 = SOM[3] )
  
  
  SS_prova  <- stode(y = Ty, time = SStime, fun = XEQ_ss, parms = Tpars_ss, positive = TRUE)
  SS_prova[[1]]
  
  #Model forcing
  soultions_SS<-SS_prova$y #(?) mgC/cm3
  #Save spinup
  list_spinup_nonOPT<-append(list_spinup_nonOPT,list(soultions_SS))
  
  soultions_SS_tC<-soultions_SS*depth/10 #convert mgC/cm3 to gC/m3 (*1000), then gC/m3 to tC/m3 (1/1e+06), then tC/m3 to tC/m2 (*depth/100 (m)), then tC/m2 to tC/ha (*1e+04)
  
  #------------------------------------
  tot_SS <- sum(soultions_SS_tC)
  print(paste(site," total C spinup (tC/ha):",tot_SS))
  
  list_SS[j]=tot_SS
  
  #Forward
  Tpars_fw <- c( forc_npp=forc_npp,
                 Kslope=Kslope,Kint=Kint,Vslope=Vslope,Vint=Vint,aV=aV,aK=aK,#Added
                 forc_st=forc_st_fw,fMET=fMET) #mean annual temperature 
  

  #forward simu
  model_out <- Model.fwd(soultions_SS,Tpars_fw,run.steps)
  
  print(model_out)
  #Annual simulated SOC
  n<-365
  yearly_model_out<-aggregate(model_out,list(rep(1:(nrow(model_out)%/%n+1),each=n,len=nrow(model_out))),mean)[-1]
  
  #print(yearly_model_out)
  yearly_SOC <- rowSums(yearly_model_out[c(2:8)])*depth/10 # from mgC/cm3 to tC/ha (*depth/10)
  
  #Add spinup
  #yearly_SOC<-c(tot_SS,yearly_SOC)
  yearly_SOC_df <- cbind(seq(1:exper_len),yearly_SOC)
  
  #Annual observed SOC
  yearly_SOC_simu_vs_obs <- cbind(yearly_SOC_df,SOC_site_na)
  print("Simulated vs observed SOC (tC/ha)")
  print(yearly_SOC_simu_vs_obs)
  
  #SOC stocks evolution model (list of sites)
  list_SOC_sdvers<-append(list_SOC_sdvers,list(yearly_SOC))
  #SOC stocks evolution obs (list of sites)
  list_SOC_obs <-append(list_SOC_obs,list(SOC_site_na))
  
  
  j=j+1
}

#print(list_SS)
#print(list_SOC_0)
SOCSim<-cbind(site_names_all,list_SOC_0,list_SS)
print(SOCSim)

##################
#OPTIMIZATION
##################
list_opt_param<-list()
list_SOC_simu <- list()
list_spinup_OPT <- list()
j<-1
for (site in site_names_all){
  #site='CHNO3'
  #j=1
  print("******************")
  print("******************")
  print(paste("SITE:",site))
  print("******************")
  print("******************")
  #Initialize
  forc_npp <-list_forc_npp[j] #gC/m2/day
  fCLAY <-list_clay[j]
  pH <-list_pH[j]
  exper_len<-list_exper_len[j]
  fMET<-fMET_prior
  #Esperiment length, temperature and moisture
  date_init_ss<-list_date_init_ss[j]
  date_end_ss<-list_date_end_ss[j]
  date_init<-list_date_init[j]
  date_end<-list_date_end[j]
  #WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
  SOC_0<-list_SOC_0[j]
  SOC_site_na<-SOC_list[[j]]
  SOC_var_na<-SOC_var_list[[j]]
  
  SOC_simu_sd <-list_SOC_sdvers[[j]]
  SOC_simu_sd_obs <- cbind(SOC_site_na,SOC_simu_sd)
  SOC_simu_sd_obs <- na.omit(SOC_simu_sd_obs)
  
  #run.steps <- seq(1,(exper_len-1)*365) #daily
  run.steps <- seq(1,exper_len*365)
  run.steps.minus.one <- run.steps[1:length(run.steps)-1]
  
  #Forcing spinup
  forc_st_ss<-mean(as.numeric(TSOI_fw[j][[1]])) #Average temperature among different experiment's years
  
  #RICALCOLA I x inizializzare i pool
  EST_LIT_in  <- forc_npp   		# ( g/m2/day, Knapp et al. Science 2001)
  EST_LIT     <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h 
  I        <- array(NA, dim=2)              #Litter inputs to MET/STR 
  I[1]     <- (EST_LIT / depth) * fMET      #partitioned to layers (mgC/cm3/h)
  I[2]     <- (EST_LIT / depth) * (1-fMET)
  
  #Spinup 
  #initialize pools
  LIT     <- I   # * 1e3
  MIC     <- I   # * 1e2
  SOM     <- rep(NA, 3) 
  SOM[1]  <- I[1]
  SOM[2]  <- I[2]
  SOM[3]  <- I[1]
  
  Ty    <- c( LIT_1 = LIT[1], LIT_2 = LIT[2], 
              MIC_1 = MIC[1], MIC_2 = MIC[2], 
              SOM_1 = SOM[1], SOM_2 = SOM[2], SOM_3 = SOM[3] )
  
  Tpars_ss <- c( forc_npp=forc_npp,
                 Kslope=Kslope,Kint=Kint,Vslope=Vslope,Vint=Vint,aV=aV,aK=aK,#Added
                 forc_st=forc_st_ss,fMET=fMET) #mean annual temperature 
  
  #Forcing forward  
  #Repeat each year 365 times
  forc_st_fw<-rep(TSOI_fw[j][[1]], each=365)
  #Remove temperature of first year  (because already simulated by spinup)
  #forc_st_fw<-list(tail(forc_st_fw,-365))
  forc_st_fw<-list(forc_st_fw)
  
  Tpars_fw <- c( forc_npp = forc_npp,
                 Kslope=Kslope,Kint=Kint,Vslope=Vslope,Vint=Vint,aV=aV,aK=aK,#Added
                 forc_st=forc_st_fw,fMET=fMET) #mean annual temperature 
  
  #OPTIM ALGO
  prior_param1<- Tpars_fw$aV
  prior_param2 <- Tpars_fw$aK
  prior_param3 <- Tpars_fw$fMET
  prior_params<-c(prior_param1,prior_param2,prior_param3)

  #For fMET bounds should be 0.013 - 0.85, but MIC1 pool at steady state too low if fMET=0.013 in Broabalk-> lower bound set to 0.15
  optim_func<-optim(prior_params,Jparam_new, method = "L-BFGS-B",
                    lower = c(prior_param1-0.000006,2,0.15), upper = c(prior_param1,20.,0.85))

  print(optim_func)
  
  print("Optim param")
  print(optim_func$par)
  print("Optim Jparam")
  print(optim_func$value)
  print("Convergence")
  print(optim_func$convergence)
  
  list_opt_param<-append(list_opt_param,list(optim_func$par))

  ####################################################################
  #OPTIMIZATION OVER, RICALCULATE SPINUP+FORWARD WITH OPTIMIZED PARAMS
  ####################################################################
  #Optimized version 
  #Lista parametri ss
  param_list_opt_ss <- Tpars_ss
  new_params<-optim_func$par
  param_list_opt_ss['aV'] <- new_params[1]
  param_list_opt_ss['aK'] <- new_params[2]
  param_list_opt_ss['fMET'] <- new_params[3]
  
  #RICALCOLA I x inizializzare i pool con fMET_optimized
  EST_LIT_in  <- forc_npp   		# ( g/m2/day, Knapp et al. Science 2001)
  EST_LIT     <- EST_LIT_in  * 1e3 / 1e4    #mgC/cm2/h 
  I        <- array(NA, dim=2)              #Litter inputs to MET/STR 
  I[1]     <- (EST_LIT / depth) * param_list_opt_ss['fMET']      #partitioned to layers (mgC/cm3/h)
  I[2]     <- (EST_LIT / depth) * (1-param_list_opt_ss['fMET'])
  
  #Spinup 
  #initialize pools
  LIT     <- I   # * 1e3
  MIC     <- I   # * 1e2
  SOM     <- rep(NA, 3) 
  SOM[1]  <- I[1]
  SOM[2]  <- I[2]
  SOM[3]  <- I[1]
  
  Ty    <- c( LIT_1 = LIT[1], LIT_2 = LIT[2], 
              MIC_1 = MIC[1], MIC_2 = MIC[2], 
              SOM_1 = SOM[1], SOM_2 = SOM[2], SOM_3 = SOM[3] )
  
  opt_SS <- stode(y = Ty, time = SStime, fun = XEQ_ss, parms = param_list_opt_ss, positive = TRUE)
  opt_ssp = sum(opt_SS$y)
  opt_ssp_tC = opt_ssp*depth/10
  
  #Forward
  #Lista parametri fw
  param_list_opt_fw <- Tpars_fw
  
  param_list_opt_fw['aV'] <- new_params[1]
  param_list_opt_fw['aK'] <- new_params[2]
  param_list_opt_fw['fMET'] <- new_params[3]
  
  print("Optimized param_list")
  print(param_list_opt_fw)  
  
  #Model forcings
  state.model_opt = opt_SS$y
  #Save spinup
  list_spinup_OPT<-append(list_spinup_OPT,list(state.model_opt))

  model_out_opt <- Model.fwd(state.model_opt,param_list_opt_fw,run.steps)
  
  #Annual simulated SOC
  n<-365
  yearly_model_out_opt<-aggregate(model_out_opt,list(rep(1:(nrow(model_out_opt)%/%n+1),each=n,len=nrow(model_out_opt))),mean)[-1]
  yearly_SOC_opt <- rowSums(yearly_model_out_opt[c(2:8)])*depth/10. # from mgC/cm3 to tC/ha (*depth/10)
  
  #Add spinup
  #yearly_SOC_opt<-c(opt_ssp_tC,yearly_SOC_opt)
  yearly_SOC_df_opt <- cbind(seq(1:exper_len),yearly_SOC_opt)
  
  #Annual observed SOC (and var)
  opt_yearly_SOC_simu_vs_obs <- cbind(yearly_SOC_df_opt,SOC_site_na)
  #selct model values where data avail
  opt_yearly_SOC_simu_vs_obs <- na.omit(opt_yearly_SOC_simu_vs_obs)
  
  #print("Standard version")
  #print(SOC_simu_sd_obs)
  
  sd_opt_SOC <- cbind("SOC_obs"=SOC_simu_sd_obs[,'SOC_site_na'],"SOC_simu_sd"=SOC_simu_sd_obs[,'SOC_simu_sd'],"SOC_simu_opt"=opt_yearly_SOC_simu_vs_obs[,'yearly_SOC_opt'])
  #sd_opt_SOC <- cbind("SOC_obs"=SOC_simu_sd_obs[,'SOC_site_na'],"SOC_simu_opt"=opt_yearly_SOC_simu_vs_obs[,'yearly_SOC_opt'])
  print(opt_yearly_SOC_simu_vs_obs)
  print(sd_opt_SOC)
  list_SOC_simu<-append(list_SOC_simu,list(sd_opt_SOC))
  
  j<-j+1
}


print(list_opt_param)
print(list_SOC_simu)



#Save optimized SOC in csv file

count_df<-1
for(i in list_SOC_simu){
  site<-site_names_all[count_df]
  i<-cbind(site,i)
  print(i)
  if(count_df==1){
    colnames<-colnames(i)
    cat(c("",colnames), file = location_dataframe,sep=",")
    cat("\n", file = location_dataframe, append = TRUE)
  }
  #cat(site, file = location_dataframe,append=TRUE)
  #cat("\n", file = location_dataframe, append = TRUE)
  write.table(i, location_dataframe, col.names=FALSE, sep=",", append=TRUE)
  count_df<-count_df+1
}


#Save optimized param in csv file

count_df<-1
for(i in list_opt_param){
  site<-site_names_all[count_df]
  if(count_df==1){
    colnames<-colnames(i)
    cat(c("",colnames), file = location_dataframe2,sep=",")
    cat("\n", file = location_dataframe2, append = TRUE)
  }
  cat(site, file = location_dataframe2,append=TRUE)
  cat("\n", file = location_dataframe2, append = TRUE)
  write.table(i, location_dataframe2, col.names=FALSE, sep=",", append=TRUE)
  count_df<-count_df+1
}


#SPINUP SAVE
#Save spinup in CSV file
#OPT
count_df<-1
for(i in list_spinup_OPT){
  print(i)
  site<-site_names_all[count_df]
  #i<-cbind(site,i)
  i<-append(site,i)
  print(i)
  trans_i <-t(i)
  if(count_df==1){
    colnames<-colnames(trans_i)
    cat(c("",colnames), file = location_spinup_OPT,sep=",")
    cat("\n", file = location_spinup_OPT, append = TRUE)
  }
  
  write.table(trans_i, location_spinup_OPT, col.names=FALSE, sep=",", append=TRUE)
  count_df<-count_df+1
}

#nonOPTI
count_df<-1
for(i in list_spinup_nonOPT){
  print(i)
  site<-site_names_all[count_df]
  #i<-cbind(site,i)
  i<-append(site,i)
  print(i)
  trans_i <-t(i)
  if(count_df==1){
    colnames<-colnames(trans_i)
    cat(c("",colnames), file = location_spinup_nonOPT,sep=",")
    cat("\n", file = location_spinup_nonOPT, append = TRUE)
  }
  
  write.table(trans_i, location_spinup_nonOPT, col.names=FALSE, sep=",", append=TRUE)
  count_df<-count_df+1
}

#--------PLOT-----------
#Print SOC fits

pdf(file = , location_plot)
for(i in c(1:length(site_names_all))){
  print("site")
  #i=1
  print(site_names_all[i])
  #SOC_fit <- cbind(list_SOC_sdvers[[i]],list_SOC_obs[[i]])
  #years<-c(1:length(list_SOC_sdvers[[i]]))
  SOC_fit <-list_SOC_simu[[i]]
  print(SOC_fit)
  SOC_obs<-SOC_fit[,1]
  SOC_sd<-SOC_fit[,2]
  SOC_opt<-SOC_fit[,3]
  years<-c(1:length(SOC_obs))
  min_y<-min(min(as.numeric(SOC_obs)),min(as.numeric(SOC_sd)),min(as.numeric(SOC_opt)),na.rm = TRUE)
  #print(min_y)
  max_y<-max(max(as.numeric(SOC_obs)),max(as.numeric(SOC_sd)),max(as.numeric(SOC_opt)),na.rm = TRUE)
  #print(max_y)
  
  plot(years,SOC_obs,col='black',
       ylab="SOC stocks (tC/ha/year)",xlab = "Years",main=site_names_all[i],
       ylim = c(as.numeric(min_y),as.numeric(max_y)))
  lines(years,SOC_sd,col='blue',type="l")
  lines(years,SOC_opt,col='red',type="l")
  legend("topright",
         c("Standard simu","Optimized simu"),
         fill=c("blue","red"))
  
}
dev.off()
