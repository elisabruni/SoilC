#Install libraries
setRepositories(addURLs =
                  c(CRANxtras = "https://ftp.igh.cnrs.fr/pub/CRAN/"))
#set to 1 /slash/ 5
#install.packages("desolve")
#install.packages("readxl")
library(readxl)
library(rootSolve)
library(deSolve)
#Define input directory
inputdir <- "/Users/ebruni/Desktop/DOTTORATO/"
loc_exp = paste0(inputdir,'SCRIPT_MODELLI/MULTIMODEL/SITES_dataset_multimodel.xlsx')

#OUTPUT FILES
location_dataframe<-"/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM5/SOC_MILLENNIALv2_optim.csv"
location_dataframe2<-"/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM5/optim_param_MILLENNIALv2.csv"
location_spinup_OPT<-"/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM5/list_SOC_spinup_opti_MILv2.csv"
location_spinup_nonOPT<-"/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/RUN_FWD5/list_SOC_spinup_nonOPTI_MILv2.csv"
location_plot<-"/Users/ebruni/Desktop/DOTTORATO/SCRIPT_MODELLI/MULTIMODEL/OPTIM_MODELS/OPTIM5/GRAFICI_FIT/FIT_MODELLI/fit_MILLENNIAL_v2.pdf"

#Upload excel data
C_input_exp = read_excel(loc_exp, range="A1:AK1539")
C_input_exp = C_input_exp[-c(1,2),]
site_names_all = unique(C_input_exp$ID.Site)
site_names_all=site_names_all[!is.na(site_names_all)]
N_sites_all=length(site_names_all)

#Control plots
site_T0_array_all=c('CHNO3_Min', 'COL_T0', 'CREC3_Min', 'FEU_T0', 'JEU2_M0', 'LAJA2_Min', 
                    'RHEU1_Min', 'RHEU2_T0','ARAZ_D0_N0', 'ULT_P0_B', 'BROAD_3_Nill', 
                    'TREV1_Min','AVRI_T1TR','BOLO_T0', 'GRAB_CP','MUNCHE_CP','RITZ_CP')
#Parameters

#Read in parameters
param_list<-c(param_pi=0.66, param_pa=0.33, kaff_pl=1e4, alpha_pl=2.5e12,
              eact_pl=64320.,rate_pa=0.02, rate_break=0.019, rate_leach=0.0015, 
              kaff_des=1, param_p1=0.186, param_p2=0.216, kaff_lb=290., alpha_lb=2.6e12,
              eact_lb=60260.,rate_bd=0.0036, rate_ma=0.02, cue_ref=0.6,
              cue_t=0.012, tae_ref=15., matpot=15., lambda=2.1e-4, 
              porosity = 0.6, kamin=0.2, param_pb=0.5, param_c1=0.86)

#Initial state
state <- c(POM = 1, LMWC = 1, AGG = 1, MIC = 1, MAOM=1)
#Set up with daily steps
SStime<-5000*365


############################################################
# param_opt ######## OBJECTIVE FUNCTION
############################################################
Jparam_new<-function(param_new){
  Delta=rep(NA,exper_len)
  #j=Current_Site_Index
  
  #new param
  param_new_t<-param_new[1]
  param_new_k<-param_new[2]
  param_new_m<-param_new[3]
  param_list['eact_lb']=param_new_t
  param_list['eact_pl']=param_new_k
  param_list['kaff_pl']=param_new_m
  print(paste('New param',param_new))
  
  #spinup
  
  param_predict_c <- stode(y = state, time = SStime, 
                           func = decomp, parms = param_list, temp_func=forc_st_ss,water_func=forc_sw_ss, positive=TRUE,rtol = 1e-4)
  ssp = sum(param_predict_c$y)
  ssp_tC = ssp/100.
  
  #print(paste("steady state optim",ssp))
  
  #Forward
  #Model forcings
  state.model = param_predict_c$y
  
  model_out_pred <- Model.fwd(param_list,run.steps)
  #Annual simulated SOC
  n<-365
  yearly_model_out_pred<-aggregate(model_out_pred,list(rep(1:(nrow(model_out_pred)%/%n+1),each=n,len=nrow(model_out_pred))),mean)[-1]
  yearly_SOC_pred <- rowSums(yearly_model_out_pred[c(2:6)])/100 #tC/ha
  
  #print("yearly_SOC_pred")
  #print(yearly_SOC_pred)
  
  #Add spinup
  yearly_SOC_pred<-c(ssp_tC,yearly_SOC_pred)
  yearly_SOC_df_pred <- cbind(seq(1:exper_len),yearly_SOC_pred)
  
  #print("yearly_SOC_df_pred")
  #print(yearly_SOC_df_pred)
  
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

#------------
#Decomp function
decomp <- function(step.num,state,parameters,temp_func,water_func) {
  with(as.list(c(state,parameters)), {
    
    forc_st_prova=temp_func[[1]]
    forc_sw_prova=water_func[[1]]
    
    #print('forc_sw_prova')
    #print(forc_sw_prova)    
    # Soil type properties
    #Equation 11
    #Mayes 2012, SSAJ
    kaff_lm = exp(-parameters["param_p1"] * 7 - parameters["param_p2"]) * parameters["kaff_des"]
    
    #print("kaff_lm")
    #print(kaff_lm)
    #Equation 12
    #Georgiou in review
    param_qmax = parameters["param_bulkd"] * parameters["param_c1"] * parameters["param_claysilt"]
    
    #print("param_qmax")
    #print(param_qmax)
    
    # Hydrological properties
    
    #Equation 5
    scalar_wd = (forc_sw_prova[step.num]/ parameters["porosity"])^0.5
    
    #print("scalar_wd")
    #print(scalar_wd)
    
    #Equation 4
    scalar_wb = exp(parameters["lambda"] * -parameters["matpot"]) * (parameters["kamin"] + (1 - parameters["kamin"]) * ((parameters["porosity"] - forc_sw_prova[step.num]) / parameters["porosity"])^0.5) * scalar_wd
    
    #print("scalar_wb")
    #print(scalar_wb)
    
    # Decomposition
    
    gas_const <- 8.31446
    
    #Equation 3
    vmax_pl = parameters["alpha_pl"] * exp(-parameters["eact_pl"] / (gas_const * (forc_st_prova[step.num] + 273.15)))
    
    #print("vmax_pl")
    #print(vmax_pl)
    
    #Equation 2
    # POM -> LMWC
    if(POM>0 && MIC>0){
      f_PO_LM = vmax_pl * scalar_wd * POM * MIC / (parameters["kaff_pl"] + MIC)
    }else{
      f_PO_LM=0
    }
    #Equation 6
    # POM -> AGG
    if(POM>0){
      f_PO_AG = parameters["rate_pa"] * scalar_wd * POM
    }else{
      f_PO_AG=0
    }
    
    #Equation 7
    # AGG -> MAOM
    if(AGG>0){
      f_AG_break = parameters["rate_break"] * scalar_wd * AGG
    }else{
      f_AG_break=0
    }
    
    
    #Equation 9
    # LMWC -> out of system leaching
    if(LMWC>0){
      f_LM_leach = parameters["rate_leach"] * scalar_wd * LMWC
    }else{
      f_LM_leach=0
    }
    
    #Equation 10
    # LMWC -> MAOM
    if(LMWC>0 && MAOM>0){
      f_LM_MA = scalar_wd * kaff_lm * LMWC * (1 - MAOM / param_qmax)
    }else{
      f_LM_MA=0
    }
    
    #Equation 13
    # MAOM -> LMWC
    if(MAOM>0){
      f_MA_LM = parameters["kaff_des"] * MAOM / param_qmax
    }else{
      f_MA_LM=0
      
    }
    
    #Equation 15
    vmax_lb = parameters["alpha_lb"] * exp(-parameters["eact_lb"] / (gas_const * (forc_st_prova[step.num] + 273.15)))
    
    #print("vmax_lb")
    #print(vmax_lb)
    
    #Equation 14
    # LMWC -> MIC
    if(LMWC>0 && MIC>0){
      f_LM_MB = vmax_lb * scalar_wb * MIC * LMWC / (parameters["kaff_lb"] + LMWC)
    }else{
      f_LM_MB=0
    }
    
    #Equation 16
    # MIC -> MAOM/LMWC
    if(MIC>0){
      f_MB_turn = parameters["rate_bd"] * MIC^2.0
    }else{
      f_MB_turn=0
    }
    
    #Equation 18
    # MAOM -> AGG
    if(MAOM>0){
      f_MA_AG = parameters["rate_ma"] * scalar_wd * MAOM
    }else{
      f_MA_AG=0
    }
    
    #Equation 22
    # microbial growth flux, but is not used in mass balance
    
    #Equation 21
    # MIC -> atmosphere
    if(MIC>0 && LMWC>0){
      f_MB_atm = f_LM_MB * (1 - (parameters["cue_ref"]- parameters["cue_t"] * (forc_st_prova[step.num] - parameters["tae_ref"]) ) )
    }else{
      f_MB_atm=0
    }
    
    
    #print("fluxes")
    #print("f_MB_atm,f_MA_AG,f_MB_turn,f_LM_MB,f_MA_LM,f_LM_MA,f_LM_leach,f_AG_break,f_PO_AG,f_PO_LM")
    #print(c(f_MB_atm,f_MA_AG,f_MB_turn,f_LM_MB,f_MA_LM,f_LM_MA,f_LM_leach,f_AG_break,f_PO_AG,f_PO_LM))
    
    # Update state variables
    
    #Equation 1
    dPOM = forc_npp * parameters["param_pi"] + f_AG_break * parameters["param_pa"] - f_PO_AG - f_PO_LM
    
    #Equation 8
    dLMWC = forc_npp * (1. - parameters["param_pi"]) - f_LM_leach + f_PO_LM - f_LM_MA - f_LM_MB + f_MB_turn * (1. - parameters["param_pb"]) + f_MA_LM
    
    #Equation 17
    dAGG = f_MA_AG + f_PO_AG - f_AG_break
    
    #Equation 19
    dMAOM = f_LM_MA - f_MA_LM + f_MB_turn * parameters["param_pb"] - f_MA_AG + f_AG_break * (1. - parameters["param_pa"])
    
    #Equation 20
    dMIC = f_LM_MB - f_MB_turn - f_MB_atm
    
    diff_eq <- c(dPOM, dLMWC, dAGG, dMIC, dMAOM )
    names(diff_eq) <- c("dPOM", "dLMWC", "dAGG", "dMIC", "dMAOM")
    
    return(list(diff_eq))
  })
}

Model.fwd <- function (parameters, times=run.steps) {
  output <- ode(y = state.SS, times=run.steps, func=decomp, parms =parameters, temp_func=forc_st_fw,water_func=forc_sw_fw, method="rk4") #solve ode, return output
  #print(length(forc_st_fw[[1]])/365)
  return(as.data.frame(cbind(time = output[run.steps.minus.one,"time"], POM = output[run.steps.minus.one,"POM"], LMWC = output[run.steps.minus.one,"LMWC"], AGG = output[run.steps.minus.one,"AGG"], MIC = output[run.steps.minus.one,"MIC"], MAOM = output[run.steps.minus.one,"MAOM"])))
}

##################################
#--------Data intialization-----
##################################

list_SOC_0 <-rep(NA,length(site_names_all))

list_forc_npp <-rep(NA,length(site_names_all))
list_clay <-rep(NA,length(site_names_all))
list_silt <-rep(NA,length(site_names_all))
list_pH <-rep(NA,length(site_names_all))
list_exper_len <-rep(NA,length(site_names_all))
list_BD <-rep(NA,length(site_names_all))

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

j<-1
for(site in site_names_all){
  #print(site)
  #site='CHNO3'
  site_df <- subset(C_input_exp, C_input_exp$ID.Site == site)
  site_T0_name<-site_T0_array_all[j]
  site_T0 <- subset(site_df,site_df$ID.Treatment==site_T0_name)
  #print(site_T0)
  #Choose years of simulations
  year_site <-subset(site_T0$Year, site_T0$ID.Site == site)
  exper_len <- length(year_site)
  year_0 <- min(year_site)
  year_end <- max(year_site)
  date_init_ss <- as.character(year_0-30)
  date_end_ss <- as.character(year_0-1)
  
  #print(date_init_ss)
  #print(date_end_ss)
  
  #Run for one year less than experiment length and add spinup up after
  #Se non vuoi -> run.steps <- seq(1,exper_len*365)
  run.steps <- seq(1,(exper_len-1)*365)
  run.steps.minus.one <- run.steps[1:length(run.steps)-1]
  
  #Read data spinup
  soil_temp_ss <- paste0(inputdir,'SCRIPT_MODELLI/',site,'/',"temp_",site,"_",date_init_ss,"_",date_end_ss,".txt")
  soil_hum_ss  <- paste0(inputdir,'SCRIPT_MODELLI/',site,'/',"hum_",site,"_",date_init_ss,"_",date_end_ss,".txt")
  forc_st_table <- read.table(soil_temp_ss) #degreeC
  forc_sw_table <- read.table(soil_hum_ss)/100. #kg(H2O)/m2(soil) to m3(H2O)/m3(soil)
  
  if(site=='RHEU1' | site == 'BOLO' | site == 'GRABOW'){
    forc_st_table=11.67
  }
  
  
  above <-mean(as.numeric(subset(site_T0$ABOVE, site_T0$ID.Site == site)), na.rm=TRUE)
  below<-mean(as.numeric(subset(site_T0$BELOW, site_T0$ID.Site == site)), na.rm=TRUE)
  forc_npp <- (above+below)*100/365 #gC/m2/day
  pH <-mean(subset(site_T0$pH, site_T0$ID.Site == site))
  clay <-mean(as.numeric(subset(site_T0$Clay, site_T0$ID.Site == site)))*100 #%
  silt <-mean(as.numeric(subset(site_T0$Silt, site_T0$ID.Site == site)))*100 #%
  bdens<-mean(as.numeric(unlist(subset(site_T0["Bulk density"], site_T0$ID.Site == site))))*1000 #kgsoil/m3
  
  
  #forc_sw <- approxfun(1:SStime, rep(mean(as.numeric(forc_sw_table$V1)),SStime)) #moisture input function
  if(site=='RHEU1'| site == 'BOLO' | site == 'GRABOW'){
    forc_st_ss <- rep(mean(as.numeric(forc_st_table)),SStime)}
  else{
    forc_st_ss<-rep(mean(as.numeric(forc_st_table$V1)),SStime)
  }
  
  forc_sw_ss <-rep(mean(as.numeric(forc_sw_table$V1)),SStime)
  
  
  #Forward
  #Read data forward
  date_init <- year_0
  date_end <- year_end
  soil_temp <- paste0(inputdir,'SCRIPT_MODELLI/',site,'/',"temp_",site,"_",date_init,"_",date_end,".txt")
  soil_hum  <- paste0(inputdir,'SCRIPT_MODELLI/',site,'/',"hum_",site,"_",date_init,"_",date_end,".txt")
  forc_st_table_fw <- read.table(soil_temp) #degreeC
  forc_sw_table_fw <- read.table(soil_hum)/100. #kg(H2O)/m2(soil) to m3(H2O)/m3(soil)
  
  
  forc_st_fw<-forc_st_table_fw$V1[1:length(run.steps)]
  forc_sw_fw<-forc_sw_table_fw$V1[1:length(run.steps)]
  
  
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
  list_silt[j] <- silt
  list_pH[j] <- pH
  list_exper_len[j] <-exper_len
  list_BD[j] <- bdens
  
  list_date_init_ss[j] <-date_init_ss
  list_date_end_ss[j] <-date_end_ss
  list_date_init[j] <-date_init
  list_date_end[j] <-date_end
  
  list_temp_ss <- append(list_temp_ss,list(forc_st_ss))
  list_water_ss <- append(list_water_ss,list(forc_sw_ss))
  list_temp_fw <- append(list_temp_fw,list(forc_st_fw))
  list_water_fw <- append(list_water_fw,list(forc_sw_fw))
  
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

rel_variance_df<-mean(relative_variances[relative_variances>0])

##################################
#--------Spinup + Forward run-----
##################################
list_SOC_sdvers<-list()
list_SS <- rep(NA,length(site_names_all))
list_spinup_nonOPT <-list()
j<-1
for(site in site_names_all){
  
  #Initialize
  forc_npp <-list_forc_npp[j]
  clay <-list_clay[j]
  silt <-list_silt[j]
  if(is.na(silt)==TRUE){
    silt<-(100.-clay)/2.
  }
  pH <-list_pH[j]
  exper_len<-list_exper_len[j]
  bdens<-list_BD[j]
  
  param_list["param_claysilt"]<-clay+silt  #% of soil in the clay and silt fraction
  param_list["param_bulkd"]<-bdens  
  
  
  #print(paste("exper_len",exper_len))
  date_init_ss<-list_date_init_ss[j]
  date_end_ss<-list_date_end_ss[j]
  date_init<-list_date_init[j]
  date_end<-list_date_end[j]
  
  forc_st_ss<-list_temp_ss[j]
  forc_sw_ss<-list_water_ss[j]
  forc_st_fw<-list_temp_fw[j]
  forc_sw_fw<-list_water_fw[j]
  #print(paste("len forc st fw", length(forc_st_fw[[1]])/365))
  
  SOC_0<-list_SOC_0[j]
  SOC_site_na<-SOC_list[[j]]
  SOC_var_na<-SOC_var_list[[j]]
  
  run.steps <- seq(1,(exper_len-1)*365)
  run.steps.minus.one <- run.steps[1:length(run.steps)-1]
  
  #STODE
  #SS_prova <- stode(y = state, time = SStime, 
  #                   func = decomp, parms = param_list, positive=TRUE,verbose = TRUE,rtol = 1e-6)
  #Spinup
  SS_prova <- stode(y = state, time = SStime, 
                    func = decomp, parms =param_list, temp_func=forc_st_ss,water_func=forc_sw_ss, positive=TRUE, rtol = 1e-4)
  
  soultions_SS<-SS_prova$y #gC/m2
  soultions_SS_tC<-soultions_SS/100 #tC/ha
  tot_SS <- sum(soultions_SS_tC)
  print(paste(site," total C spinup (tC/ha):",tot_SS))
  list_SS[j]=tot_SS
  
  
  #Model forcing
  state.SS = SS_prova$y
 
  #Save spinup
  list_spinup_nonOPT<-append(list_spinup_nonOPT,list(soultions_SS))
   
  #forward simu
  model_out <- Model.fwd(param_list,run.steps)
  
  print(model_out)
  #Annual simulated SOC
  n<-365
  yearly_model_out<-aggregate(model_out,list(rep(1:(nrow(model_out)%/%n+1),each=n,len=nrow(model_out))),mean)[-1]
  
  #print(yearly_model_out)
  yearly_SOC <- rowSums(yearly_model_out[c(2:6)])/100 #tC/ha
  
  #Add spinup
  yearly_SOC<-c(tot_SS,yearly_SOC)
  yearly_SOC_df <- cbind(seq(1:exper_len),yearly_SOC)
  
  #Annual observed SOC
  yearly_SOC_simu_vs_obs <- cbind(yearly_SOC_df,SOC_site_na)
  print("Simulated vs observed SOC (tC/ha)")
  print(yearly_SOC_simu_vs_obs)
  
  list_SOC_sdvers<-append(list_SOC_sdvers,list(yearly_SOC))
  
  j=j+1
}

#print(list_SS)
#print(list_SOC_0)
SOCSim<-cbind(site_names_all,list_SOC_0,list_SS)
print(SOCSim)

#print(list_SOC_sdvers)

##################
#OPTIMIZATION
##################
list_opt_param<-list()
list_SOC_simu <- list()
list_spinup_OPT <-list()
j<-1
for (site in site_names_all){
  print("******************")
  print("******************")
  print(paste("SITE:",site))
  print("******************")
  print("******************")
  #Initialize
  forc_npp <-list_forc_npp[j]
  clay <-list_clay[j]
  silt <- list_silt[j]
  if(is.na(silt)==TRUE){
    silt<-(100.-clay)/2.
  }
  pH <-list_pH[j]
  bdens <- list_BD[j]
  exper_len<-list_exper_len[j]
  date_init_ss<-list_date_init_ss[j]
  date_end_ss<-list_date_end_ss[j]
  date_init<-list_date_init[j]
  date_end<-list_date_end[j]
  SOC_0<-list_SOC_0[j]
  SOC_site_na<-SOC_list[[j]]
  SOC_var_na<-SOC_var_list[[j]]
  
  SOC_simu_sd <-list_SOC_sdvers[[j]]
  SOC_simu_sd_obs <- cbind(SOC_site_na,SOC_simu_sd)
  SOC_simu_sd_obs <- na.omit(SOC_simu_sd_obs)
  
  run.steps <- seq(1,(exper_len-1)*365)
  run.steps.minus.one <- run.steps[1:length(run.steps)-1]
  
  forc_st_ss<-list_temp_ss[j]
  forc_sw_ss<-list_water_ss[j]
  forc_st_fw<-list_temp_fw[j]
  forc_sw_fw<-list_water_fw[j]
  
  
  #Add site specific parameters to param_list
  param_list["param_claysilt"]<-clay+silt  #% of soil in the clay and silt fraction
  param_list["param_bulkd"]<-bdens  
  
  #OPTIM ALGO
  
  prior_param1<- param_list['eact_lb']
  prior_param2 <- param_list['eact_pl']
  prior_param3<-param_list['kaff_pl']
  prior_params<-c(prior_param1,prior_param2,prior_param3)
  
  optim_func<-optim(prior_params,Jparam_new, method = "L-BFGS-B",
                    lower = c(prior_param1-1e4,prior_param2-1e4,prior_param3-1e3), upper = c(prior_param1+1e4,prior_param2+1e4,prior_param3+1e3))
  
  #optim_func<-optim(prior_params,Jparam_new, method = "BFGS")
  print(optim_func)
  
  print("Optim param")
  print(optim_func$par)
  print("Optim Jparam")
  print(optim_func$value)
  print("Convergence")
  print(optim_func$convergence)
  
  list_opt_param<-append(list_opt_param,list(optim_func$par))
  
  #Optimized version 
  param_list_opt <- param_list
  new_params<-optim_func$par
  param_list_opt['eact_lb'] <- new_params[1]
  param_list_opt['eact_pl'] <- new_params[2]
  param_list_opt['kaff_pl'] <- new_params[3]
  
  print("Optimized param_list")
  print(param_list_opt)  
  opt_SS <- stode(y = state, time = SStime, 
                  func = decomp, parms = param_list_opt, temp_func=forc_st_ss,water_func=forc_sw_ss, positive=TRUE,rtol = 1e-4)
  opt_ssp = sum(opt_SS$y)
  opt_ssp_tC = opt_ssp/100.
  

  #Model forcings
  state.model_opt = opt_SS$y
  #Save spinup
  list_spinup_OPT<-append(list_spinup_OPT,list(state.model_opt))
  #Forward
  model_out_opt <- Model.fwd(param_list_opt,run.steps)
  
  #Annual simulated SOC
  n<-365
  yearly_model_out_opt<-aggregate(model_out_opt,list(rep(1:(nrow(model_out_opt)%/%n+1),each=n,len=nrow(model_out_opt))),mean)[-1]
  yearly_SOC_opt <- rowSums(yearly_model_out_opt[c(2:6)])/100 #tC/ha
  
  #Add spinup
  yearly_SOC_opt<-c(opt_ssp_tC,yearly_SOC_opt)
  yearly_SOC_df_opt <- cbind(seq(1:exper_len),yearly_SOC_opt)
  
  #Annual observed SOC (and var)
  opt_yearly_SOC_simu_vs_obs <- cbind(yearly_SOC_df_opt,SOC_site_na)
  #selct model values where data avail
  opt_yearly_SOC_simu_vs_obs <- na.omit(opt_yearly_SOC_simu_vs_obs)
  
  #print("Standard version")
  #print(SOC_simu_sd_obs)
  
  sd_opt_SOC <- cbind("SOC_obs"=SOC_simu_sd_obs[,'SOC_site_na'],"SOC_simu_sd"=SOC_simu_sd_obs[,'SOC_simu_sd'],"SOC_simu_opt"=opt_yearly_SOC_simu_vs_obs[,'yearly_SOC_opt'])
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

#--------STOP-----------
#Print SOC fits
#Plot
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

