## FEISTY functions
##
## Authors: Daniel Ottmann 
## Email: daniel.ottmann.riera@gmail.com
##
## Date created: January 2024
## Last update:  March 2025
##
## ---------------------------
##
## Readme:
##
## This modifies functions from the FEISTY model to extract values carbon fluxes through functional groups
## and to run a sensitivity analysis.
##
## ---------------------------

# load libraries:
library(tidyverse)
library(FEISTY)


# ---
paramAddPhysiology_sensitivity <- function (p, 
                                            ac = 20,          # Max. consumption coefficient  [g^(-n)/yr]
                                            bc = -0.25 * 1,       # Max. consumption exponent     [-]
                                            am = 0.2*20 * 1,      # Metabolism coefficient        [g^(-p)/yr] 0.2*ac
                                            bm = -0.175 * 1,      # Metabolism exponent           [-]
                                            ae = 70,          # Coef. for clearance rate      [m2*g^(q-1)/yr] encounter slope
                                            be = -0.2,        # Clearance rate exponent
                                            epsRepro = 0.01,  # reproduction * recruitment efficiency )
                                            epsAssim = 0.7    # Assimilation efficiency  
)
{
  # index pointing to fish stages and the fish weights:
  
  ix = p$ixFish
  m  = p$mc[ix]  # size of fish
  
  # maximum consumption rate, [/yr] 
  p$Cmax[ix]       = ac*m^bc          
  
  # standard metabolism (basal respiration), [/yr]
  p$metabolism[ix] = am*m^bm                   
  
  # clearance rate, search rate [m2/g/yr]  
  p$V[ix]          = ae*m^be                       
  
  # reproduction * recruitment efficiency , [-]
  p$epsRepro = rep(epsRepro, times=p$nGroups)  
  
  # Assimilation efficiency, [-]
  p$epsAssim = epsAssim                   
  
  # remove the NAs
  p$mc       [is.na(p$mc)]        <- 0 
  p$u0       [is.na(p$u0)]        <- 0 
  p$z        [is.na(p$z)]         <- 0 
  p$psiMature[is.na(p$psiMature)] <- 0 
  p$mortF    [is.na(p$mortF)]     <- 0 
  p$mort0    [is.na(p$mort0)]     <- 0 
  p$Cmax     [is.na(p$Cmax)]      <- 0 
  p$metabolism[is.na(p$metabolism)] <- 0 
  p$V        [is.na(p$V)]         <- 0 
  
  #Oct 2023 add for Temperature effects
  p$Cmaxsave=p$Cmax
  p$Vsave=p$V
  p$metabolismsave=p$metabolism
  
  names(p$u0) =  p$stagenames  
  iF = p$ixFish
  p$fishes = data.frame(mc=p$mc[iF], mLower=p$mLower[iF], mUpper=p$mUpper[iF], 
                        z=p$z[iF], psiMature=p$psiMature[iF], mortF=p$mortF[iF], 
                        mort0=p$mort0[iF], Cmax=p$Cmax[iF], 
                        metabolism=p$metabolism[iF], V=p$V[iF] )
  row.names(p$fishes)= p$stagenames[iF]  
  p$groups = data.frame(epsRepro=p$epsRepro, epsAssim=p$epsAssim)  
  row.names(p$groups) = p$groupnames[-p$ixR]
  return(p)
}


# ----
my_setupvertical2 <- function(szprod = 80, # small zoo production
                              lzprod = 80, # large zoo production
                              bprodin = NA, # benthos production
                              dfbot  = NA, # detrital flux reaching the bottom
                              dfpho  = NA, # detrital flux out of photic zone
                              nStages = 9, # No. of size groups
                              Tp = NA, # Average T of top 100 m (up to 100 m). Default 10 Celsius.
                              Tm = NA, # Average T of 500 - 1500 m (up to 1500 m). Default 10 Celsius. Keep it as NA, if no Tm data. Tm = Tb.
                              Tb = NA, # Bottom T (last layer value). Default 10 Celsius.
                              depth = 800, # Bottom depth
                              photic = 150, # Photic zone depth
                              shelfdepth = 250, # shelf region depth
                              visual = 1.5, # >1 visual predation primarily during the day, = 1 equal day and night
                              etaMature = 0.25, # Size of matureation relative to
                              # asymptotic size. Different from
                              # van Denderen (2021), where it is 0.002
                              Fmax = 0,
                              etaF = 0.05,
                              sensitivity_demmig = 0.9 ) { # Sensitivity of demersal fish daily upward migration (0.9 = 10% shallower;1.1 = 10% deeper)
  # benthic production calc
  if (is.na(bprodin) & is.na(dfbot) & is.na(dfpho)){ # if all benthic arguments are NA, assign bprod to 5
    bprodin = -1; dfbot = -1; dfpho = 350
    bprod=0.1*(dfpho*(depth/photic)^-0.86)
    if(bprod>=0.1*dfpho) bprod=0.1*dfpho
  } else {
    if (sum(!is.na(c(bprodin, dfbot, dfpho)))>1) stop('Please check "bprod" and "dfbot" input. Only one of them should be assigned values, others should be kept as "NA".')
    if (!is.na(bprodin)) {bprod = bprodin} else {bprodin = -1}
    if (!is.na(dfbot)) {bprod = dfbot*0.1} else {dfbot = -1}
    if (!is.na(dfpho)) {bprod=0.1*(dfpho*(depth/photic)^-0.86); if(bprod>=0.1*dfpho) bprod=0.1*dfpho} else {dfpho = -1}
  }
  
  # Temperature initialize
  if (is.na(Tp)) Tp = 10
  if (is.na(Tb)) Tb = 10
  if (is.na(Tm)) Tm = Tb # if Tm is not provided, Tm = Tb
  
  #------------------  
  # Initialize the parameters:
  # habitat and small benthos
  #------------------  
  
  param = paramInit(bottom=depth, szprod=szprod, lzprod=lzprod, photic=photic,
                    shelfdepth=shelfdepth, visual=visual, bprodin=bprodin, dfbot=dfbot, dfpho=dfpho, bprod=bprod,
                    etaMature=etaMature, Tp=Tp, Tm=Tm, Tb=Tb)
  
  #------------------  
  # Setup resource groups:
  #------------------  
  
  param = paramAddResource(
    param, 
    names= c("smallZoo", "largeZoo", "benthos", "Spare_position"),
    K    = c(szprod, lzprod, bprod, 0),  # g ww/m2  - maximum resource concentration
    r    = c(1, 1, 1, 1),              # [/yr] nudging coefficient
    mc   = c(2e-06*sqrt(500), 0.001*sqrt(500), 1e-04*sqrt(250000), 0.25*sqrt(500)),
    mLower = c(2e-06,0.001, 0.5e-03, 0.25), # weight lower limit
    mUpper = c(0.001, 0.5, 125, 125),
    u0     = c(0.5,0.5,0.5,0))
  
  #------------------  
  # Add fish groups:
  #------------------  
  nSmall = round(0.66*nStages)
  # mMature=NA overrides the generic psiMature-> only adult classes 50% mature
  u0  = 0.0001
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, u0=u0,
                        mortF=0,      nStages=nSmall, name="smallPel")
  
  
  u0M = u0  # initial condition = 0 if no mesopelagic zone
  if (param$bottom <= param$shelfdepth) u0M <- 0
  
  param = paramAddGroup(param, mMin=0.001, mMax=   250, mMature=etaMature*250, u0=u0M,
                        mortF=0,   nStages=nSmall, name="mesoPel")
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0, 
                        mortF=0, nStages=nStages, name="largePel") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0M, 
                        mortF=0, nStages=nStages, name="midwPred") 
  
  param = paramAddGroup(param, mMin=0.001, mMax=125000, mMature=etaMature*125000, u0=u0,
                        mortF=0, nStages=nStages, name="demersals")
  #param$mortF[length(param$mortF)]=0.5
  
  #------------------  
  # Setup physiology:
  #------------------  
  param = paramAddPhysiology_sensitivity(param)
  
  # Add fishing mortality
  param=setFishing(param, Fmax=Fmax, etaF=etaF)
  
  #------------------  
  # theta (preferences):
  #------------------  
  
  param$vertover   = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  # calculate size-preference matrix
  param$sizeprefer=paramSizepref(p=param,           # parameter settings 
                                 beta = 400,  # preferred predator/prey mass ratio
                                 sigma = 1.3, # width of size preference for feeding
                                 type = 1)
  
  #------------------  
  # overlap from depth distribution
  #------------------  
  ssigma = 10    # width of initial distribution
  tau    = 10    # increase in width
  
  sigmap = ssigma + tau*log10(param$mc/param$mc[1]) # width for each size class
  xrange = 0 : param$bottom
  param$dvm = param$photic + 500 # 650
  
  if (param$bottom < (param$photic + 500)) 
    param$dvm = param$bottom   # migration to bottom in intermediate habitats
  
  if (param$bottom <= param$shelfdepth) 
    param$dvm = 0              # no migration in shallow habitats
  
  ixmedium = which.min(abs(param$mLower[param$ix[[5]]] - 0.5))# which.min(abs(param$mLower[param$ix[[5]]] - etaMature*250))
  ixlarge = which.min(abs(param$mLower[param$ix[[5]]] - 250))# which.min(abs(param$mLower[param$ix[[5]]] - etaMature*125000))
  
  # a function to generate vertical distributions (a normal distribution)
  VertDist <- function(sigma, xloc){
    xloc = rep(xloc, length.out=length(sigma))
    zp_n = matrix(nrow=length(xrange), ncol=length(sigma), data=0) 
    for (i in 1: length(sigma)){      
      zp_n[,i] = (1/(sqrt(2*pi*sigma[i]^2)))* 
        exp(-(((xrange - xloc[i])^2)/(2*sigma[i]^2)))
    }
    zp_n = zp_n %*% diag(1/colSums(zp_n))
    zp_n  
  }
  
  ## zooplankton : small zoo & large zoo
  # at night: zooplankton is close to surface
  zp_n = VertDist(sigmap[1:2], xloc=0)
  
  # zooplankton day (half at surface, half at dvm depth
  zp_d = VertDist(sigmap[1:2], xloc=param$dvm)
  zp_d = (zp_n + zp_d)/2
  
  ## benthos small and large (at bottom with width ssigma)
  bent_dn = VertDist(c(ssigma, ssigma), xloc=param$bottom)
  
  ## small pelagic fish (day + night) always at surface
  ix = param$ix[[1]]
  spel_dn = VertDist(sigmap[ix], xloc=0)
  
  ## meso pelagic night   at surface  
  mpel_n = spel_dn
  
  # meso pelagic day (all at dvm)
  ix = param$ix[[2]]
  mpel_d = VertDist(sigmap[ix], xloc=param$dvm)
  
  ## large pelagic fish night (all at surface)
  ix = param$ix[[3]]
  lpel_n = VertDist(sigmap[ix], xloc=0)
  
  # large pelagic fish day (non-large at surface   large half at surface half at dvm)
  xlocvec = rep(0,length(ix)) 
  xlocvec[ixlarge:length(xlocvec)] = param$dvm 
  lpel_d = VertDist(sigmap[ix], xloc=xlocvec)
  lpel_d = (lpel_d + lpel_n)/2
  
  ## bathypelagic night (large in midwater, others at surface)
  ix = param$ix[[4]]
  xlocvec = rep(0,length(ix)) # initialization
  xlocvec[ixlarge:length(xlocvec)] = param$dvm # non-large at surface   large at dvm
  bpel_n = VertDist(sigmap[ix], xloc=xlocvec)
  
  # bathypelagic day (all at dvm)
  bpel_d = VertDist(sigmap[ix], xloc=param$dvm)
  
  ## demersal fish night
  ix = param$ix[[5]]
  xlocvec = rep(0,length(ix)) # initialization
  xlocvec[ixmedium:length(xlocvec)] = param$bottom # small at surface   medium and large at bottom
  dem_n = VertDist(sigmap[ix], xlocvec)
  
  # demersal fish day; small at surface/ medium at bottom/ large at middle
  demmig = param$dvm * sensitivity_demmig# ? from matlab
  if ((param$bottom - param$dvm) >= 1200) 
    demmig = (param$dvm + param$bottom-param$dvm-1200) * sensitivity_demmig
  if ((param$bottom - param$dvm) >= 1500)
    demmig = param$bottom 
  
  dem_d= matrix(nrow=length(xrange), ncol=length(param$ix[[5]]), data=0)
  
  xlocvec[ixlarge:length(xlocvec)] = demmig #param$dvm ### or demmig???
  dem_d =  VertDist(sigmap[ix], xlocvec)
  
  #if shallower than euphotic depth, large demersals feed across-habitats
  if (param$bottom <= param$photic) {
    dem_d = (dem_d + dem_n)/2
    dem_n = dem_d
  }
  
  # calculate overlap during day
  param$depthDay = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  test     = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  param$dayout = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  param$depthDay[, 1:2] = zp_d # resources
  param$depthDay[, 3:4] = bent_dn # resources
  param$depthDay[, param$ix[[1]]] = spel_dn
  param$depthDay[, param$ix[[2]]] = mpel_d
  param$depthDay[, param$ix[[3]]] = lpel_d
  param$depthDay[, param$ix[[4]]] = bpel_d
  param$depthDay[, param$ix[[5]]] = dem_d
  
  for (i in 1: param$nStages) {
    for ( j in 1: param$nStages) {
      test[, j] = pmin(param$depthDay[, i], param$depthDay[, j])
    }
    param$dayout[, i] = colSums(test)
  }
  
  # calculate overlap during night
  param$depthNight = matrix(nrow=length(xrange), ncol=param$nStages, data=0)
  # test will be overwritten
  param$nightout = matrix(nrow=param$nStages, ncol=param$nStages, data=0)
  
  param$depthNight[, 1:2] = zp_n # resources
  param$depthNight[, 3:4] = bent_dn # resources
  param$depthNight[, param$ix[[1]]] = spel_dn
  param$depthNight[, param$ix[[2]]] = mpel_n
  param$depthNight[, param$ix[[3]]] = lpel_n
  param$depthNight[, param$ix[[4]]] = bpel_n
  param$depthNight[, param$ix[[5]]] = dem_n
  
  for (i in 1: param$nStages) {
    for ( j in 1: param$nStages) {
      test[, j] = pmin(param$depthNight[, i], param$depthNight[, j])
    }
    param$nightout[, i] = colSums(test)
  }
  
  #------------------  
  # visual ability
  #------------------  
  
  # visual predatars: good at light, bad in the dark
  visualpred = c(param$ix[[1]], # small palegic 5 6 always at surface
                 param$ix[[3]]) # large pelagic 9 10 11
  param$dayout[visualpred,]   = param$dayout[visualpred,]*param$visual       # predation enhanced during day
  param$nightout[visualpred,] = param$nightout[visualpred,]*(2-param$visual) # predation decreased at night 
  
  # pelagic predators: limited vision in twilight zone during day
  pelpred = param$ix[[3]]                    # large pelagic   9 10 11
  pelpred = pelpred[ixlarge:length(pelpred)] # large large pelagic  11  at dvm during day
  preytwi = c(param$ix[[2]], param$ix[[4]])  # mesopelagic 7 8   bathypelagic 12 13 14
  param$dayout[pelpred, preytwi] = param$dayout[pelpred, preytwi]/param$visual*(2 - param$visual)  # /1.5 to restore  then *0.5 
  
  # average overlap during the whole day
  param$vertover = (param$dayout + param$nightout)*0.5
  
  # calculate combined feeding preference matrix
  param$theta = param$sizeprefer*param$vertover
  
  #  specific revision of feeding preference
  idx_be = param$ixFish[1]: (param$ix[[5]][1] + (ixmedium - 2)) # all pelagic and small demersals
  param$theta[idx_be, 3:4] = 0 # all pelagic and small demersals do not eat benthos,
  # only medium & large demersals eat benthos
  
  # medium demersals are less preyed on
  idx_smd = (param$ix[[5]][1] + (ixmedium - 1)): (param$ix[[5]][1] + (ixlarge - 2)) #
  param$theta[idx_be, idx_smd] = param$theta[idx_be, idx_smd]*0.25
  
  # small & large demersals do not eat zooplankton
  param$theta[(param$ix[[5]][1] + (ixmedium - 1)) : param$ix[[5]][length(param$ix[[5]])], 1:2] = 0
  
  # provide benefit to forage and mesopelagic fish (predator avoidance)
  pred1 = (param$ix[[3]][1]+ (ixlarge - 1)) : param$ix[[3]][length(param$ix[[3]])]
  pred2 = (param$ix[[4]][1]+ (ixlarge - 1)) : param$ix[[4]][length(param$ix[[4]])]
  pred3 = (param$ix[[5]][1]+ (ixlarge - 1)) : param$ix[[5]][length(param$ix[[5]])]
  prey1 = (param$ix[[1]][1]+ (ixmedium   - 1)) : param$ix[[1]][length(param$ix[[1]])]
  prey2 = (param$ix[[2]][1]+ (ixmedium   - 1)) : param$ix[[2]][length(param$ix[[2]])]
  idx_predat = c(pred1, pred2, pred3)
  idx_prey   = c(prey1, prey2)
  param$theta[idx_predat,idx_prey] = param$theta[idx_predat,idx_prey]*0.5
  
  # update temperature
  Q10=1.88
  Q10m=1.88
  
  # initialize effective T vector (all resoources and fish)
  Teff = rep(0, length(param$u0))
  
  # zooplankton (no use)
  Tday = (param$Tp + param$Tm)/2 # half surface half dvm  param$dvm = param$photic + 500
  if (param$dvm == param$bottom) { # when param$bottom < (param$photic + 500)
    Tday = (param$Tp + param$Tb)/2
  }
  if (param$dvm == 0) Tday = param$Tp # when param$bottom <= param$shelfdepth
  Tnight = param$Tp # all surface
  Teff[1:2] = (Tday+Tnight)/2
  # benthos (no use)
  Teff[3:4] = param$Tb
  # small pelagics
  ix = param$ix[[1]]
  Teff[ix] = param$Tp # always surface
  # mesopelagics 
  ix = param$ix[[2]]
  Tday = param$Tm # dvm
  if (param$dvm == param$bottom) Tday = param$Tb
  if (param$dvm == 0) Tday = param$Tp
  Tnight = param$Tp # surface
  Teff[ix] = (Tday+Tnight)/2
  # large pelagics
  ix = param$ix[[3]]
  # daytime large half at surface half at dvm
  Tdaylarge = (param$Tp+param$Tm)/2
  if (param$dvm == param$bottom) Tdaylarge = (param$Tp+param$Tb)/2
  if (param$dvm == 0) Tdaylarge = param$Tp
  Tdaynonlarge = param$Tp # non-large at surface at daytime
  Tnight = param$Tp     # all at surface at night
  Teff[ix[ixlarge:length(ix)]]  = (Tdaylarge+Tnight)/2      # large
  Teff[ix[-(ixlarge:length(ix))]] = (Tdaynonlarge+Tnight)/2 # non-large
  # bathypelagics    
  ix = param$ix[[4]]
  Tday = param$Tm # all at dvm at daytime
  if (param$dvm == param$bottom) Tday = param$Tb
  if (param$dvm == 0)            Tday = param$Tp
  Tnightlarge = param$Tm # large at dvm
  if (param$dvm == param$bottom) Tnightlarge = param$Tb
  if (param$dvm == 0)            Tnightlarge = param$Tp
  Tnightnonlarge = param$Tp # non-large at surface at night
  Teff[ix[ixlarge:length(ix)]]  = (Tday+Tnightlarge)/2      # large
  Teff[ix[-(ixlarge:length(ix))]] = (Tday+Tnightnonlarge)/2 # non-large
  # demersals    
  ix = param$ix[[5]]
  Tsmall  = param$Tp # small at surface
  Tmedium = param$Tb # medium at bottom
  # large
  # daytime
  Tdaylarge = param$Tm # large at middle
  # if the water is very deep large demersals always stay at the bottom
  if ((param$bottom - param$dvm) >= 1500) Tdaylarge = param$Tb
  # if the water is very shallow large demersals migrate over the whole water column both day and night
  if (param$bottom <= param$photic) {
    Tdaylarge = (param$Tp + param$Tb)/2
  }
  # nighttime
  Tnightlarge  = param$Tb # large at bottom if water is deep enough
  # if the water is very shallow large demersals migrate over the whole water column both day and night
  if (param$bottom <= param$photic) {
    Tnightlarge = (param$Tp + param$Tb)/2
  }
  
  Teff[ix[-(ixmedium:length(ix))]] = Tsmall # small
  Teff[ix[ixmedium:(ixlarge-1)]]   = Tmedium # medium
  Teff[ix[ixlarge:length(ix)]]  = (Tdaylarge+Tnightlarge)/2 # large
  
  scTemp =  Q10^((Teff-10)/10)
  scTempm =  Q10m^((Teff-10)/10)
  
  param$Cmax = scTemp* param$Cmax # maximum consumption rate 
  param$V= scTemp* param$V # clearance rate 
  param$metabolism = scTempm* param$metabolism
  
  param$setup="setupVertical2"
  
  return(param)  
}


#----------------------
run_feisty <- function(data) {
  # Create output files
  suppressMessages({
    out <- data.frame()
    out2 <- data.frame()
    id <- data$id
    
    #--------------------------------------------------
    # For loop for each station:
    
    for(i in 1:nrow(data)) {
      
      # Customize FEISTY setup:
      my_setup <- setupVertical2(szprod = data[i, 4], lzprod = data[i, 5], # Pelagic productivities
                                 dfbot = data[i, 6] * 1, # detrital flux reaching the bottom 
                                 nStages = 9, # No. of size groups
                                 Tp = as.numeric(data[i, 9]), # Temperature of the surface
                                 Tb = as.numeric(data[i, 10]), # Temperature of the bottom
                                 Tm = as.numeric(data[i, 11]), # Temperature of the midwater
                                 depth = as.numeric(data[i, 7]), # Bottom depth
                                 photic = as.numeric(data[i, 8])) # photic depth
      
      # # Sensitivity to maturation size of demersals:
      # etaMature_new <- 0.002 # your new eatMature factor
      # my_setup$mMature[5] <- etaMature_new * 125000 # new size with 50% maturity
      # ix <- my_setup$ix[[5]] # demersal fish index in the vector
      # my_setup$psiMature[ix] <- ( 1 + (my_setup$mc[ix]/my_setup$mMature[5])^(-5) )^(-1)
      
      
      
      # Run FEISTY:
      sim <- simulateFEISTY(bCust    = T,
                            p      = my_setup, 
                            tEnd   = 200,
                            USEdll = TRUE) 
      
      #--------------------------------------------------
      # Take average biomass of each functional group of the last 80 years:
      
      totBiomass <- sim$totBiomass %>%
        as.data.frame() %>%
        tail(80) %>%
        summarise(across(everything(), mean)) %>%
        rename(smallPel = totBiomass.smallPel,
               mesoPel = totBiomass.mesoPel,
               largePel = totBiomass.largePel,
               midwPred = totBiomass.midwPred,
               demersals = totBiomass.demersals) %>%
        t() %>%
        as.data.frame() %>%
        rename(totBiomass = V1) %>%
        rownames_to_column(var = "funGroup") %>%
        mutate(propBiomass = totBiomass / sum(totBiomass)) %>%
        mutate(id = id[i])
      
      totResources <- sim$R %>%
        as.data.frame() %>%
        tail(80) %>%
        summarise(across(everything(), mean)) %>%
        mutate(Zoo = smallZoo + largeZoo) %>%
        dplyr::select(Zoo, benthos) %>%
        t() %>%
        as.data.frame() %>%
        rename(totBiomass = V1) %>%
        rownames_to_column(var = "funGroup") %>%
        mutate(propBiomass = NA) %>%
        mutate(id = id[i])
      
      out <- rbind(out, totResources, totBiomass)
      
      #--------------------------------------------------
      # Stomach content of functional groups:
      
      # Extract relevant information from the simulation:
      p=sim$p
      u=sim$u
      p$nstage <-lengths <- max(sapply(p$ix, length)) #maximum number of stages for one group
      biomass <- u
      Bin <- round(0.6 * nrow(biomass), digits = 0) 
      biomassend <- colMeans(biomass[Bin:nrow(biomass),]) # mean value of the last 40% time 
      biomassstage <- p$ixFish[length(p$ixFish)]
      biomasssmall <- p$nstage - round(2/3*p$nstage, digits = 0)
      Enc = p$V * (p$theta %*% biomassend)
      f   = Enc / (p$Cmax + Enc)
      f[is.na(f)] = 0  
      bom <- t(t(p$theta[5:biomassstage, ]) * colMeans(biomass[Bin:nrow(biomass),])) 
      fbom <- f[5:biomassstage] / rowSums(bom)
      output <- bom * fbom

      # Calculate "ingestion" of juvenile and adult demersals (smaller stages are pelagic) - total amount of prey consumed by demersals:
      ingestion <- f * p$Cmax * biomassend
      ingestion <- ingestion %>%
        as.data.frame() %>%
        rename(tot_ingestion = V1) %>%
        rownames_to_column(var = "stage") %>%
        filter(str_detect(stage, "demersals")) %>%
        mutate(stage = as.numeric(str_replace_all(stage, "demersals_", "")),
               id = id[i]) %>%
        filter(stage > ((1/3) * p$nstage)) %>%
        relocate(id)
      
      # Define carbon paths to demersals depending on food source:
      pelagic_spp <- c("smallZoo", "largeZoo", "smallPel", "largePel")
      benthic_spp <- c("smallBenthos", "largeBenthos")
      demersal_spp <- "demersals"
      midwater_spp <- c("mesoPel", "midwPred")
      
      # Differentiate systems with/without mesopelagic and midwPred fish:
      if (length(p$ix)==5){
        
        p$SpId <- c('smallPel','mesoPel','largePel', 'midwPred', 'demersals')
        SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
                  rep(p$SpId[1], length(p$ix[[1]])),
                  rep(p$SpId[2], length(p$ix[[2]])),
                  rep(p$SpId[3], length(p$ix[[3]])),
                  rep(p$SpId[4], length(p$ix[[4]])),
                  rep(p$SpId[5], length(p$ix[[5]])))
        
        p$RSpName <- c("Small zooplankton", "Big zooplankton", "Benthos", "Small pelagics",
                       "Mesopelagics", "Large pelagics", "Midwater predators", "Demersals")
        
        # Filter demersals and add flux pathway: 
        demers <- t(output[(p$ix[[5]][1] - length(p$ixR)):(p$ix[[5]][length(p$ix[[5]])] - length(p$ixR)), ])
        demers <- data.frame(val = c(demers), 
                             stage = rep(1:p$nstage, each = nrow(demers)), 
                             SpId = as.factor(rep(SpId, p$nstage))) %>%
          mutate(flux_pathway = case_when(SpId %in% pelagic_spp ~ "Epipelagic",
                                          SpId %in% benthic_spp ~ "Benthic",
                                          SpId %in% demersal_spp ~ "Demersal",
                                          T ~ "Midwater")) %>%
          filter(stage > ((1/3) * p$nstage))
        
        # Group by SpId and stage, keep flux pathway:
        demers <- demers %>%
          group_by(stage, SpId, flux_pathway) %>%
          summarise(feeding = sum(val))
        
      } else {
        
        # Differentiate systems with/without mesopelagic and midwater predators fish:
        p$SpId <- c('smallPel','largePel', 'demersals')
        SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
                  rep(p$SpId[1], length(p$ix[[1]])),
                  rep(p$SpId[2], length(p$ix[[2]])),
                  rep(p$SpId[3], length(p$ix[[3]])))
        
        p$RSpName <- c("Small zooplankton", "Large zooplankton", "Benthos", "Small pelagics",
                       "Large pelagics", "Demersals")   
        
        # Filter demersals and add flux pathway: 
        demers <- t(output[(p$ix[[3]][1] - length(p$ixR)):(p$ix[[3]][length(p$ix[[3]])] - length(p$ixR)), ])
        demers <- data.frame(val = c(demers), 
                             stage = rep(1:p$nstage, each = nrow(demers)), 
                             SpId = as.factor(rep(SpId, p$nstage))) %>%
          mutate(flux_pathway = case_when(SpId %in% pelagic_spp ~ "Epipelagic",
                                          SpId %in% benthic_spp ~ "Benthic",
                                          SpId %in% demersal_spp ~ "Demersal",
                                          T ~ "Midwater")) %>%
          filter(stage > ((1/3) * p$nstage))
        
        # Group by SpId and stage, keep flux pathway:
        demers <- demers %>%
          group_by(stage, SpId, flux_pathway) %>%
          summarise(feeding = sum(val))
        
      }
      
      # Add ingestion and calculate prey-speciffic ingestion:
      demers <- demers %>%
        mutate(id = id[i]) %>%
        relocate(id) %>%
        group_by(stage) %>%
        mutate(sum_feeding = sum(feeding),
               rel_feeding = feeding/sum(feeding)) %>%
        left_join(ingestion, by = c("id", "stage")) %>%
        mutate(ingestion = tot_ingestion * rel_feeding)
      
      
      out2 <- rbind(out2, demers)
    }
  })  
  
  return(list(Biomass = out,
              diet = out2))
}




#----------------------
run_feisty_sensitivity <- function(data) {
  # Create output files
  suppressMessages({
    out <- data.frame()
    out2 <- data.frame()
    id <- data$id
    
    #--------------------------------------------------
    # For loop for each station:
    
    for(i in 1:nrow(data)) {
      
      # Customize FEISTY setup:
      my_setup <- my_setupvertical2(szprod = data[i, 4], lzprod = data[i, 5], # Pelagic productivities
                                 dfbot = data[i, 6]*1, # detrital flux reaching the bottom # <-- sensitivity *2 to compare with Clive's results and trawl surveys
                                 nStages = 9, # No. of size groups
                                 Tp = as.numeric(data[i, 9]), # Temperature of the surface
                                 Tb = as.numeric(data[i, 10]), # Temperature of the bottom
                                 Tm = as.numeric(data[i, 11]), # Temperature of the midwater
                                 depth = as.numeric(data[i, 7]), # Bottom depth
                                 photic = as.numeric(data[i, 8])) # photic depth
      
      # # Sensitivity to maturation size of demersals:
      my_setup$metabolism  <- my_setup$metabolism * 1
      my_setup$Cmax <- my_setup$Cmax * 1
      my_setup$epsAssim <- my_setup$epsAssim * 1
      # Sensitivity to vertical migration of demersal fishes is coded in my_setupvertical2 function
      
      # etaMature_new <- 0.002 # your new eatMature factor
      # my_setup$mMature[5] <- etaMature_new * 125000 # new size with 50% maturity
      # ix <- my_setup$ix[[5]] # demersal fish index in the vector
      # my_setup$psiMature[ix] <- ( 1 + (my_setup$mc[ix]/my_setup$mMature[5])^(-5) )^(-1)
      

      # Run FEISTY:
      sim <- simulateFEISTY(bCust    = T,
                            p      = my_setup, 
                            tEnd   = 200,
                            USEdll = TRUE) 
      
      #--------------------------------------------------
      # Take average biomass of each functional group of the last 80 years:
      
      totBiomass <- sim$totBiomass %>%
        as.data.frame() %>%
        tail(80) %>%
        summarise(across(everything(), mean)) %>%
        rename(smallPel = totBiomass.smallPel,
               mesoPel = totBiomass.mesoPel,
               largePel = totBiomass.largePel,
               midwPred = totBiomass.midwPred,
               demersals = totBiomass.demersals) %>%
        t() %>%
        as.data.frame() %>%
        rename(totBiomass = V1) %>%
        rownames_to_column(var = "funGroup") %>%
        mutate(propBiomass = totBiomass / sum(totBiomass)) %>%
        mutate(id = id[i])
      
      totResources <- sim$R %>%
        as.data.frame() %>%
        tail(80) %>%
        summarise(across(everything(), mean)) %>%
        mutate(Zoo = smallZoo + largeZoo) %>%
        dplyr::select(Zoo, benthos) %>%
        t() %>%
        as.data.frame() %>%
        rename(totBiomass = V1) %>%
        rownames_to_column(var = "funGroup") %>%
        mutate(propBiomass = NA) %>%
        mutate(id = id[i])
      
      out <- rbind(out, totResources, totBiomass)
      
      #--------------------------------------------------
      # Stomach content of functional groups:
      
      # Extract relevant information from the simulation:
      p=sim$p
      u=sim$u
      p$nstage <-lengths <- max(sapply(p$ix, length)) #maximum number of stages for one group
      biomass <- u
      Bin <- round(0.6 * nrow(biomass), digits = 0) 
      biomassend <- colMeans(biomass[Bin:nrow(biomass),]) # mean value of the last 40% time 
      biomassstage <- p$ixFish[length(p$ixFish)]
      biomasssmall <- p$nstage - round(2/3*p$nstage, digits = 0)
      Enc = p$V * (p$theta %*% biomassend)
      f   = Enc / (p$Cmax + Enc)
      f[is.na(f)] = 0  
      bom <- t(t(p$theta[5:biomassstage, ]) * colMeans(biomass[Bin:nrow(biomass),])) 
      fbom <- f[5:biomassstage] / rowSums(bom)
      output <- bom * fbom
      
      # Calculate "ingestion" of juvenile and adult demersals (smaller stages are pelagic) - total amount of prey consumed by demersals:
      ingestion <- f * p$Cmax * biomassend
      ingestion <- ingestion %>%
        as.data.frame() %>%
        rename(tot_ingestion = V1) %>%
        rownames_to_column(var = "stage") %>%
        filter(str_detect(stage, "demersals")) %>%
        mutate(stage = as.numeric(str_replace_all(stage, "demersals_", "")),
               id = id[i]) %>%
        filter(stage > ((1/3) * p$nstage)) %>%
        relocate(id)
      
      # Define carbon paths to demersals depending on food source:
      pelagic_spp <- c("smallZoo", "largeZoo", "smallPel", "largePel")
      benthic_spp <- c("smallBenthos", "largeBenthos")
      demersal_spp <- "demersals"
      midwater_spp <- c("mesoPel", "midwPred")
      
      # Differentiate systems with/without mesopelagic and midwPred fish:
      if (length(p$ix)==5){
        
        p$SpId <- c('smallPel','mesoPel','largePel', 'midwPred', 'demersals')
        SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
                  rep(p$SpId[1], length(p$ix[[1]])),
                  rep(p$SpId[2], length(p$ix[[2]])),
                  rep(p$SpId[3], length(p$ix[[3]])),
                  rep(p$SpId[4], length(p$ix[[4]])),
                  rep(p$SpId[5], length(p$ix[[5]])))
        
        p$RSpName <- c("Small zooplankton", "Big zooplankton", "Benthos", "Small pelagics",
                       "Mesopelagics", "Large pelagics", "Midwater predators", "Demersals")
        
        # Filter demersals and add flux pathway: 
        demers <- t(output[(p$ix[[5]][1] - length(p$ixR)):(p$ix[[5]][length(p$ix[[5]])] - length(p$ixR)), ])
        demers <- data.frame(val = c(demers), 
                             stage = rep(1:p$nstage, each = nrow(demers)), 
                             SpId = as.factor(rep(SpId, p$nstage))) %>%
          mutate(flux_pathway = case_when(SpId %in% pelagic_spp ~ "Epipelagic",
                                          SpId %in% benthic_spp ~ "Benthic",
                                          SpId %in% demersal_spp ~ "Demersal",
                                          T ~ "Midwater")) %>%
          filter(stage > ((1/3) * p$nstage))
        
        # Group by SpId and stage, keep flux pathway:
        demers <- demers %>%
          group_by(stage, SpId, flux_pathway) %>%
          summarise(feeding = sum(val))
        
      } else {
        
        # Differentiate systems with/without mesopelagic and midwater predators fish:
        p$SpId <- c('smallPel','largePel', 'demersals')
        SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
                  rep(p$SpId[1], length(p$ix[[1]])),
                  rep(p$SpId[2], length(p$ix[[2]])),
                  rep(p$SpId[3], length(p$ix[[3]])))
        
        p$RSpName <- c("Small zooplankton", "Large zooplankton", "Benthos", "Small pelagics",
                       "Large pelagics", "Demersals")   
        
        # Filter demersals and add flux pathway: 
        demers <- t(output[(p$ix[[3]][1] - length(p$ixR)):(p$ix[[3]][length(p$ix[[3]])] - length(p$ixR)), ])
        demers <- data.frame(val = c(demers), 
                             stage = rep(1:p$nstage, each = nrow(demers)), 
                             SpId = as.factor(rep(SpId, p$nstage))) %>%
          mutate(flux_pathway = case_when(SpId %in% pelagic_spp ~ "Epipelagic",
                                          SpId %in% benthic_spp ~ "Benthic",
                                          SpId %in% demersal_spp ~ "Demersal",
                                          T ~ "Midwater")) %>%
          filter(stage > ((1/3) * p$nstage))
        
        # Group by SpId and stage, keep flux pathway:
        demers <- demers %>%
          group_by(stage, SpId, flux_pathway) %>%
          summarise(feeding = sum(val))
        
      }
      
      # Add ingestion and calculate prey-speciffic ingestion:
      demers <- demers %>%
        mutate(id = id[i]) %>%
        relocate(id) %>%
        group_by(stage) %>%
        mutate(sum_feeding = sum(feeding),
               rel_feeding = feeding/sum(feeding)) %>%
        left_join(ingestion, by = c("id", "stage")) %>%
        mutate(ingestion = tot_ingestion * rel_feeding)
      
      
      out2 <- rbind(out2, demers)
    }
  })  
  
  return(list(Biomass = out,
              diet = out2))
}



#----------------------
run_feisty_sensitivity2 <- function(data) {
  # Create output files
  suppressMessages({
    out <- data.frame()
    out2 <- data.frame()
    id <- data$id
    
    #--------------------------------------------------
    # For loop for each station:
    
    for(i in 1:nrow(data)) {
      
      # Customize FEISTY setup:
      my_setup <- setupVertical2(szprod = data[i, 4], lzprod = data[i, 5], # Pelagic productivities
                                 dfbot = data[i, 6], # detrital flux reaching the bottom # <-- sensitivity *2 to compare with Clive's results and trawl surveys
                                 nStages = 9, # No. of size groups
                                 Tp = as.numeric(data[i, 9]), # Temperature of the surface
                                 Tb = as.numeric(data[i, 10]), # Temperature of the bottom
                                 Tm = as.numeric(data[i, 11]), # Temperature of the midwater
                                 depth = as.numeric(data[i, 7]), # Bottom depth
                                 photic = as.numeric(data[i, 8])) # photic depth
      
      # # Sensitivity to maturation size of demersals:
      # etaMature_new <- 0.002 # your new eatMature factor
      # my_setup$mMature[5] <- etaMature_new * 125000 # new size with 50% maturity
      # ix <- my_setup$ix[[5]] # demersal fish index in the vector
      # my_setup$psiMature[ix] <- ( 1 + (my_setup$mc[ix]/my_setup$mMature[5])^(-5) )^(-1)
      
      
      
      # Run FEISTY:
      sim <- simulateFEISTY(bCust    = T,
                            p      = my_setup, 
                            tEnd   = 200,
                            USEdll = TRUE) 
      
      #--------------------------------------------------
      # Take average biomass of each functional group of the last 80 years:
      
      totBiomass <- sim$totBiomass %>%
        as.data.frame() %>%
        tail(80) %>%
        summarise(across(everything(), mean)) %>%
        rename(smallPel = totBiomass.smallPel,
               mesoPel = totBiomass.mesoPel,
               largePel = totBiomass.largePel,
               midwPred = totBiomass.midwPred,
               demersals = totBiomass.demersals) %>%
        t() %>%
        as.data.frame() %>%
        rename(totBiomass = V1) %>%
        rownames_to_column(var = "funGroup") %>%
        mutate(propBiomass = totBiomass / sum(totBiomass)) %>%
        mutate(id = id[i])
      
      totResources <- sim$R %>%
        as.data.frame() %>%
        tail(80) %>%
        summarise(across(everything(), mean)) %>%
        mutate(Zoo = smallZoo + largeZoo) %>%
        dplyr::select(Zoo, benthos) %>%
        t() %>%
        as.data.frame() %>%
        rename(totBiomass = V1) %>%
        rownames_to_column(var = "funGroup") %>%
        mutate(propBiomass = NA) %>%
        mutate(id = id[i])
      
      out <- rbind(out, totResources, totBiomass)
      
      #--------------------------------------------------
      # Stomach content of functional groups:
      
      # Extract relevant information from the simulation:
      p=sim$p
      u=sim$u
      p$nstage <-lengths <- max(sapply(p$ix, length)) #maximum number of stages for one group
      biomass <- u
      Bin <- round(0.6 * nrow(biomass), digits = 0) 
      biomassend <- colMeans(biomass[Bin:nrow(biomass),]) # mean value of the last 40% time 
      biomassstage <- p$ixFish[length(p$ixFish)]
      biomasssmall <- p$nstage - round(2/3*p$nstage, digits = 0)
      Enc = p$V * (p$theta %*% biomassend)
      f   = Enc / (p$Cmax + Enc)
      f[is.na(f)] = 0  
      bom <- t(t(p$theta[5:biomassstage, ]) * colMeans(biomass[Bin:nrow(biomass),])) 
      fbom <- f[5:biomassstage] / rowSums(bom)
      output <- bom * fbom
      
      # Calculate "ingestion" of juvenile and adult demersals (smaller stages are pelagic) - total amount of prey consumed by demersals:
      ingestion <- f * p$Cmax * biomassend
      ingestion <- ingestion %>%
        as.data.frame() %>%
        rename(tot_ingestion = V1) %>%
        rownames_to_column(var = "stage") %>%
        filter(str_detect(stage, "demersals")) %>%
        mutate(stage = as.numeric(str_replace_all(stage, "demersals_", "")),
               id = id[i]) %>%
        filter(stage > ((1/3) * p$nstage)) %>%
        relocate(id)
      
      # Define carbon paths to demersals depending on food source:
      pelagic_spp <- c("smallZoo", "largeZoo", "smallPel", "largePel")
      benthic_spp <- c("smallBenthos", "largeBenthos")
      demersal_spp <- "demersals"
      midwater_spp <- c("mesoPel", "midwPred")
      
      # Differentiate systems with/without mesopelagic and midwPred fish:
      if (length(p$ix)==5){
        
        p$SpId <- c('smallPel','mesoPel','largePel', 'midwPred', 'demersals')
        SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
                  rep(p$SpId[1], length(p$ix[[1]])),
                  rep(p$SpId[2], length(p$ix[[2]])),
                  rep(p$SpId[3], length(p$ix[[3]])),
                  rep(p$SpId[4], length(p$ix[[4]])),
                  rep(p$SpId[5], length(p$ix[[5]])))
        
        p$RSpName <- c("Small zooplankton", "Big zooplankton", "Benthos", "Small pelagics",
                       "Mesopelagics", "Large pelagics", "Midwater predators", "Demersals")
        
        # Filter demersals and add flux pathway: 
        demers <- t(output[(p$ix[[5]][1] - length(p$ixR)):(p$ix[[5]][length(p$ix[[5]])] - length(p$ixR)), ])
        demers <- data.frame(val = c(demers), 
                             stage = rep(1:p$nstage, each = nrow(demers)), 
                             SpId = as.factor(rep(SpId, p$nstage))) %>%
          mutate(flux_pathway = case_when(SpId %in% pelagic_spp ~ "Epipelagic",
                                          SpId %in% benthic_spp ~ "Benthic",
                                          SpId %in% demersal_spp ~ "Demersal",
                                          T ~ "Midwater")) %>%
          filter(stage > ((1/3) * p$nstage))
        
        # Group by SpId and stage, keep flux pathway:
        demers <- demers %>%
          group_by(stage, SpId, flux_pathway) %>%
          summarise(feeding = sum(val))
        
      } else {
        
        # Differentiate systems with/without mesopelagic and midwater predators fish:
        p$SpId <- c('smallPel','largePel', 'demersals')
        SpId <- c("smallZoo", "largeZoo", "smallBenthos", "largeBenthos", 
                  rep(p$SpId[1], length(p$ix[[1]])),
                  rep(p$SpId[2], length(p$ix[[2]])),
                  rep(p$SpId[3], length(p$ix[[3]])))
        
        p$RSpName <- c("Small zooplankton", "Large zooplankton", "Benthos", "Small pelagics",
                       "Large pelagics", "Demersals")   
        
        # Filter demersals and add flux pathway: 
        demers <- t(output[(p$ix[[3]][1] - length(p$ixR)):(p$ix[[3]][length(p$ix[[3]])] - length(p$ixR)), ])
        demers <- data.frame(val = c(demers), 
                             stage = rep(1:p$nstage, each = nrow(demers)), 
                             SpId = as.factor(rep(SpId, p$nstage))) %>%
          mutate(flux_pathway = case_when(SpId %in% pelagic_spp ~ "Epipelagic",
                                          SpId %in% benthic_spp ~ "Benthic",
                                          SpId %in% demersal_spp ~ "Demersal",
                                          T ~ "Midwater")) %>%
          filter(stage > ((1/3) * p$nstage))
        
        # Group by SpId and stage, keep flux pathway:
        demers <- demers %>%
          group_by(stage, SpId, flux_pathway) %>%
          summarise(feeding = sum(val))
        
      }
      
      # Add ingestion and calculate prey-speciffic ingestion:
      demers <- demers %>%
        mutate(id = id[i]) %>%
        relocate(id) %>%
        group_by(stage) %>%
        mutate(sum_feeding = sum(feeding),
               rel_feeding = feeding/sum(feeding)) %>%
        left_join(ingestion, by = c("id", "stage")) %>%
        mutate(ingestion = tot_ingestion * rel_feeding)
      
      
      out2 <- rbind(out2, demers)
    }
  })  
  
  return(list(Biomass = out,
              diet = out2))
}




