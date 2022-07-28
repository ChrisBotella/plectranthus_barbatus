nIterations=100000
AccRateSpan=100000

load(file="data_for_model")
load(file="data_full")

nYear = dim(Y)[2]
nDS = dim(Y)[3]
lili = strsplit(getwd(),split = "/")[[1]]
foldy = lili[length(lili)]

allParams = c('lM','theta','matur',
              'd_s','d_l','lphi',
              'Beta','rho',
              'pdetec','popIni','avgAgeRatioIni')

######
# Functions
#####

getIniPop=function(popIni,avgAgeRatioIni,maxAge){
  t(sapply(1:length(popIni),function(i){
    if(popIni[i]>0){
      tmp = popIni[i]*
        dbinom(0:(maxAge-1),
               size=maxAge-1,
               prob = avgAgeRatioIni[i])
      tmp[which.max(tmp)]=ceiling(tmp[which.max(tmp)]) # At least one
      return(round(tmp))
    }else{return(rep(0,maxAge))}
  }))
}

dlunif = function(x,a,b){1/(log(b/a)*x)}

rlunif = function(n,a,b){exp(runif(n,log(a),log(b)))}

calculate.logProb=function(momo,params){
  for(par in params){
    eval(parse(text=paste('momo$logProb_',par,'=momo$prior_',par,'(momo$',par,')',sep="")))
  }
  return(momo)
}

calculate.model=function(mod){
  theta = mod$theta
  M = exp(mod$lM)
  matur = mod$matur
  d_s = mod$d_s
  d_l = mod$d_l
  phi = exp(mod$lphi)
  Beta = mod$Beta
  pdetec = mod$pdetec
  rho = mod$rho
  mod$iniPop = getIniPop(mod$popIni,mod$avgAgeRatioIni,maxAge)
  nCell = dim(mod$iniPop)[1]
  mod$fecund = sapply(1:maxAge,function(k) floor( M*k^theta / ((M*(matur-1)^theta)+k^theta)) ) 
  mod$n = array(NA,dim = c(nCell,nYear,dim(mod$iniPop)[2]))
  mod$n[1:nCell,1,1:maxAge] <- mod$iniPop[1:nCell,1:maxAge]
  mod$dispI = 1+ (d_s-d_l*urban) * as.vector(Adj %*% matrix(1,nCell,1)) + d_l*urban*(nCell-1)
  mod$disp_to_from = t( (diag(nCell) + d_s*Adj + d_l* urban* (1-Adj)*(1-diag(nCell)))/mod$dispI )
  mod$s_prod =matrix(NA,nCell,dim(x)[2]);mod$s=mod$s_prod;mod$p=mod$s_prod;mod$COMP=mod$s_prod
  mod$E_y = array(NA,dim=c(nCell,dim(x)[2],nDS))
  # Following years
  for(t in 2:nYear){
    mod$s_prod[,t-1] = mod$n[,t-1,]%*%matrix(mod$fecund,maxAge,1)
    mod$s[,t-1] = mod$disp_to_from%*%mod$s_prod[,t-1,drop=F]
    pre_p = cbind(rep(1,dim(x)[1]),x[,1,])%*%matrix(Beta,dim(x)[3]+1,1)
    mod$p[,t] = exp(pre_p)/(1+exp(pre_p))
    survive = round((1-rho)*mod$n[,t-1,1:(maxAge-1)])
    mod$COMP[,t-1] = sapply(1:nCell,function(i){
      nOccupied = max(mod$p[i,t]*area[i]*phi-sum(survive[i,]),0)
      nSlots = max(mod$p[i,t]*area[i]*phi,mod$s[i,t-1])
      return(nOccupied/nSlots)
    })
    mod$n[,t,1] = round(mod$s[,t-1]*mod$COMP[,t-1]*mod$p[,t])
    mod$n[,t,2:maxAge] = survive
  }
  for(t in 1:nYear){
    for(i in 1:nCell){
      for(d in 1:nDS){
        mod$E_y[i,t,d] <- pdetec[d]*TG[i,t,d]*sum(mod$n[i,t,])
      }
    }
  }
  return(mod)
}

calculate.logLik = function(mod){
  pred = mod$E_y[]+1e-50
  mod$logLik = sum( Y[]*log(pred) - pred - log(factorial(Y[])))
  return(mod)
}

check.model.na = function(momo,toCheck){
  nNAs = list()
  for(par in toCheck){
    eval(parse(text=paste('nNAs[[which(par==toCheck)]] = sum(is.na(as.vector(momo$',par,')))',sep="")))
  }
  names(nNAs) = toCheck
  print(nNAs)
}

discrete.unif.dis = function(b,a,span=2,m=1,M=2){
  minBound = max(a-span,m-.4999)
  maxBound = min(a+span,M+.4999)
  if(b==round(b) & b<=round(maxBound) & b>=round(minBound)){
    bUpBound = min(maxBound,b+.5) 
    bLowBound = max(minBound,b-.5)
    return((bUpBound-bLowBound)/(maxBound-minBound))
  }else{
    return(0)
  }
}

discrete.unif.sampler = function(prev,span,m=1,M=2){
  round(runif(1,max(prev-span,m-.4999),
              min(prev+span,M+.4999)))
}

#####
# Define prior and samplers
#####

minRho=1-(1/8.4e6)^(1/maxAge)
### Priors
Cmodel=list(
  prior_lM=function(val)log(dnorm(val,mean = log(2e5),sd=.5)),
  prior_theta=function(val)log(dunif(val,1e-5,100)),
  prior_matur=function(val)log(dbinom(val-1,prob=.3,6)),
  prior_d_s=function(val)log(dunif(val,0,1)),
  prior_d_l=function(val)log(dunif(val,0,1)),
  prior_lphi=function(val)log(dnorm(val,mean = log(8.4e6),sd=.7)),
  # CUSTOMIZED FOR PLECTRANTHUS CASE
  prior_Beta=function(vals){
    sum(log(dnorm(vals[1:3],mean = 0,sd = 30)))+ # linear bioclim terms
    sum(log(2*dnorm(vals[4:5],mean = 0,sd = 30)*vals[3:4]<0))+ # quadratic bioclim terms must be negative
    sum(log(dnorm(vals[6:8],mean = 0,sd = 30)))}, # land cover rate terms 
  prior_rho=function(val)log(dunif(val,minRho,1)),
  prior_pdetec=function(vals)sum(log(sapply(vals,function(val)dunif(val,0,1)))),
  prior_popIni=function(vals)sum(log(sapply(1:length(vals),function(i)if(initN[i]>0){dunif(vals[i],0,1000)}else if(vals[i]==0){1}else{0}))),
  prior_avgAgeRatioIni=function(vals)sum(log(dunif(vals,0,1)))
)

### Samplers
gSampler = list()
### Probability functions
gDistri = list()
### Sampler scale parameters 
gScale = list()

gScale$lM=.7
gSampler$lM=function(prev,sca)rnorm(1,prev,sd=sca)
gDistri$lM=list()
gDistri$lM[[1]]=function(b,a,sca)dnorm(b,mean=a,sd=sca)

gScale$theta=1
gSampler$theta=function(prev,sca)runif(1,prev-min(c(sca,prev)),prev+min(c(sca,100-prev)))
gDistri$theta=list()
gDistri$theta[[1]]=function(b,a,sca)dunif(b,a-min(sca,a),a+min(sca,100-a))

gScale$matur=1
gSampler$matur=function(prev,sca){discrete.unif.sampler(prev,span=sca,m=1,M=7)}
gDistri$matur=list()
gDistri$matur[[1]]=function(b,a,sca)discrete.unif.dis(b,a,span=sca,m=1,M=7)

gScale$d_s = .02
gSampler$d_s=function(prev,sca){
  MaxSca= .99 * min(prev,1-prev);sca = min(MaxSca,sca)
  expo = exp(2*sca/prev)
  b = 2*sca*expo/(expo-1);a = b-2*sca
  return(rlunif(1,a,b))
}

gDistri$d_s=list()
gDistri$d_s[[1]]=function(b,a,sca){
  MaxSca= .99 * min(a,1-a);sca = min(MaxSca,sca)
  expo = exp(2*sca/a)
  M = 2*sca*expo/(expo-1);m = M-2*sca
  return(dlunif(b,a = m,b = M))
}

gScale$d_l = .01
gSampler$d_l=gSampler$d_s
gDistri$d_l=list()
gDistri$d_l[[1]]=gDistri$d_s[[1]]

gScale$lphi = .7
gSampler$lphi=gSampler$lM
gDistri$lphi = list()
gDistri$lphi[[1]]=gDistri$lM[[1]]

# CUSTOMIZED FOR PLECTRANTHUS CASE
gScale$Beta = .7
gSampler$Beta=function(prev,sca){
  c(sapply(prev[1:3],function(pre)rnorm(1,pre,sd=sca)),
    sapply(prev[4:5],function(pre)runif(1,pre+max(-sca,pre),pre-max(-sca,pre))),
    sapply(prev[6:8],function(pre)rnorm(1,pre,sd=sca)))
  }
gDistri$Beta = list()
for(i in 1:(dim(x)[3]+1)){
  if(i%in%c(1:3,6:8)){
    gDistri$Beta[[i]]=function(b,a,sca)dnorm(b,mean=a,sd=sca)
  }else{
    gDistri$Beta[[i]]=function(b,a,sca){
      if(a!=b | b!=0){
        dunif(b,a+max(-sca,a),a-max(-sca,a))
      }else{1}}
  }
}

gScale$rho=.04
gSampler$rho=function(prev,sca)runif(1,prev-min(sca,prev-minRho,1-prev+1e-2),prev+min(sca,prev-minRho+1e-2,1-prev))

gDistri$rho=list()
gDistri$rho[[1]]=function(b,a,sca)dunif(b,a-min(sca,a-minRho,1-a+1e-2),a+min(sca,a-minRho+1e-2,1-a))

gScale$pdetec = .01
gSampler$pdetec=function(prev,sca){
  sapply(prev,function(pre)gSampler$d_s(pre,sca))}
gDistri$pdetec = list()
for(i in 1:dim(Y)[3]){gDistri$pdetec[[i]]=gDistri$d_s[[1]]}

gScale$popIni=maxAge
gSampler$popIni=function(prev,sca)sapply(1:length(prev),function(i)if(initN[i]>0){
  runif(1,max(prev[i]-sca,1),min(prev[i]+sca,1000))}else{0})
gDistri$popIni = list()
for(i in 1:dim(Y)[1]){gDistri$popIni[[i]]=function(b,a,sca)if(a!=0){dunif(b,max(a-sca,1),min(a+sca,1000))}else if(b==0){1}else{0}}

gScale$avgAgeRatioIni=.05
gSampler$avgAgeRatioIni=function(prev,sca)sapply(prev,function(pre)runif(1,pre-min(sca,pre),pre+min(sca,1-pre)))
gDistri$avgAgeRatioIni= list()
for(i in 1:dim(Y)[1]){gDistri$avgAgeRatioIni[[i]]=function(b,a,sca)dunif(b,a-min(sca,a),a+min(sca,1-a))} 

######
# MCMC with Metropolis Hastings sampling 
######

popInito = rep(0,length(initN))
popInito[initN>0]=c(123,183,97,184,538,358,72,205,332,228,551,587,274,73,108,443,546,523,131,90,75,385,124,81)
avgAgeRatioInito = rep(0.1,length(initN))
avgAgeRatioInito[initN>0]=c(0.101254,0.100395,0.154849,0.237527,0.689562,0.731742,0.468439,0.607656,0.33309,0.876167,0.028901,0.361551,0.196713,0.241745,0.278403,0.384715,0.280448,0.638099,0.394353,0.73292,0.656702,0.167059,0.921603,0.303481)

for(seed in 1:10){
  # Run several chains starting from different random seeds
  set.seed(seed)
  saveFileName = paste('samples_',foldy,'_seed',seed,sep="")

  # Initial values
  Cmodel = c(Cmodel,
             list(lM=13.03496,
                  theta=6.647526,
                  matur=3,
                  d_s=1.623841e-06,
                  d_l=1.625358e-07,
                  lphi=17.50612,
                  Beta=c(0.7688628,-0.03306837,-0.07731294,-1e-300,-1e-300,
                         10.79504,-1.505942,6.881799),
                  rho=0.4937154,
                  pdetec=c(5.958113e-11, 2.178828e-11,2.368315e-10),
                  popIni=popInito,
                  avgAgeRatioIni=avgAgeRatioInito)
  )
  
  deb = Sys.time()
  Cmodel = calculate.logProb(Cmodel,allParams)
  Cmodel = calculate.model(Cmodel)
  # Compute p(y|parameters)
  Cmodel = calculate.logLik(Cmodel)
  
  log_f_a = Cmodel$logLik
  # Evaluate parameters prior likelihood
  for(i in 1:length(allParams)){
    log_f_a = log_f_a + Cmodel[[paste('logProb_',allParams[i],sep="")]]
  }
  cat('posterior logLikelihood:',log_f_a,' \n',sep = "")
  
  Cmodel2 = Cmodel
  
  cat('\n MCMC-MH initialized for seed ',seed,' \n')
  
  ### Start MCMC
  toSave = list()
  accepto=NULL
  iteration = 1
  likDf=data.frame(iteration=NA,logLik=NA)
  likDf=likDf[-1,,drop=F]
  while(iteration<nIterations){
    cat('iteration ',iteration,'\n',sep = "")
    ### Sample new propositions of parameter values 
    for(i in 1:length(allParams)){
      Cmodel2[[allParams[i]]] = gSampler[[allParams[i]]](Cmodel[[allParams[i]]],sca=gScale[[allParams[i]]])
    }
    ### Update model state with candidate parameters
    # Compute candidate parameters log probabilities 
    Cmodel2 = calculate.logProb(Cmodel2,allParams)
    # Compute all deterministic nodes given candidate parameters
    Cmodel2 = calculate.model(Cmodel2)
    # Initialize candidate logLikelhood p(y,parameters) with logLikelhood p(y|candidate parameters)
    Cmodel2 = calculate.logLik(Cmodel2)
    log_f_b = Cmodel2$logLik
    # Compute acceptance probability of all candidate together 
    # by multiplying over parameters
    A_b_given_a = 1
    for(i in 1:length(allParams)){
      logProbPar = Cmodel2[[paste('logProb_',allParams[i],sep="")]]
      log_f_b = log_f_b + logProbPar
      for(j in 1:length(Cmodel2[[allParams[i]]])){
        g_b_a = gDistri[[allParams[i]]][[j]](
          Cmodel2[[allParams[i]]][j],
          Cmodel[[allParams[i]]][j],
          sca = gScale[[allParams[i]]])  
        if(is.na(g_b_a)){cat('Error: g_b_a of ',allParams[i],'[',j,'] is NA \n',sep = "")}
        g_a_b = gDistri[[allParams[i]]][[j]](
          Cmodel2[[allParams[i]]][j],
          Cmodel[[allParams[i]]][j],
          sca = gScale[[allParams[i]]])  
        if(is.na(g_a_b)){cat('Error: g_a_b of ',allParams[i],'[',j,'] is NA \n',sep = "")}
        if((!is.infinite(g_b_a) | !is.infinite(g_a_b)) & (g_b_a!=0 | g_a_b!=0)){
          A_b_given_a = A_b_given_a*(g_b_a/g_a_b)
        }
        if(is.na(A_b_given_a)){
          cat('Error: A_b_given_a of ',allParams[i],'[',j,'] is NA \n',sep = "")
          cat('g_b_a:',g_b_a,' \n')
          cat('g_a_b:',g_a_b,' \n')
          cat('logProbPar:',logProbPar,' \n')
        } 
      }
    }
    cat('log_f_b:',log_f_b,' \n',sep = "")
    # Finally multiply by f_b/f_a = p(y|candidate parameters) / p(y|previous parameters)
    A_b_given_a = min( A_b_given_a* exp(log_f_b - log_f_a) , 1)
    cat('Acceptance proba=',A_b_given_a,' ...',sep = "")
    acceptationDraw = runif(1,0,1)
    if(acceptationDraw<=A_b_given_a){
      # If the candidate parameters set is accepted (probability A_b_given_a) 
      # we update Cmodel with this one
      Cmodel = Cmodel2
      log_f_a = log_f_b
      cat('candidates accepted \n')
    }else{cat('candidates rejected \n')}# Otherwise, we keep the previous value
    likDf = rbind(likDf,data.frame(iteration=iteration,logLik=log_f_a))
    
    # Acceptance rate
    if(acceptationDraw<=A_b_given_a){accepto=c(accepto,T)}else{accepto=c(accepto,F)}
    if(length(accepto)>AccRateSpan){accepto=accepto[(length(accepto)-AccRateSpan):length(accepto)]}
    cat('Acceptance rate over last ',AccRateSpan,' proposals: ',sum(accepto)/length(accepto),' \n')
    
    gc(reset=T)
    # Save after each new iteration
    toSave[[iteration]] = lapply(1:length(allParams),function(i)Cmodel[[allParams[i]]])
    names(toSave[[iteration]])=allParams
    names(toSave)[iteration]=paste('iteraction_',iteration,sep="")
    if(iteration/30==round(iteration/30)){
      saveRDS(object = toSave,file = saveFileName)
      saveRDS(object = likDf,file = paste(saveFileName,'_logLik',sep=""))
    }
    iteration = iteration + 1
  }  
  
}



