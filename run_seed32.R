seed=32
nIterations=1000
thin=1
AccRateSpan=50
set.seed(seed)
load(file="data_for_model")
nYear = dim(Y)[2]
nDS = dim(Y)[3]
lili = strsplit(getwd(),split = "/")[[1]]
foldy = lili[length(lili)]
saveFileName = paste('samples_',foldy,'_seed',seed,sep="")

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
    pre_p = x[,t,]%*%matrix(Beta,dim(x)[3],1)
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


#####
# Initialize model
#####

minRho=1-(1/8.4e6)^(1/maxAge)
# Priors
Cmodel=list(
  prior_lM=function(val)log(dnorm(val,mean = log(2e5),sd=.5)),
  prior_theta=function(val)log(dunif(val,1e-5,100)),
  prior_matur=function(val)log(dpois(val-1,lambda = 6)),
  prior_d_s=function(val)log(dunif(val,0,1)),
  prior_d_l=function(val)log(dunif(val,0,1)),
  prior_lphi=function(val)log(dnorm(val,mean = log(8.4e6),sd=.7)),
  prior_Beta=function(vals)sum(log(sapply(vals,function(val)dnorm(val,mean = 0,sd = 30)))),
  prior_rho=function(val)log(dunif(val,minRho,1)),
  prior_pdetec=function(vals)sum(log(sapply(vals,function(val)dunif(val,0,1)))),
  prior_popIni=function(vals)sum(log(sapply(1:length(vals),function(i)if(initN[i]>0){dunif(vals[i],0,1000)}else if(vals[i]==0){1}else{0}))),
  prior_avgAgeRatioIni=function(vals)sum(log(dunif(vals,0,1)))
  )

# Init values
Cmodel = c(Cmodel,
           list(lM=rnorm(1,log(2e5),sd=.5),
                theta=max(rnorm(1,6,sd=2),0.5),
                matur=round(max(rnorm(1,7,sd=1.3),1.5)),
                d_s=runif(1,5e-4,1e-2),
                d_l=runif(1,5e-6,1e-4),
                lphi=log(8.4e6),
                Beta=rnorm(dim(x)[3],0,1),
                rho=minRho+.02,
                pdetec=runif(dim(Y)[3],1e-4,1e-2),
                popIni=initN*100,
                avgAgeRatioIni=rep(.2,dim(x)[1]))
)

deb = Sys.time()
Cmodel = calculate.logProb(Cmodel,allParams)
Cmodel = calculate.model(Cmodel)
print(Sys.time()-deb)

######
# Define samplers and their distribs
######

### Samplers
gSampler = list()
### Probability functions
gDistri = list()

span_lM_and_lphi=1
gSampler$lM=function(prev)rnorm(1,prev,sd=span_lM_and_lphi)
gDistri$lM=list()
gDistri$lM[[1]]=function(b,a)dnorm(b,mean = a,sd = span_lM_and_lphi)

spanTheta=8
gSampler$theta=function(prev)runif(1,prev-min(c(spanTheta,prev,100-prev+1e-3)),prev+min(c(spanTheta,prev+1e-3,100-prev)))
gDistri$theta=list()
gDistri$theta[[1]]=function(b,a)dunif(b,a-min(spanTheta,a,100-a+1e-3),a+min(spanTheta,a+1e-3,100-a))

gSampler$matur=function(prev)1+rpois(1,prev-.51)
gDistri$matur=list()
gDistri$matur[[1]]=function(b,a)dpois(b-1,lambda = a-.51)

spanD_s=.1
gSampler$d_s=function(prev)runif(1,prev-min(spanD_s,prev,1-prev+1e-10),prev+min(spanD_s,prev+1e-10,1-prev))
gDistri$d_s=list()
gDistri$d_s[[1]]=function(b,a)dunif(b,a-min(spanD_s,a,1-a+1e-10),a+min(spanD_s,a+1e-10,1-a))

gSampler$d_l=gSampler$d_s
gDistri$d_l=list()
gDistri$d_l[[1]]=gDistri$d_s[[1]]

gSampler$lphi=gSampler$lM
gDistri$lphi = list()
gDistri$lphi[[1]]=gDistri$lM[[1]]

spanBeta = 1
gSampler$Beta=function(prev){
  sapply(prev,function(pre)rnorm(1,pre,sd=spanBeta))}
gDistri$Beta = list()
for(i in 1:dim(x)[3]){gDistri$Beta[[i]]=function(b,a)dnorm(b,mean=a,sd=spanBeta)}

spanRho = .2
gSampler$rho=function(prev)runif(1,prev-min(spanRho,prev-minRho,1-prev+1e-2),prev+min(spanRho,prev-minRho+1e-2,1-prev))
gDistri$rho=list()
gDistri$rho[[1]]=function(b,a)dunif(b,a-min(spanRho,a-minRho,1-a+1e-2),a+min(spanRho,a-minRho+1e-2,1-a))

gSampler$pdetec=function(prev){
  sapply(prev,function(pre)gSampler$d_s(pre))}
gDistri$pdetec = list()
for(i in 1:dim(Y)[3]){gDistri$pdetec[[i]]=gDistri$d_s[[1]]}

spanPopIni=3*maxAge
gSampler$popIni=function(prev)sapply(1:length(prev),function(i)if(initN[i]>0){
  runif(1,max(prev[i]-spanPopIni,1),min(prev[i]+spanPopIni,1000))}else{0})
gDistri$popIni = list()
for(i in 1:dim(Y)[1]){gDistri$popIni[[i]]=function(b,a)if(a!=0){dunif(b,max(a-spanPopIni,1),min(a+spanPopIni,1000))}else if(b==0){1}else{0}}

spanAvgRatio=.2
gSampler$avgAgeRatioIni=function(prev)sapply(prev,function(pre)runif(1,pre-min(spanAvgRatio,pre,1-pre+1e-2),pre+min(spanAvgRatio,pre+1e-2,1-pre)))
gDistri$avgAgeRatioIni= list()
for(i in 1:dim(Y)[1]){gDistri$avgAgeRatioIni[[i]]=function(b,a)dunif(b,a-min(spanAvgRatio,a,1-a+1e-2),a+min(spanAvgRatio,a+1e-2,1-a))} 

######
# MCMC with Metropolis Hastings sampling 
######
### Initialization
# Compute p(y|parameters)
#log_f_a = sum(Cmodel$logProb_Y[!is.na(Cmodel$logProb_Y[])])
Cmodel = calculate.logLik(Cmodel)
log_f_a = Cmodel$logLik
# Evaluate parameters prior likelihood
for(i in 1:length(allParams)){
  log_f_a = log_f_a + Cmodel[[paste('logProb_',allParams[i],sep="")]]
}
Cmodel2 = Cmodel

cat('MCMC-MH initialized \n')

### Iterate
toSave = list()
accepto=NULL
iteration = 1
while(iteration<nIterations){
  cat('iteration ',iteration,'\n',sep = "")
  ### Sample new propositions of parameter values 
  for(i in 1:length(allParams)){
    Cmodel2[[allParams[i]]] = gSampler[[allParams[i]]](Cmodel[[allParams[i]]])
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
    #print(i)
    #cat('log_f_b:',log_f_b,'\n')
    for(j in 1:length(Cmodel2[[allParams[i]]])){
      g_b_a = gDistri[[allParams[i]]][[j]](
        Cmodel2[[allParams[i]]][j],
        Cmodel[[allParams[i]]][j])  
      if(is.na(g_b_a)){cat('Error: g_b_a of ',allParams[i],'[',j,'] is NA \n',sep = "")}
      g_a_b = gDistri[[allParams[i]]][[j]](
        Cmodel2[[allParams[i]]][j],
        Cmodel[[allParams[i]]][j])  
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
    iteration = iteration + 1
    Cmodel = Cmodel2
    log_f_a = log_f_b
    cat('candidates accepted \n')
  }else{cat('candidates rejected \n')}# Otherwise, we keep the previous value
  
  # Acceptance rate
  if(acceptationDraw<=A_b_given_a){accepto=c(accepto,T)}else{accepto=c(accepto,F)}
  if(length(accepto)>AccRateSpan){accepto=accepto[(length(accepto)-AccRateSpan):length(accepto)]}
  cat('Acceptance rate over last ',AccRateSpan,' proposals: ',sum(accepto)/length(accepto),' \n')
  
  gc(reset=T)
  # Save after each thinning cycle
  if(acceptationDraw<=A_b_given_a & iteration/thin==round(iteration/thin)){
    toSave[[round(iteration/thin)]] = lapply(1:length(allParams),function(i)Cmodel[[allParams[i]]])
    names(toSave[[round(iteration/thin)]])=allParams
    names(toSave)[round(iteration/thin)]=paste('iteraction_',iteration,sep="")
    saveRDS(object = toSave,file = saveFileName)
  }
}  

