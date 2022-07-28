require(corrplot)
require(raster)
require(data.table)
require(ggplot2)
require(ggmap)

chDir = getwd()
prefix = 'samples_pt33_seed'
validPrefix = 'samples_pt34_seed'

shinyDir = "Plectranthus_barbatus_SA_maps"
burnin = 15000
thin = 450
setwd(chDir)
seeds = c(1:5,7:10)
load(file="toAnalyse.Rdata")
load(file="data_for_model")
load(file="data_full")
load(file="validation_kit")
xm = range(ctab$x)[1];xM = range(ctab$x)[2]
xRg = c(xm-(xM-xm)/10,xM+(xM-xm)/10)
ym = range(ctab$y)[1];yM = range(ctab$y)[2]
yRg = c(ym-(yM-ym)/10,yM+(yM-ym)/10)
spatial_ext = extent(c(xRg,yRg))

mapo = get_stamenmap( bbox = c(left=spatial_ext[1],
                               bottom=spatial_ext[3],
                               right=spatial_ext[2],
                               top = spatial_ext[4]), 
                      zoom = 8, maptype = "toner-lite")

params = c('M',"theta","matur",'d_s','d_l','phi',
           'Beta','rho','pdetec','iniPop')
allParams = c('lM','theta','matur',
              'd_s','d_l','lphi',
              'Beta','rho',
              'pdetec','popIni','avgAgeRatioIni')

nDS = dim(Y)[3]
nCell = dim(x)[1]
nYear = dim(x)[2]

#####
# Functions
#####

panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE,breaks="fd")
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts;y=y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y,...)
}

minRho=1-(1/8.4e6)^(1/maxAge)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Priors
modelPriors=list(
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

# Priors
sepPriors=list(
  prior_lM=function(val)log(dnorm(val,mean = log(2e5),sd=.5)),
  prior_theta=function(val)log(dunif(val,1e-5,100)),
  prior_matur=function(val)log(dbinom(val-1,prob=.3,6)),
  prior_d_s=function(val)log(dunif(val,0,1)),
  prior_d_l=function(val)log(dunif(val,0,1)),
  prior_lphi=function(val)log(dnorm(val,mean = log(8.4e6),sd=.7)),
  prior_Beta=function(vals){
    c(log(dnorm(vals[1:3],mean = 0,sd = 30)),
    log(2*dnorm(vals[4:5],mean = 0,sd = 30)*vals[3:4]<0),
      log(dnorm(vals[6:8],mean = 0,sd = 30))) }, 
  prior_rho=function(val)log(dunif(val,minRho,1)),
  prior_pdetec=function(vals)c(log(sapply(vals,function(val)dunif(val,0,1)))),
  prior_popIni=function(vals)sum(log(sapply(1:length(vals),function(i)if(initN[i]>0){dunif(vals[i],0,1000)}else if(vals[i]==0){1}else{0}))),
  prior_avgAgeRatioIni=function(vals)sum(log(dunif(vals,0,1)))
)

prior.sample = list(
  prior_popIni=function(n)sapply(1:dim(x)[1],function(i)if(initN[i]>0){runif(n,0,1000)}else{rep(0,n)}),
  prior_avgAgeRatioIni=function(n)matrix(runif(n*dim(x)[1],0,1),n,dim(x)[1])
)

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

calculate.model = function(mod){
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

calculate.logLik = function(mod,params){
  pred = mod$E_y[]+1e-50
  mod$logLikY = sum( Y[]*log(pred) - pred - log(factorial(Y[])) ,na.rm=T)
  logL = mod$logLikY 
  for(i in 1:length(params)){
    logProbPar = mod[[paste('logProb_',params[i],sep="")]]
    logL = logL + logProbPar
  }
  mod$logLikFull = logL
  saturated_logLikY = sum(Y[]*log(Y[])-Y[]-log(factorial(Y[])),na.rm=T)
  mod$relativeLogLik = (mod$logLikFull - saturated_logLikY)
  mod$relativeLogLikY = (mod$logLikY - saturated_logLikY)
  return(mod)
}

multiplot <- function(plots=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  numPlots = length(plots)
  print(numPlots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

get.model=function(proto,params){
  model = modelPriors
  model$theta = proto$theta
  model$lM = proto$lM
  model$matur = proto$matur
  model$d_s = proto$d_s
  model$d_l = proto$d_l
  model$lphi = proto$lphi
  model$Beta = proto$Beta
  model$pdetec = proto$pdetec
  model$rho = proto$rho
  model$iniPop = proto$iniPop
  model$popIni=proto$popIni
  model$avgAgeRatioIni=proto$avgAgeRatioIni
  model = calculate.model(model)
  model = calculate.logProb(model,params)
  model = calculate.logLik(model,params)
  return(model)
}

get.dynamic = function(mod){
  nYear = dim(mod$n)[2]
  ss = colSums(mod$s)
  nn = sapply(1:nYear,function(yearo)sum(mod$n[,yearo,]))
  nt = colSums(mod$n[,,1])
  nCell = sapply(1:nYear,function(yearo)sum(rowSums(mod$n[,yearo,])>0) )
  sl = rep(NA,nYear)
  ssh = rep(NA,nYear)
  ntl = rep(NA,nYear)
  ntsh = rep(NA,nYear)
  Shanon = rep(NA,nYear)
  Disp = mod$disp_to_from
  diag(Disp) = 0
  Displ = (1 - Adj)*Disp
  Dispsh = Adj*Disp
  DisplFrom = colSums(Displ)
  DispshFrom = colSums(Dispsh)
  for(year in 1:nYear){
    sl[year] = sum(mod$s_prod[,year] * DisplFrom)
    ssh[year] = sum(mod$s_prod[,year] * DispshFrom)
    if(year>1){
      newTreesPerCellLD = mod$p[,year]*
                        mod$COMP[,year-1]*
                        colSums(mod$s_prod[,year-1] * t(Displ)) 
      ntl[year] = sum(newTreesPerCellLD)
      newTreesPerCellSD = mod$p[,year]*
                            mod$COMP[,year-1]*
                            colSums(mod$s_prod[,year-1] * t(Dispsh))
      ntsh[year] = sum(newTreesPerCellSD)
    }
    ny = rowSums(mod$n[,year,]) + 1e-20
    ny = ny/sum(ny)
    Shanon[year] = -sum(ny*log(ny))
  }
  res=data.frame(year=1979+(1:nYear),
                 pop=nn,
                 seeds=ss,
                 nCell=nCell,
                 seedsSD=ssh,
                 seedsLD=sl,
                 newTree=nt,
                 newTreeSD=ntsh,
                 newTreeLD=ntl,
                 sha=Shanon)
  return(res)
}

refineSample = function(df){
  nNeg = min(2*sum(df$y>0),sum(df$y==0))
  # We chose the minimum sampling effort threshold so that there is at least nNeg negative samples
  tgThres = max(df$tg[df$y==0])
  nSamp = sum(df$tg[df$y==0]>=tgThres)
  while(nSamp<nNeg){
    tgThres = tgThres - 1
    nSamp = sum(df$tg[df$y==0]>=tgThres)
  }
  return(df[df$y>0 | df$tg>=tgThres,])
}

auc_wmw <- function(presences, scores){
  pos <- scores[presences]
  neg <- scores[!presences]
  U <- as.numeric(wilcox.test(pos, neg)$statistic)
  U/(length(pos) * length(neg))
}

######
# Figure S4.1 - 3 Parameter sample vs prior density
######

setwd(chDir)
load(file = 'toAnalyse.Rdata')

pars = c('matur','theta','rho','lM',
  'd_s','d_l','pdetec','lphi','Beta')
setwd(chDir)
seqLen = 50
itNames=paste0('iteraction_',post$iteration)
k = 1
pList=list()
for(par in pars){
  print(par)
  if(par%in%c('pdetec','Beta')){
    sample = post[,regexpr(par,colnames(post))>0]
    priorX = sapply(1:dim(sample)[2],
      function(i)seq(min(sample[,i])-sd(sample[,i])/3,
                  max(sample[,i])+sd(sample[,i])/3,length.out=seqLen))
    colnames(priorX) = colnames(sample)
    prior=t(sapply(1:seqLen,function(j)exp(sepPriors[paste0('prior_',par)][[1]](priorX[j,]))))
    colnames(prior)=colnames(sample)
    if(par=="Beta"){
      sample=sample[,!colnames(sample)%in%paste0('Beta_',4:5)]
      prior = prior[,!colnames(prior)%in%paste0('Beta_',4:5)]
      priorX = priorX[,!colnames(priorX)%in%paste0('Beta_',4:5)]
    }
  }else{
    vec=sapply(1:dim(post)[1],function(i)chs[[as.numeric(post$chain[i])]][itNames[i]][[1]][par][[1]])
    sample = matrix(vec,dim(post)[1],1)
    colnames(sample)=par
    if(par=="matur"){
      custom = seq(0,10,1)
      priorX = matrix(custom,length(custom),1)
    }else{
      priorX = matrix(seq(min(sample[,1])-sd(sample[,1])/3,
                          max(sample[,1])+sd(sample[,1])/3,length.out=seqLen),seqLen,1)
    }
    colnames(priorX)=par
    prior = matrix(exp(modelPriors[paste('prior_',par,sep="")][[1]](priorX[,1])),
                   dim(priorX)[1],1)
    colnames(prior)=par
  }
  for(col in colnames(sample)){
    breaks = hist(sample[,col],plot=F,breaks="fd")$breaks
    if(length(breaks)>2*dim(sample)[1]){breaks=hist(sample[,col],plot=F)$breaks}
    forPrior = data.frame(x=priorX[,col],y=prior[,col])
    forPost = data.frame(x=sample[,col])
    pList[[k]] = ggplot()+
      geom_line(data=forPrior,aes(x=x,y=y),col="blue",size=1)+
      geom_histogram(data=forPost,breaks=breaks,aes(x=x,y=..density..),fill="red",alpha=.5)+
      scale_x_continuous(limits=c(min(priorX[,col]),max(priorX[,col])))+ggtitle(col)+theme_bw()
    k=k+1
  }
}

png('Z_post_priors.png',height=1000,width=666)
multiplot(pList,cols=4)
dev.off()

par=c('popIni')
k = 1
pList=list()
print(par)
sample = t(sapply(1:dim(post)[1],function(i)chs[[as.numeric(post$chain[i])]][itNames[i]][[1]][par][[1]][initN>0]))
colnames(sample)=paste0(par,'_cell_',which(initN>0))
forPrior=data.frame(x=c(.01,999.99),y=dunif(c(.01,999.99),min = 0,max = 1000))
for(col in colnames(sample)){
  breaks = hist(sample[,col],plot=F,breaks="fd")$breaks
  if(length(breaks)>2*dim(sample)[1]){breaks=hist(sample[,col],plot=F)$breaks}
  forPost = data.frame(x=sample[,col],type='posterior')
  pList[[k]] = ggplot()+
    geom_histogram(data=forPost,breaks=breaks,aes(x=x,y=..density..),alpha=.5,fill='red')+
    geom_line(data=forPrior,aes(x=x,y=y),color='blue')+
    scale_x_continuous(limits=c(0,1000))+
    ggtitle(col)+theme_bw()+theme(legend.position = 'None')
  k=k+1
}

setwd(chDir)
png('Z_post_priors_popIni.png',height=1000,width=666)
multiplot(pList,cols=4)
dev.off()

par='avgAgeRatioIni'
k = 1
pList=list()
print(par)
sample = t(sapply(1:dim(post)[1],function(i)chs[[as.numeric(post$chain[i])]][itNames[i]][[1]][par][[1]][initN>0]))
colnames(sample)=paste0(par,'_cell_',which(initN>0))
forPrior=data.frame(x=c(.01,.99),y=dunif(c(.01,.99),min = 0,max = 1))
for(col in colnames(sample)){
  breaks = hist(sample[,col],plot=F,breaks="fd")$breaks
  if(length(breaks)>2*dim(sample)[1]){breaks=hist(sample[,col],plot=F)$breaks}
  forPost = data.frame(x=sample[,col],type='posterior')
  pList[[k]] = ggplot()+
    geom_histogram(data=forPost,breaks=breaks,aes(x=x,y=..density..),alpha=.5,fill='red')+
    geom_line(data=forPrior,aes(x=x,y=y),color='blue')+
    scale_x_continuous(limits=c(0,1))+
    ggtitle(col)+theme_bw()+theme(legend.position = 'None')
  k=k+1
}

setwd(chDir)
png('Z_post_priors_avgAgeRatioIni.png',height=1000,width=666)
multiplot(pList,cols=4)
dev.off()
  
######
# Figure S4.4 Parameters correlation matrix (Identifiability)
######

setwd(chDir)
load(file = 'toAnalyse.Rdata')

#reOrd = c('matur','theta','rho','M','d_s','d_l','pdetec_1','pdetec_2','pdetec_3',
#          'phi','nTot','Beta_1','Beta_2','Beta_5','Beta_6','Beta_7')
reOrd = c('matur','theta','rho','M','d_s','d_l','pdetec_1','pdetec_2','pdetec_3',
          'phi','Beta_1','Beta_2','Beta_3','Beta_6','Beta_7','Beta_8')

parVals = post
for(col in c('theta','M','d_s','d_l','phi','rho',
          'pdetec_1','pdetec_2','pdetec_3')){
  eval(parse(text=paste('parVals$',col,'=','log10(post$',col,')',sep="")))
}

corMat=matrix(NA,length(reOrd),length(reOrd))
for(i in 1:length(reOrd)){
  for(j in 1:length(reOrd)){
    cd = complete.cases(parVals[,c(reOrd[c(i,j)])])
    corMat[i,j] = cor(parVals[cd,reOrd[i]],parVals[cd,reOrd[j]])
  }
}
rownames(corMat)=reOrd;colnames(corMat)=reOrd
corMat[is.na(corMat)]=0

png('parameter_corrplot.png',height=1000,width=1200)
print(corrplot(corMat[reOrd,reOrd],type="upper", tl.col="black", tl.srt=45, diag=F,cl.cex=1.5,cl.lim=c(-1,1),addCoef.col = "grey80",number.cex=1))
dev.off()

formu = as.formula(paste("~",paste(reOrd,collapse="+")))
png("parameter_pairs.png",height=1100,width=1100)
pairs( formu ,data=post , 
       diag.panel = panel.hist, 
       upper.panel = NULL , 
       labels= reOrd,
       pch=3 ,
       col=  rep('black',dim(post)[1]) ,
       cex.labels = 1.4 , 
       cex.axis = 1)
dev.off()

SVD=svd(corMat)
axes=SVD$d[1:2]*t(SVD$v)[1:2,,drop=F]
print(sum(SVD$d[1:2])/sum(SVD$d))
colnames(axes)=reOrd
print(axes,digits=1)

######
# Figure 7. Fecundity curves
######

setwd(chDir)
load(file = 'toAnalyse.Rdata')

for(id in 1:dim(post)[1]){
  it = post$iteration[id]
  ich = post$chain[id]
  mod = chs[[ich]][paste0('iteraction_',it)][[1]]
  M = exp(mod$lM)
  theta = mod$theta
  matur = mod$matur
  tmp = data.frame(age=1:maxAge,iteration=it,chain=ich,label=paste('ch',ich,'_it',it,sep=""),chain=as.character(ich))
  tmp$relFecundity = sapply(1:maxAge,function(k) floor( M*k^theta / ((M*(matur-1)^theta)+k^theta))/M ) 
  if(id==1){toPlot=tmp}else{toPlot=rbind(toPlot,tmp)}
  cat('\r Processed ',round(1000*id/dim(post)[1])/10,'%')
}
toPlot$label=factor(toPlot$label)
toPlot$chain = factor(toPlot$chain)

toPlot2 = aggregate(list(relFecundity=toPlot$relFecundity),
                   by=list(age=toPlot$age),mean) 
tmp = aggregate(list(minF=toPlot$relFecundity),
                by=list(age=toPlot$age),function(x)quantile(x,prob=.025,na.rm=T)) 
toPlot2 = cbind(toPlot2,tmp[,c('minF'),drop=F])
tmp = aggregate(list(maxF=toPlot$relFecundity),
                by=list(age=toPlot$age),function(x)quantile(x,prob=.975,na.rm=T)) 
toPlot2 = cbind(toPlot2,tmp[,c('maxF'),drop=F])

pFec = ggplot()+geom_line(data=toPlot2,aes(x=age,y=relFecundity))+
  geom_point(data=toPlot2,aes(x=age,y=relFecundity))+
  geom_ribbon(data=toPlot2,aes(x=age,ymin=minF,ymax=maxF),fill="grey50",alpha=.4)+
  ylab('Fecundity / M')+
  theme_bw()

if(F){pFec = ggplot(toPlot,aes(x=age,y=relFecundity,group=label,colour=chain))+geom_line(alpha=.5)+geom_point(alpha=.5)+ylab('Fecundity / M')+theme_bw()}

png('Fecundity_vs_age_relative.png')
print(pFec)
dev.off()

#####
# Figure 4. Population, cells, seeds vs time 
#####

setwd(chDir)
load(file = 'toAnalyse.Rdata')

for(id in 1:dim(post)[1]){
  flush.console()
  cat('\r sample ',id,' over ',dim(post)[1])
  it = post$iteration[id]
  ich = post$chain[id]
  modo = get.model(chs[[ich]][paste0('iteraction_',it)][[1]],allParams)
  tmp=get.dynamic(modo)
  tmp$iteration=post$iteration[id]
  tmp$chain = post$chain[id]
  tmp$label = paste('ch',ich,'_it',it,sep="")
  if(id==1){toPlot=tmp}else{toPlot=rbind(toPlot,tmp)}
}

toPlot$label=factor(toPlot$label)
toPlot$chain=factor(toPlot$chain)

setwd(chDir)
saveRDS(toPlot,file = 'pop_and_seed_vs_year.Rdata')


toPlot2 = readRDS(file='pop_and_seed_vs_year.Rdata')
toPlot = aggregate(list(pop=log10(toPlot2$pop),sha=toPlot2$sha,seeds=log10(toPlot2$seeds)),
                   by=list(year=toPlot2$year),mean) 
tmp = aggregate(list(minPop=log10(toPlot2$pop),minSha=toPlot2$sha,minSeeds=log10(toPlot2$seeds)),
                         by=list(year=toPlot2$year),function(x)quantile(x,prob=.025,na.rm=T)) 
toPlot = cbind(toPlot,tmp[,c('minPop','minSha','minSeeds')])
tmp = aggregate(list(maxPop=log10(toPlot2$pop),maxSha=toPlot2$sha,maxSeeds=log10(toPlot2$seeds)),
                  by=list(year=toPlot2$year),function(x)quantile(x,prob=.975,na.rm=T)) 
toPlot = cbind(toPlot,tmp[,c('maxPop','maxSha','maxSeeds')])

pPop = ggplot(toPlot)+
  geom_line(aes(x=year,y=pop),size=2)+
  geom_ribbon(aes(x=year,ymin=minPop,ymax=maxPop),col='grey50',alpha=.4)+
  ylab('Log10-Total population')+
  theme_bw()+theme(text=element_text(size=20))
pSha = ggplot(toPlot)+
  geom_line(aes(x=year,y=sha),size=2)+
  geom_ribbon(aes(x=year,ymin=minSha,ymax=maxSha),col='grey50',alpha=.4)+
  ylab('Shanon entropy across cells')+
  theme_bw()+theme(text=element_text(size=20))
pSee = ggplot(toPlot)+
  geom_line(aes(x=year,y=seeds),size=2)+
  geom_ribbon(aes(x=year,ymin=minSeeds,ymax=maxSeeds),col='grey50',alpha=.4)+
  ylab('Log10-Number of seeds produced')+
  theme_bw()+theme(text=element_text(size=20))

png('pop_and_seed_vs_year.png',height=1000,width=666)
multiplot(list(pPop,pSha,pSee),cols=1)
dev.off()

#####
# Figure XX. Invasion syndromes Legend
#####

tp = as.data.frame(
  expand.grid(b=seq(0,1,.05),
              r=seq(0,1,.05)))
tp$color= sapply(1:dim(tp)[1],
                 function(i)rgb(tp$r[i],0,tp$b[i]))
pts = data.frame(b=c(1,0,1,.3,0),r=c(0,0,1,1,.6))
pts$col = sapply(1:dim(pts)[1],
                 function(i)rgb(pts$r[i],0,pts$b[i]))
tp$color = factor(tp$color)
pts$col=factor(pts$col)
sizo = 1000
p=ggplot()+
  geom_tile(data=tp,aes(x=r,y=b,fill=color),alpha=.6)+
  scale_fill_manual(values=levels(tp$color))+
  geom_point(data=pts,aes(x=r,y=b,colour=col),size=round(.05*sizo))+
  scale_colour_manual(values=levels(pts$col))+
  scale_x_continuous(limits=c(-.3,1.3))+
  scale_y_continuous(limits=c(-.3,1.3))+
  theme_void()+theme(legend.position = "None")
setwd('C:/Users/user/pCloud local/boulot/data/Invasions SA et FR/pt33/')
jpeg('Figure_invasionSyndrome_legend.jpeg',height=sizo,width=sizo,quality=75)
print(p)
dev.off()


#######
# Figure 2 & 3. Main historical maps (MS)
#######

setwd(chDir)
load(file = 'toAnalyse.Rdata')
thre = .65
years = c(1980,1984,1988,1996,2008,2020)

### Plot pop growth 
statusLevels = c('uncertain',
                 'no disp & decline',
                 'no disp but growth',
                 'disp & growth',
                 'disp',
                 'disp & decline')
colors = c('gray50',
              rgb(0,0,0),# black
              rgb(0,0,1),# blue
              rgb(1,0,1),# lila
              rgb(1,0,.2),# Red
              rgb(.6,0,0))# rust
#write.table(ctab,'ctab.csv',sep=";",row.names=F,col.names=T)
pList=list()
k=1
subo=500
sampo=sort(sample(1:dim(post)[1],subo,replace = F))
for(year in years-1979){
  cat('\n year ',year+1979,'\n \n')
  tmp=as.data.frame(
    expand.grid(cell=ctab$cell[order(ctab$arrayID)],
                model=1:dim(post)[1]))
  ### Plot pop growth
  tmp$decline = NA
  tmp$growth = NA
  for(i in sampo){
    flush.console()
    cat('\r Processed ',round(100*i/dim(post)[1]),'% iterations')
    it = post$iteration[i];ich = post$chain[i]
    mod = get.model(chs[[ich]][paste0('iteraction_',it)][[1]],allParams)
    nn = rowSums(mod$n[,year,])
    nn_ = rowSums(mod$n[,year+1,])
    phi = exp(mod$lphi)
    tmp$growth[tmp$model==i] = nn_>nn 
    tmp$decline[tmp$model==i] = nn_<nn
    Disp = mod$disp_to_from
    ndie = nn[maxAge] + round(sum(mod$rho*mod$n[,year,1:(maxAge-1)]))
    diag(Disp) = 0
    to_from = mod$p[,year+1]*mod$COMP[,year]*t(mod$s_prod[,year] * t(Disp))
    tmp$spread[tmp$model==i] = colSums(to_from)>0
  }
  if(sum(is.na(tmp$growth))>0){tmp=tmp[!is.na(tmp$growth),]}
  toPlot = aggregate(list(
    decline=tmp$decline,
    growth=tmp$growth,
    locGrowth=tmp$locGrowth,
    spread=tmp$spread),by=list(cell=tmp$cell),
    FUN=function(bol)sum(bol)/subo)
  toPlot = merge(toPlot,ctab,by="cell",all.x=T)
  toPlot$status = NA
  toPlot$status[toPlot$decline>thre & toPlot$spread<=thre] = 'no disp & decline'
  toPlot$status[toPlot$growth>thre & toPlot$spread<=thre] = 'no disp but growth'
  toPlot$status[toPlot$growth>thre & toPlot$spread>thre] = 'disp & growth'
  toPlot$status[toPlot$growth<=thre & toPlot$decline<=thre & toPlot$spread>thre] = "disp"
  toPlot$status[toPlot$decline>thre & toPlot$spread>thre] = 'disp & decline'
  toPlot$status[is.na(toPlot$status)]="uncertain"
  cd = statusLevels%in%toPlot$status
  toPlot$status = factor(toPlot$status,levels=statusLevels[cd])
  
  pList[[k]]=ggmap(mapo)+geom_tile(data=toPlot,aes(x=x,y=y,fill=status),alpha=.6)+
    scale_fill_manual(values=colors[cd])+ggtitle(as.character(1979+year))+
    xlab('Longitude')+ylab('Latitude')+
    scale_x_continuous(limits=c(spatial_ext@xmin,spatial_ext@xmax))+
    scale_y_continuous(limits=c(spatial_ext@ymin,spatial_ext@ymax))+
    theme_bw()+
    theme(legend.position="None",text=element_text(size=30))
  
  if(year==years[max(1,round(length(years)/2))]-1979){
    pLeg = ggmap(mapo)+geom_tile(data=toPlot,aes(x=x,y=y,fill=status),alpha=.6)+
      scale_fill_manual(values=colors[cd])+ggtitle(as.character(1979+year))+
      xlab('Longitude')+ylab('Latitude')+theme(text=element_text(size=30))
  }
  k=k+1
}
setwd(chDir)
png('map_growth_Main.png',width=2.5*1200,height=2.5*700)
multiplot(pList,cols=2)
dev.off()

setwd(chDir)
png('map_growth_Main_Legend.png',width=2.5*1200,height=2.5*700)
print(pLeg)
dev.off()

### Plot pop state
bLabels = c(']0,q40]',']q40,q70]',']q70,q90]',']q90,q100]')
qtLevels = c(bLabels,'uncertain')
popCol = colorRampPalette(colors = c('goldenrod','darkorchid4'))(length(bLabels))
colors = c(popCol,'gray50') 

Obs = sapply(1:nYear,function(yea)rowSums(Y[,yea,])>0);Obs[!Obs]=NA;Obs[!is.na(Obs)]="record"
tgObs = sapply(1:nYear,function(yea)rowSums(TG[,yea,])>0)
obsLev = c('no record but TG records','record','not any record');obsCol = c('firebrick2','deepskyblue2','gray50')

k = 1
pList = list()
for(year in years-1979){
  tmp=as.data.frame(
    expand.grid(cell=ctab$cell[order(ctab$arrayID)],
                model=1:dim(post)[1]))
  tmp$qt = NA
  for(i in sampo){
    it = post$iteration[i]
    ich = post$chain[i]
    mod = get.model(chs[[ich]][paste0('iteraction_',it)][[1]],allParams)
    nn = sapply(1:nYear,function(yea) rowSums(mod$n[,yea,]) )
    nn = nn/max(nn,na.rm=T)
    breakos = c(-.1,quantile(as.vector(nn),prob=c(.4,.7,.9)),1)
    nf = data.frame(num=cut(nn[,year],breaks=breakos))
    leg = data.frame(num=levels(nf$num),qt=as.character(c(bLabels)))
    nf$qt = sapply(1:dim(nf)[1],function(ii)leg$qt[leg$num==nf$num[ii]])
    nf$qt=factor(nf$qt,levels=leg$qt[leg$qt%in%unique(nf$qt)])
    tmp$qt[tmp$model==i] = as.character(nf$qt)
  }
  toPlot = aggregate(list(qt=tmp$qt),
                     by=list(cell=tmp$cell),
                     FUN=function(fac){
                       prop=table(fac)/subo
                       if(max(as.numeric(prop))>=thre){
                         return(names(prop)[which.max(as.numeric(prop))[1]])
                       }else{return('uncertain')}})
  toPlot = merge(toPlot,ctab,by="cell",all.x=T)
  cd = qtLevels%in%as.character(toPlot$qt)
  toPlot$qt = factor(toPlot$qt,levels=qtLevels[cd])
  toPlot$obs = Obs[,year];toPlot$obs[is.na(toPlot$obs) & tgObs[,year]>0]='no record but TG records';toPlot$obs[is.na(toPlot$obs)]="not any record"
  cdCol=obsLev%in%toPlot$obs
  toPlot$obs=factor(toPlot$obs,levels=obsLev[cdCol])
  pList[[k]]=ggmap(mapo)+geom_tile(data=toPlot,aes(x=x,y=y,fill=qt,colour=obs),alpha=.6)+
    scale_fill_manual(values=colors[cd])+
    scale_x_continuous(limits=c(spatial_ext@xmin,spatial_ext@xmax))+
    scale_y_continuous(limits=c(spatial_ext@ymin,spatial_ext@ymax))+
    scale_colour_manual(values=obsCol[cdCol])+
    ggtitle(as.character(1979+year))+
    xlab('Longitude')+ylab('Latitude')+theme_bw()+
      theme(legend.position="None",text=element_text(size=30))
  k=k+1
}
setwd(chDir)
png('map_pop_Main.png',width=2.5*1200,height=2.5*700)
multiplot(pList,cols=2)
dev.off()


#####
# Images for Shiny app
#####

thre = .65
setwd(chDir)
if(!shinyDir%in%list.files()){
  dir.create(shinyDir)
}

# no disp & decline (marine = half blue)
# no disp but growth (blue)
# disp & growth (lila = red + blue)
# disp (magenta = red + half blue)
# disp & decline (rusty = half red)

### Plot pop growth 
statusLevels = c('uncertain',
                 'no disp & decline',
                 'no disp but growth',
                 'disp & growth',
                 'disp',
                 'disp & decline')
colorsGro = c('gray50',
              rgb(0,0,0),# black
              rgb(0,0,1),# blue
              rgb(1,0,1),# lila
              rgb(1,0,.25),# red
              rgb(.6,0,0))# rust

bLabels = c('[0,q40]',']q40,q70]',']q70,q90]',']q90,q100]')
qtLevels = c(bLabels,'uncertain')
popCol = colorRampPalette(colors = c('goldenrod','darkorchid4'))(length(bLabels))
colors = c(popCol,'gray50') 
Obs = sapply(1:nYear,function(yea)rowSums(Y[,yea,])>0);Obs[!Obs]=NA;Obs[!is.na(Obs)]="record"
tgObs = sapply(1:nYear,function(yea)rowSums(TG[,yea,])>0)
obsLev = c('no record but TG records','record','not any record');obsCol = c('firebrick2','deepskyblue2','gray50')

#write.table(ctab,'ctab.csv',sep=";",row.names=F,col.names=T)
subo=200
sampo=sort(sample(1:dim(post)[1],subo,replace = F))
for(year in 1:nYear){
  cat('\n year ',year,' \n')
  print(year)
  tmp=as.data.frame(
    expand.grid(cell=ctab$cell[order(ctab$arrayID)],
                model=sampo,
                growth=NA,decline=NA,locGrowth=NA,spread=NA))
  if(year<nYear){
    ### Plot pop growth
    for(i in sampo){
      flush.console()
      cat('\r Processed ',round(100*i/dim(post)[1]),'% iterations')
      it = post$iteration[i];ich = post$chain[i]
      mod = get.model(chs[[ich]][paste0('iteraction_',it)][[1]],allParams)
      nn = rowSums(mod$n[,year,])
      nn_ = rowSums(mod$n[,year+1,])
      phi = exp(mod$lphi)
      tmp$growth[tmp$model==i] = nn_>nn 
      tmp$decline[tmp$model==i] = nn_<nn
      Disp = mod$disp_to_from
      seeds_local = diag(Disp) * mod$s_prod[,year]
      nnew = round(seeds_local * mod$COMP[,year] * mod$p[,year+1])
      ndie = nn[maxAge] + round(sum(mod$rho*mod$n[,year,1:(maxAge-1)]))
      tmp$locGrowth[tmp$model==i] = nnew>ndie
      diag(Disp) = 0
      to_from = mod$p[,year+1]*mod$COMP[,year]*t(mod$s_prod[,year] * t(Disp))
      tmp$spread[tmp$model==i] = colSums(to_from)>0
    }
    toPlotGro = aggregate(list(
      decline=tmp$decline,
      growth=tmp$growth,
      locGrowth=tmp$locGrowth,
      spread=tmp$spread),by=list(cell=tmp$cell),
      FUN=function(bol)sum(bol)/subo)
    toPlotGro = merge(toPlotGro,ctab,by="cell",all.x=T)
    toPlotGro$status = NA
    toPlotGro$status[toPlotGro$decline>thre & toPlotGro$spread<=thre] = 'no disp & decline'
    toPlotGro$status[toPlotGro$growth>thre & toPlotGro$spread<=thre] = 'no disp but growth'
    toPlotGro$status[toPlotGro$growth>thre & toPlotGro$spread>thre] = 'disp & growth'
    toPlotGro$status[toPlotGro$growth<=thre & toPlotGro$decline<=thre & toPlotGro$spread>thre] = "disp"
    toPlotGro$status[toPlotGro$decline>thre & toPlotGro$spread>thre] = 'disp & decline'
    toPlotGro$status[is.na(toPlotGro$status)]="uncertain"
    cd = statusLevels%in%toPlotGro$status
    toPlotGro$status = factor(toPlotGro$status,levels=statusLevels[cd])
    
    pGro= ggmap(mapo)+geom_tile(data=toPlotGro,aes(x=x,y=y,fill=status),alpha=.7)+
      scale_fill_manual('growth syndrome',values=colorsGro[cd])+
      ggtitle(paste('Map of population growth syndrome in',year+1979))+
      xlab('Longitude')+ylab('Latitude')+
      theme(text=element_text(size=18))
  }else{
    pGro=list()
  }
  tmp=as.data.frame(
    expand.grid(cell=ctab$cell[order(ctab$arrayID)],
                model=sampo))
  # Plot pop
  tmp$qt = NA
  for(i in sampo){
    it = post$iteration[i]
    ich = post$chain[i]
    mod = get.model(chs[[ich]][paste0('iteraction_',it)][[1]],allParams)
    nn = sapply(1:nYear,function(yea) rowSums(mod$n[,yea,]) )
    nn = nn/max(nn,na.rm=T)
    breakos = c(-.1,quantile(as.vector(nn),prob=c(.4,.7,.9)),1)
    nf = data.frame(num=cut(nn[,year],breaks=breakos))
    leg = data.frame(num=levels(nf$num),qt=as.character(bLabels))
    nf$qt = sapply(1:dim(nf)[1],function(ii)leg$qt[leg$num==nf$num[ii]])
    nf$qt=factor(nf$qt,levels=leg$qt[leg$qt%in%unique(nf$qt)])
    tmp$qt[tmp$model==i] = as.character(nf$qt)
  }
  toPlotPop = aggregate(list(qt=tmp$qt),
                        by=list(cell=tmp$cell),
                        FUN=function(fac){
                          prop=table(fac)/subo
                          if(max(as.numeric(prop))>=thre){
                            return(names(prop)[which.max(as.numeric(prop))[1]])
                          }else{return('uncertain')}})
  toPlotPop = merge(toPlotPop,ctab,by="cell",all.x=T)
  cd = qtLevels%in%as.character(toPlotPop$qt)
  toPlotPop$qt = factor(toPlotPop$qt,levels=qtLevels[cd])
  toPlotPop$obs = Obs[,year];toPlotPop$obs[is.na(toPlotPop$obs) & tgObs[,year]>0]='no record but TG records';toPlotPop$obs[is.na(toPlotPop$obs)]="not any record"
  cdCol=obsLev%in%toPlotPop$obs
  toPlotPop$obs=factor(toPlotPop$obs,levels=obsLev[cdCol])
  
  pPop = ggmap(mapo)+geom_tile(data=toPlotPop,aes(x=x,y=y,fill=qt,colour=obs),alpha=.6,size=.3)+
    ggtitle(paste('Map of relative population in',year+1979))+
    scale_fill_manual('pop. percentile          ',values=colors[cd])+
    scale_colour_manual('data',values=obsCol[cdCol])+
    xlab('Longitude')+ylab('Latitude')+
    theme(text=element_text(size=18))
  
  setwd(paste(chDir,'/',shinyDir,'/',sep=""))
  png(paste('double_map_',year+1979,'.png',sep=""),width=1200,height=800)
  multiplot(list(pGro,pPop))
  dev.off()
  
  gc(reset=T)
}



#####
# Figure 6. Relative pop and cell decrease under dispersal ablation
#####

setwd(chDir)
load(file = 'toAnalyse.Rdata')
post$label = paste0(post$chain,'_',post$iteration)

subo = dim(post)[1]
sampo =  sample(1:dim(post)[1],subo)
metricos = as.data.frame(expand.grid(
  year=1:nYear,
  label=post$label[sampo],
  dispMode=c('all','onlyShort','onlyLong'),
  sha=NA,pop=NA))
metricos = merge(metricos,post[,c('label','iteration','chain')],by="label",all.x=T)
metricos$chain=factor(metricos$chain)
for(lab in unique(metricos$label)){
  it = unique(metricos$iteration[metricos$label==lab])
  ich = unique(metricos$chain[metricos$label==lab])
  mod = get.model(chs[[ich]][paste0('iteraction_',it)][[1]],allParams)
  dlo=mod$d_l
  dso=mod$d_s
  for(dispersal in c('all','onlyShort','onlyLong')){
    if(dispersal=='onlyLong'){
      mod$d_s=0
      mod$d_l=dlo
      mod = calculate.model(mod)
    }else if(dispersal=='onlyShort'){
      mod$d_l=0
      mod$d_s=dso
      mod = calculate.model(mod)
    }
    for(year in 1:nYear){
      cd = metricos$year==year & metricos$label==lab & metricos$dispMode==dispersal
      metricos$pop[cd] = sum(mod$n[,year,])
      ny = rowSums(mod$n[,year,]) + 1e-20
      ny = ny/sum(ny)
      metricos$sha[cd] = -sum(ny*log(ny))
    }
  }
}
tmpo = aggregate(list(pop=metricos$pop,sha=metricos$sha),by=list(year=metricos$year,dispMode=metricos$dispMode),mean)


popRatio = as.data.frame(expand.grid(year=1:nYear,label=unique(metricos$label),dispMode=c('onlyShort','onlyLong')))
popRatio$popDiff = NA
popRatio$shaDiff = NA
for(i in 1:dim(popRatio)[1]){
  yearo = popRatio$year[i]
  labo = popRatio$label[i]
  dispo = as.character(popRatio$dispMode[i])
  cdNow = metricos$year==yearo & metricos$label==labo 
  #cdBef = metricos$year==(yearo-1) & metricos$label==labo
  popRatio$popDiff[i] = metricos$pop[cdNow & metricos$dispMode==dispo] - metricos$pop[cdNow & metricos$dispMode=="all"]
  popRatio$popDiff[i] = popRatio$popDiff[i] / metricos$pop[cdNow & metricos$dispMode=="all"]
  
  popRatio$shaDiff[i] = metricos$sha[cdNow & metricos$dispMode==dispo] - metricos$sha[cdNow & metricos$dispMode=="all"]
  popRatio$shaDiff[i] = popRatio$sha[i] / metricos$sha[cdNow & metricos$dispMode=="all"]
  if(i/100==round(i/100)){
    flush.console()
    cat('\r Processed...',round(1000*i/dim(popRatio)[1])/10,'%')
  }
}
popRatio$grouping = paste(popRatio$dispMode,"_",popRatio$label,sep="")


tmp= as.data.frame(expand.grid(year=1:nYear,
                               label="MEAN",
                               dispMode=c('onlyShort','onlyLong')))
Alpha=.1
for(i in 1:dim(tmp)[1]){
  vals = popRatio$popDiff[popRatio$year==tmp$year[i] & popRatio$dispMode==tmp$dispMode[i]]
  tmp$popDiff[i] = mean(vals)
  tmp$lowerPopDiff[i] = quantile(vals,probs=Alpha/2)
  tmp$upperPopDiff[i] = quantile(vals,probs=1-Alpha/2)
  
  vals = popRatio$shaDiff[popRatio$year==tmp$year[i] & popRatio$dispMode==tmp$dispMode[i]]
  tmp$shaDiff[i] = mean(vals)
  tmp$lowerShaDiff[i] = quantile(vals,probs=Alpha/2)
  tmp$upperShaDiff[i] = quantile(vals,probs=1-Alpha/2)
}
tmp$year = tmp$year+1979

pPopDiff = ggplot(tmp,aes(x=year,y=popDiff,group=dispMode,colour=dispMode,fill=dispMode),size=1)+
  geom_line(size=1)+
  geom_ribbon(aes(ymin=lowerPopDiff,ymax=upperPopDiff),alpha=0.3,size=.3)+
  geom_point(size=2)+
  scale_y_continuous(limits=c(-1,max(tmp$upperPopDiff)))+ylab('Relative population difference (vs full model)')+
  theme_bw()+theme(text=element_text(size=20))

pShaDiff = ggplot(tmp,aes(x=year,y=shaDiff,colour=dispMode,fill=dispMode))+
  geom_line(size=1)+
  geom_ribbon(aes(ymin=lowerShaDiff,ymax=upperShaDiff),alpha=0.3,size=.3)+
  geom_point(size=2)+
  scale_y_continuous(limits=c(-1,max(tmp$upperShaDiff)))+ylab('Relative Shanon entropy difference (vs full model)')+
  theme_bw()+theme(text=element_text(size=20))

png('relative_decrease_ablation.png',width=666,height=1000)
multiplot(list(pPopDiff,pShaDiff),cols=1)
dev.off()

#####
# Figure 5. Number of new plants from long-distance dispersal vs year
#####


setwd(chDir)
load(file = 'toAnalyse.Rdata')
post$label = paste0(post$chain,'_',post$iteration)

seedSpread = as.data.frame(expand.grid(year=1:(nYear-1),label=post$label))
seedSpread$nSeeds = NA
for(lab in as.character(post$label)){
  it = post$iteration[post$label==lab]
  ich = post$chain[post$label==lab]
  mod = get.model(chs[[ich]][paste0('iteraction_',it)][[1]],allParams)
  phi = exp(mod$lphi)
  Disp = mod$disp_to_from
  diag(Disp) = 0
  Disp[Adj[]>0]=0
  prodPerYear = rep(NA,nYear-1)
  for(year in 1:(nYear-1)){
    print(year)
    to_from = mod$p[,year+1]*mod$COMP[,year]*t(mod$s_prod[,year] * t(Disp))
    prodPerYear[year]=sum(to_from)
  }
  seedSpread$nSeeds[seedSpread$label==lab] = prodPerYear
  flush.console()
  ii = which(post$label==lab)
  cat('\r Processed ',round(1000*ii/dim(post)[1])/10,'%')
}

saveRDS(seedSpread,file = "seedSpread")

alpha = .05
oneCurve = aggregate(list(mean=seedSpread$nSeeds),
                     by=list(year=seedSpread$year),function(x)mean(log10(x)))
toAdd = aggregate(list(qt2.5=seedSpread$nSeeds),
                  by=list(year=seedSpread$year),function(x)quantile(log10(x),prob=alpha/2))
toAdd2 = aggregate(list(qt97.5=seedSpread$nSeeds),
                  by=list(year=seedSpread$year),function(x)quantile(log10(x),prob=1-alpha/2))
oneCurve=cbind(oneCurve,toAdd[,'qt2.5',drop=F],toAdd2[,'qt97.5',drop=F])
oneCurve$Year=oneCurve$year+1979

p=ggplot(oneCurve)+
  geom_line(aes(x=Year,y=mean),size=2)+
  geom_ribbon(aes(x=Year,ymin=qt2.5,ymax=qt97.5),col="grey20",alpha=0.4)+
  ylab('log10-Number of new plants from L.D. dispersal')+
  theme_bw()+theme(text=element_text(size=20))

png('LDD_log10plants_vs_year.png',height=600,width=900)
print(p)
dev.off()


#####
# Figure 1. Initial Populations
#####

idIni = which(initN>0)
iniAge = as.data.frame(expand.grid(cell=idIni,
                                   label=post$label,
                                   meanAge=NA,
                                   maxAge=NA,
                                   popIni=NA))
for(lab in as.character(post$label)){
  it = post$iteration[post$label==lab]
  ich = post$chain[post$label==lab]
  flush.console()
  cat('\r sample ',which(post$label==lab),' over ',dim(post)[1])
  mod = chs[[ich]][paste0('iteraction_',it)][[1]]
  iniAge$meanAge[iniAge$label==lab] = maxAge * mod$avgAgeRatioIni[idIni]
  if(F){iniAge$maxAge[iniAge$label==lab] = sapply(idIni,function(id){
    probas = choose(maxAge,0:maxAge)*(mod$avgAgeRatioIni[id]^c(0:maxAge)) * (1-mod$avgAgeRatioIni[id])^(maxAge:0)
    max(which(round(probas*mod$popIni[id])>0),0)
  })}
  iniAge$popIni[iniAge$label==lab] = mod$popIni[idIni]
}

toPlot = aggregate(list(meanAge=iniAge$meanAge,popIni=iniAge$popIni),
                   by=list(cell=iniAge$cell),mean)
sdPop = aggregate(list(sdPop=iniAge$popIni),
                  by=list(cell=iniAge$cell),sd)
toPlot$sdPop = sdPop$sdPop
#toPlot$popIni[toPlot$sdPop/toPlot$popIni>.6] = NA

toPlot = merge(ctab,toPlot,by.x='arrayID',by.y='cell')

p2=ggmap(mapo)+geom_tile(data=toPlot,aes(x=x,y=y,fill=meanAge),alpha=.6)+
  xlab('Longitude')+ylab('Latitude')+
  scale_x_continuous(limits=c(spatial_ext@xmin,spatial_ext@xmax))+
  scale_y_continuous(limits=c(spatial_ext@ymin,spatial_ext@ymax))+
  theme_bw()+
  theme(text=element_text(size=30))

p1=ggmap(mapo)+geom_tile(data=toPlot,aes(x=x,y=y,fill=popIni),alpha=.6)+
  xlab('Longitude')+ylab('Latitude')+
  scale_fill_gradient2(mid='white',high="red")+
  scale_x_continuous(limits=c(spatial_ext@xmin,spatial_ext@xmax))+
  scale_y_continuous(limits=c(spatial_ext@ymin,spatial_ext@ymax))+
  theme_bw()+
  theme(text=element_text(size=30))

png('Figure_initial_pop.png',height=1000,width=1000)
multiplot(list(p1,p2),cols=1)
dev.off()

#####
# Figure S5.1 Environmental suitability
#####

toPlot=as.data.frame(expand.grid(id=1:dim(post)[1],
    parameter=c('Itcpt',"I(svd1)","I(svd2)","I(svd1^2)",
    "I(svd2^2)","forest","crop","urb"),
    value=NA,iteration=NA,chain=NA,label=NA))

for(id in 1:dim(post)[1]){
  cd1 = toPlot$id==id
  it = post$iteration[id]
  ich = post$chain[id]
  mod=chs[[ich]][paste0('iteraction_',it)][[1]]
  Beta = mod$Beta
  names(Beta) = c('Itcpt',dimnames(x)[[3]])
  toPlot$iteration[cd1]=it
  toPlot$chain[cd1] = ich
  toPlot$label[cd1] = paste('ch',ich,'_it',it,sep="")
  for(par in names(Beta)){
    toPlot$value[cd1 & toPlot$parameter==par]= Beta[par]
  }
}

toPlot$chain = factor(toPlot$chain)

p=ggplot(toPlot,aes(x=parameter,y=value))+
  geom_boxplot()+
  geom_point(aes(colour=chain))+theme_bw()+theme(text=element_text(size=20))

png('Figure_Beta.png',width=650,height=450)
print(p)
dev.off()

#####
# Figure S7.1 & 2 Validation
#####

load(file = 'toAnalyse_Validation.Rdata')

valid$legend = 'data deficient'
valid$legend[valid$nNoDetec>0] = 'training cell'
valid$legend[valid$valid] = "validation cell"
valid$legend = factor(valid$legend,levels=c('data deficient','training cell',"validation cell"))
p = ggmap(mapo)+
  geom_tile(data=valid,
            aes(x=x,y=y,fill=legend),alpha=.6)+
  scale_fill_manual('Cell data type',values=c("gray50",'firebrick1','chartreuse2'))+
  xlab('Longitude')+ylab('Latitude')+
  theme(text=element_text(size=20))

setwd(chDir)
png('map_validation_cells.png',width=1000,height=400)
print(p)
dev.off()

# Compute correlation prediction vs nEst (ground truth) per model
yearBreaks= c(1999,2015,2021)
volumes$period = cut(volumes$year+1979,breaks=yearBreaks)
postV$label = paste('ch',postV$chain,'_it',postV$iteration,sep="")
validCells = valid$cell[valid$valid]

metrics = as.data.frame(expand.grid(
  period=levels(volumes$period),
  label=postV$label,
  validation=c(F,T),
  cor=NA,
  auc=NA))

for(lab in as.character(postV$label)){
  it = postV$iteration[postV$label==lab]
  ich = postV$chain[postV$label==lab]
  mod = get.model(chs[[ich]][paste0('iteraction_',it)][[1]],params = allParams)
  for(per in levels(volumes$period)){
    for(valido in c(F,T)){
      if(valido){
        tmp = volumes[!is.na(volumes$period) & volumes$period==per & volumes$cell%in%valid$cell[valid$valid] & volumes$tg>0,]
      }else{
        tmp = volumes[!is.na(volumes$period) & volumes$period==per & !volumes$cell%in%valid$cell[valid$valid] & volumes$tg>0,]
      }
      tmp=refineSample(tmp)
      tmp = merge(tmp,valid[,c('cell','arrayID')],by="cell",all.x=T)
      tmp$nPred = sapply(1:dim(tmp)[1],function(j)sum(mod$n[tmp$arrayID[j],tmp$year[j],]))
      metrics$cor[metrics$label==lab & metrics$period==per & metrics$validation==valido] = cor(tmp$nPred,tmp$y/tmp$tg)
      metrics$auc[metrics$label==lab & metrics$period==per & metrics$validation==valido] = auc_wmw(1.*tmp$y>0,tmp$nPred)
    }
  }
  i=which(postV$label==lab)
  if(i/5==round(i/5)){
    flush.console()
    cat('\r Processed ',round(1000*i/dim(postV)[1])/10,'%')
  }
}

metrics = merge(metrics,postV[,c('label','chain')],by="label",all.x=T)

metrics$grouping = paste(metrics$period,metrics$validation)
metrics$grouping =sapply(1:dim(metrics)[1],function(i)if(metrics$validation[i]){paste(metrics$period[i],'valid.')}else{paste(metrics$period[i],'train')})
means = aggregate(list(auc=metrics$auc,cor=metrics$cor),
                  by=list(grouping=metrics$grouping),mean)


pValid = ggplot(metrics,aes(x=grouping,y=auc))+geom_violin(aes(fill=validation),alpha=.4)+
  geom_point(aes(group=chain),colour='grey50',position=position_dodge(width=.15))+
  geom_point(data=means,aes(x=grouping,y=auc),size=4)+
  xlab('Time period - validation/training')+
  ylab('AUC')+
  theme_bw()+theme(text=element_text(size=23))

png('auc_valid_vs_train_periods.png',height=600,width=800)
print(pValid)
dev.off()

