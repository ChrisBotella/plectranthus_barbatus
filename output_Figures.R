require(corrplot)
require(raster)
require(data.table)
require(ggplot2)
require(ggmap)
require(shiny)
require(rsconnect)

chDir = getwd()
prefix = 'samples_pt21_seed'

shinyDir = "Plectranthus_barbatus_SA_maps"
minLogLik = -542.18
burnin = 20
thin = 5
setwd(chDir)
seeds = c(32,34,38,40)
chs = list()
for(seed in seeds){
  chs[[which(seeds==seed)]]=readRDS(paste(prefix,seed,sep=""))
}

load(file="data_for_model")

params = c('M',"theta","matur",'d_s','d_l','phi',
           'Beta','rho','pdetec','iniPop')
allParams = c('lM','theta','matur',
              'd_s','d_l','lphi',
              'Beta','rho',
              'pdetec','popIni','avgAgeRatioIni')

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

# Priors
modelPriors=list(
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

prior.sample = list(
  prior_lM=function(n)rnorm(n,log(2e5),sd=.5),
  prior_theta=function(val)runif(val,1e-5,100),
  prior_matur=function(val)1+rpois(val,lambda = 6),
  prior_d_s=function(val)runif(val,0,1),
  prior_d_l=function(val)dunif(val,0,1),
  prior_lphi=function(val)rnorm(val,mean = log(8.4e6),sd=.7),
  prior_Beta=function(n)sapply(1:length(models[[1]][[1]]$Beta),function(i)rnorm(n,mean = 0,sd = 30)),
  prior_rho=function(val)runif(val,minRho,1),
  prior_pdetec=function(n)sapply(1:length(models[[1]][[1]]$pdetec),function(i)runif(n,0,1)),
  prior_popIni=function(n)sapply(1:length(models[[1]][[1]]$popIni),function(i)if(initN[i]>0){runif(n,0,1000)}else{rep(0,n)}),
  prior_avgAgeRatioIni=function(n)matrix(runif(n*length(models[[1]][[1]]$avgAgeRatioIni),0,1),n,length(models[[1]][[1]]$avgAgeRatioIni))
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

calculate.logProb=function(momo,params){
  for(par in params){
    eval(parse(text=paste('momo$logProb_',par,'=momo$prior_',par,'(momo$',par,')',sep="")))
  }
  return(momo)
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

get.models.chain.list = function(chain,params,its){
  models = list()
  k= 1
  for(j in its){
    flush.console()
    cat('\r iteration ',j-1)
    models[[k]] = modelPriors
    models[[k]]$theta = chain[[j]]$theta
    models[[k]]$lM = chain[[j]]$lM
    models[[k]]$matur = chain[[j]]$matur
    models[[k]]$d_s = chain[[j]]$d_s
    models[[k]]$d_l = chain[[j]]$d_l
    models[[k]]$lphi = chain[[j]]$lphi
    models[[k]]$Beta = chain[[j]]$Beta
    models[[k]]$pdetec = chain[[j]]$pdetec
    models[[k]]$rho = chain[[j]]$rho
    models[[k]]$iniPop = chain[[j]]$iniPop
    models[[k]]$popIni=chain[[j]]$popIni
    models[[k]]$avgAgeRatioIni=chain[[j]]$avgAgeRatioIni
    models[[k]] = calculate.model(models[[k]])
    models[[k]] = calculate.logProb(models[[k]],params)
    models[[k]] = calculate.logLik(models[[k]],params)
    k=k+1
  }
  names(models)=its
  return(models)
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
  Disp = mod$disp_to_from
  diag(Disp) = 0
  Displ = (1 - Adj)*Disp
  Dispsh = Adj*Disp
  DisplFrom = colSums(Displ)
  DispshFrom = colSums(Dispsh)
  for(year in 1:nYear){
    sl[year] = round(sum(mod$s_prod[,year] * DisplFrom))
    ssh[year] = round(sum(mod$s_prod[,year] * DispshFrom))
    if(year>1){
      newTreesPerCellLD = round(mod$p[,year]*
                        mod$COMP[,year-1]*
                        colSums(mod$s_prod[,year-1] * t(Displ)) )
      ntl[year] = sum(newTreesPerCellLD)
      newTreesPerCellSD = round(mod$p[,year]*
                            mod$COMP[,year-1]*
                            colSums(mod$s_prod[,year-1] * t(Dispsh)))
      ntsh[year] = sum(newTreesPerCellSD)
    }
  }
  res=data.frame(year=1979+(1:nYear),
                 pop=nn,
                 seeds=ss,
                 nCell=nCell,
                 seedsSD=ssh,
                 seedsLD=sl,
                 newTree=nt,
                 newTreeSD=ntsh,
                 newTreeLD=ntl)
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

#####
# Compute quantities
#####

firsto = 2
nDS = length(chs[[1]][[firsto]]$pdetec)
nCell = dim(x)[1]
nYear = dim(x)[2]

models = list()
for(i in 1:length(chs)){
  cat('\n chain ',i,'\n')
  its = seq(length(chs[[i]])-thin*((length(chs[[i]])-burnin-1)%/%thin),length(chs[[i]]),by=thin)
  models[[i]] = 
    get.models.chain.list(
      chs[[i]],
      allParams,
      its=its)
}


######
# Parameters vs iterations
######

#ich = 1
#i = length(models[[ich]])
#models[[ich]][[i]]$avgAgeRatioIni[initN>0]
#models[[ich]][[i]]$popIni[initN>0]
#sum(models[[ich]][[i]]$avgAgeRatioIni[initN>0]<1e-3)


for(ich in 1:length(chs)){
  its = as.numeric(names(models[[ich]]))
  tmp = data.frame(iteration=its,chain=ich)
  tmp$totIniPop = sapply(1:length(models[[ich]]),function(i)sum(models[[ich]][[i]]$popIni))
  tmp$maxAgeIni = sapply(1:length(models[[ich]]),function(i) max(which(colSums(getIniPop(models[[ich]][[i]]$popIni,models[[ich]][[i]]$avgAgeRatioIni,maxAge))>0))  )  
  tmp$M = sapply(1:length(models[[ich]]),function(i)exp(models[[ich]][[i]]$lM))
  tmp$phi = sapply(1:length(models[[ich]]),function(i)exp(models[[ich]][[i]]$lphi))
  for(j in 1:nDS){
    tmp[,paste('pdetec_',j,sep="")] = 
      sapply(1:length(models[[ich]]),
        function(i)models[[ich]][[i]]$pdetec[j])}
  for(j in 1:length(chs[[1]][[2]]$Beta)){
    tmp[,paste('Beta_',j,sep="")] =
      sapply(1:length(models[[ich]]),
             function(i)models[[ich]][[i]]$Beta[j])
  }
  for(col in c('matur','theta','d_s','d_l','rho','logLikFull','logLikY','relativeLogLikY','lphi','lM')){
    eval(parse(text=
                 paste('tmp$',col,'=sapply(1:length(models[[ich]]),function(i)models[[ich]][[i]]$',col,')',sep="")))
  }
  # Summary statistics
  tmp$mean_p = NA
  tmp$nTot = NA
  for(i in 1:length(models[[ich]])){
    tmp$mean_p[i] = mean(models[[ich]][[i]]$p,na.rm=T)
    tmp$nTot[i] = sum(as.vector(models[[ich]][[i]]$n[,,1]))+sum(as.vector(models[[ich]][[i]]$iniPop))
  }
  
  range(models[[ich]][[i]]$avgAgeRatioIni)
  
  if(ich==1){post=tmp}else{post=rbind(post,tmp)}
}

head(post,digits=2)

#####
# Samples cutoff
#####
postKept = post[post$logLikFull>minLogLik,]
cat('\n',dim(postKept)[1],' kept models for min logLikelihood: ',minLogLik)
postKept$label= paste('ch',postKept$chain,'_it',postKept$iteration,sep="")

######
# Figure S3.2 Parameter sample vs prior density
######


pars = c('matur','theta','rho','lM',
  'd_s','d_l','pdetec','lphi','Beta')
#pars=c('popIni','avgAgeRatioIni')
setwd(chDir)
seqLen = 50
itNames=paste('iteraction_',postKept$iteration,sep="")
k = 1
pList=list()
for(par in pars){
  print(par)
  if(par%in%c('pdetec','Beta')){
    sample = t(sapply(1:dim(postKept)[1],function(i)chs[[as.numeric(postKept$chain[i])]][itNames[i]][[1]][par][[1]]))
    colnames(sample)=paste(par,'_',1:dim(sample)[2],sep="")
    
    priorX = sapply(1:dim(sample)[2],function(i)seq(min(sample[,i])-sd(sample[,i])/3,
                                                    max(sample[,i])+sd(sample[,i])/3,length.out=seqLen))
    colnames(priorX) = colnames(sample)
    prior=sapply(1:dim(sample)[2],function(i){
      priorVals = NULL
      for(j in 1:seqLen){
        priorVals[j] = exp(modelPriors[paste('prior_',par,sep="")][[1]](priorX[j,i]))
      }
      return(priorVals)
      })
    colnames(prior)=colnames(sample)
  }else if(par%in%c('popIni','avgAgeRatioIni')){
    sample = t(sapply(1:dim(postKept)[1],function(i)chs[[as.numeric(postKept$chain[i])]][itNames[i]][[1]][par][[1]][initN>0]))
    colnames(sample)=paste(par,'_',1:dim(sample)[2],sep="")
    prior=prior.sample[paste('prior_',par,sep="")][[1]](nSamp)[,initN>0]
    colnames(prior)=colnames(sample)
  }else{
    vec=sapply(1:dim(postKept)[1],function(i)chs[[as.numeric(postKept$chain[i])]][itNames[i]][[1]][par][[1]])
    sample = matrix(vec,dim(postKept)[1],1)
    colnames(sample)=par
    if(par=="matur"){
      custom = seq(0,6,.05)
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
    pList[[k]] = ggplot()+geom_line(data=forPrior,aes(x=x,y=y),col="blue",size=1)+
      geom_histogram(data=forPost,breaks=breaks,aes(x=x,y=..density..),fill="red",alpha=.5)+
      scale_x_continuous(limits=c(min(priorX[,col]),max(priorX[,col])))+ggtitle(col)+theme_bw()
    k=k+1
  }
}

png('Z_post_priors.png',height=1000,width=666)
multiplot(pList,cols=4)
dev.off()

######
# Figure S3.3 Parameters correlation matrix (Identifiability)
######

#reOrd = c('theta','matur','M','d_s','d_l','phi','rho','mean_p',
#          'pdetec_1','pdetec_2','pdetec_3','nTot')
reOrd = c('matur','theta','rho','M','d_s','d_l','pdetec_1','pdetec_2','pdetec_3',
          'phi','nTot')

parVals = postKept
for(col in c('theta','M','d_s','d_l','phi','rho','mean_p',
          'pdetec_1','pdetec_2','pdetec_3','nTot')){
  eval(parse(text=paste('parVals$',col,'=','log10(postKept$',col,')',sep="")))
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

png('parameter_corrplot.png',height=800,width=900)
print(corrplot(corMat[reOrd,reOrd],type="upper", tl.col="black", tl.srt=45, diag=F,cl.cex=1.5,cl.lim=c(-1,1),addCoef.col = "grey80",number.cex=1))
dev.off()

formu = as.formula(paste("~",paste(reOrd,collapse="+")))
png("parameter_pairs.png",height=800,width=800)
pairs( formu ,data=postKept , 
       diag.panel = panel.hist, 
       upper.panel = NULL , 
       labels= reOrd,
       pch=3 ,
       col=  rep('black',dim(postKept)[1]) ,
       cex.labels = 1.4 , 
       cex.axis = 1)
dev.off()



######
# Figure 7. Fecundity curves
######

for(id in 1:dim(postKept)[1]){
  it = postKept$iteration[id]
  ich = postKept$chain[id]
  mod = models[[ich]][as.character(it)][[1]]
  M = exp(mod$lM)
  theta = mod$theta
  matur = mod$matur
  tmp = data.frame(age=1:maxAge,iteration=it,chain=ich,label=paste('ch',ich,'_it',it,sep=""),chain=as.character(ich))
  tmp$relFecundity = sapply(1:maxAge,function(k) floor( M*k^theta / ((M*(matur-1)^theta)+k^theta)) ) /M
  if(id==1){toPlot=tmp}else{toPlot=rbind(toPlot,tmp)}
}
toPlot$label=factor(toPlot$label)
toPlot$chain = factor(toPlot$chain)
pFec = ggplot(toPlot,aes(x=age,y=relFecundity,group=label,colour=chain))+geom_line(alpha=.5)+geom_point(alpha=.5)+ylab('Fecundity / M')+scale_y_continuous(limits=c(0,1))+theme_bw()

png('Fecundity_vs_age.png')
print(pFec)
dev.off()

#####
# Figure 5. Population, cells, seeds vs time 
#####
setwd(chDir)

for(id in 1:dim(postKept)[1]){
  it = postKept$iteration[id]
  ich = postKept$chain[id]
  tmp=get.dynamic(models[[ich]][as.character(it)][[1]])
  tmp$iteration=postKept$iteration[id]
  tmp$chain = postKept$chain[id]
  tmp$label = paste('ch',ich,'_it',it,sep="")
  if(id==1){toPlot=tmp}else{toPlot=rbind(toPlot,tmp)}
}

toPlot$label=factor(toPlot$label)
toPlot$chain=factor(toPlot$chain)
pPop = ggplot(toPlot,aes(x=year,y=log10(pop),group=label,colour=chain))+geom_line(alpha=.6)+ylab('log10-total population')+theme_bw()+theme(legend.position = "None")
pCell = ggplot(toPlot,aes(x=year,y=nCell,group=label,colour=chain))+geom_line(alpha=.6)+ylab('Number of cells colonized')+theme_bw()
pSee = ggplot(toPlot,aes(x=year,y=log10(seeds),group=label,colour=chain))+geom_line(alpha=.6)+ylab('log10-seeds produced')+theme_bw()+theme(legend.position = "None")
png('pop_and_seed_vs_year.png',height=1000,width=666)
multiplot(list(pPop,pCell,pSee),cols=1)
dev.off()


#######
# Figure 3 & 4. Main historical maps (MS)
#######

thre = .65
# Color code for population status:
# Darkorchid: Certainly absent
# Grey : Not certainly absent but not certainly increasing either
# Blue : Certainly increasing 
# Green : Certainly increasing only with local recruitment
# Yellow : Certainly increasing and disseminate at short distance
# Red : Certainly increasing and disseminate at short and long distance
setwd(chDir)
load(file = "preliminary_data")
mapo = get_stamenmap( bbox = c(left=spatial_ext[1],
                               bottom=spatial_ext[3]+1,
                               right=spatial_ext[2],
                               top = spatial_ext[4]), 
                      zoom = 8, maptype = "toner-lite")
rRef = crop(rRef,spatial_ext)
cellsCov1 = rRef[]
coos = as.data.frame(rasterToPoints(rRef))
area = areas[areas$cell%in%cellsCov1,]
cellsCov2 = area$cell[area$propLand>0]
ctab = data.frame(arrayID=1:length(cellsCov2),cell=cellsCov2[order(cellsCov2)])
coos = coos[coos[,3]%in%ctab$cell,]
ctab = merge(ctab,coos,by.x="cell",by.y="layer",all.x=T)
ctab = merge(ctab,area[area$propLand>0,],by="cell",all.x=T)
urbo = aggregate(list(urb=areasFull$urb),by=list(cell=areasFull$cell),mean)
ctab = merge(ctab,urbo,by="cell",all.x=T)
load(file="data_for_model")

years = c(1980,1981,1984,1994,2005,2020)

### Plot pop growth
statusLevels = c('absence','uncertain','growth','intrinsic growth','spread')
colors = c('white','gray50','royalblue2','springgreen3','firebrick1')
#write.table(ctab,'ctab.csv',sep=";",row.names=F,col.names=T)
k=1
for(year in years-1979){
  tmp=as.data.frame(
    expand.grid(cell=ctab$cell[order(ctab$arrayID)],
                model=1:dim(postKept)[1]))
  ### Plot pop growth
  tmp$isAbs = NA
  tmp$Growth = NA
  for(i in 1:dim(postKept)[1]){
    it = postKept$iteration[i];ich = postKept$chain[i]
    mod = models[[ich]][as.character(it)][[1]]
    nn = rowSums(mod$n[,year,])
    nn_ = rowSums(mod$n[,year+1,])
    phi = exp(mod$lphi)
    tmp$isAbs[tmp$model==i] = nn==0
    tmp$growth[tmp$model==i] = nn_>nn 
    Disp = mod$disp_to_from
    seeds_local = diag(Disp) * mod$s_prod[,year]
    survive = round((1-mod$rho)*mod$n[,year,1:(maxAge-1)])
    COMP = sapply(1:nCell,function(j){
      nAvail = max(mod$p[j,year+1]*area[j]*phi-sum(survive[j,]),0)
      nSlots = max(mod$p[j,year+1]*area[j]*phi,mod$s[j,year])
      return(nAvail/nSlots)
    })
    nnew = round(seeds_local * COMP * mod$p[,year+1])
    ndie = nn[maxAge] + sum(round(mod$rho*mod$n[,year,1:(maxAge-1)]))
    tmp$locGrowth[tmp$model==i] = nnew>ndie
    diag(Disp) = 0
    to_from = round(mod$p[,year+1]*COMP*t(mod$s_prod[,year] * t(Disp)))
    tmp$spread[tmp$model==i] = colSums(to_from)>0
  }
  toPlot = aggregate(list(
    isAbs=tmp$isAbs,
    growth=tmp$growth,
    locGrowth=tmp$locGrowth,
    spread=tmp$spread),by=list(cell=tmp$cell),
    FUN=function(bol)sum(bol)/dim(postKept)[1])
  toPlot = merge(toPlot,ctab,by="cell",all.x=T)
  
  toPlot$status = NA
  toPlot$status[toPlot$growth>thre] = "growth"
  toPlot$status[toPlot$locGrowth>thre] = "intrinsic growth"
  toPlot$status[toPlot$spread>thre] = "spread"
  toPlot$status[is.na(toPlot$status) & toPlot$isAbs>thre] = "absence"
  toPlot$status[is.na(toPlot$status)]="uncertain"
  cd = statusLevels%in%toPlot$status
  toPlot$status = factor(toPlot$status,levels=statusLevels[cd])
  
  pList[[k]]=ggmap(mapo)+geom_tile(data=toPlot,aes(x=x,y=y,fill=status),alpha=.6)+
    scale_fill_manual(values=colors[cd])+ggtitle(as.character(1979+year))+
    xlab('Longitude')+ylab('Latitude')+theme(legend.position="None",text=element_text(size=30))
  
  if(year+1979==1982){
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
qtLevels = c('{0}',bLabels,'uncertain')
popCol = colorRampPalette(colors = c('goldenrod','darkorchid4'))(length(bLabels))
colors = c('white',popCol,'gray50') 

Obs = sapply(1:nYear,function(yea)rowSums(Y[,yea,])>0);Obs[!Obs]=NA;Obs[!is.na(Obs)]="record"
tgObs = sapply(1:nYear,function(yea)rowSums(TG[,yea,])>0)
obsLev = c('no record but TG records','record','not any record');obsCol = c('firebrick2','deepskyblue2','gray50')

k = 1
pList = list()
for(year in years-1979){
  tmp=as.data.frame(
    expand.grid(cell=ctab$cell[order(ctab$arrayID)],
                model=1:dim(postKept)[1]))
  tmp$qt = NA
  for(i in 1:dim(postKept)[1]){
    it = postKept$iteration[i]
    ich = postKept$chain[i]
    mod = models[[ich]][as.character(it)][[1]]
    nn = sapply(1:nYear,function(yea) rowSums(mod$n[,yea,]) )
    nn = nn/max(nn,na.rm=T)
    breakos = c(-.1,0,quantile(as.vector(nn),prob=c(.4,.7,.9)),1)
    nf = data.frame(num=cut(nn[,year],breaks=breakos))
    leg = data.frame(num=levels(nf$num),qt=as.character(c('{0}',bLabels)))
    nf$qt = sapply(1:dim(nf)[1],function(ii)leg$qt[leg$num==nf$num[ii]])
    nf$qt=factor(nf$qt,levels=leg$qt[leg$qt%in%unique(nf$qt)])
    tmp$qt[tmp$model==i] = as.character(nf$qt)
  }
  toPlot = aggregate(list(qt=tmp$qt),
                     by=list(cell=tmp$cell),
                     FUN=function(fac){
                       prop=table(fac)/dim(postKept)[1]
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

if(!shinyDir%in%list.files()){
  dir.create(shinyDir)
}

### Plot pop growth 
statusLevels = c('absence','uncertain','growth','intrinsic growth','spread')
colorsGro = c('white','gray50','royalblue2','springgreen3','firebrick1')


bLabels = c(']0,q40]',']q40,q70]',']q70,q90]',']q90,q100]')
qtLevels = c('{0}',bLabels,'uncertain')
popCol = colorRampPalette(colors = c('goldenrod','darkorchid4'))(length(bLabels))
colors = c('white',popCol,'gray50') 
Obs = sapply(1:nYear,function(yea)rowSums(Y[,yea,])>0);Obs[!Obs]=NA;Obs[!is.na(Obs)]="record"
tgObs = sapply(1:nYear,function(yea)rowSums(TG[,yea,])>0)
obsLev = c('no record but TG records','record','not any record');obsCol = c('firebrick2','deepskyblue2','gray50')

#write.table(ctab,'ctab.csv',sep=";",row.names=F,col.names=T)
for(year in 1:nYear){
  print(year)
  tmp=as.data.frame(
    expand.grid(cell=ctab$cell[order(ctab$arrayID)],
                model=1:dim(postKept)[1]))
  if(year<nYear){
    ### Plot pop growth
    tmp$isAbs = NA
    tmp$Growth = NA
    for(i in 1:dim(postKept)[1]){
      it = postKept$iteration[i];ich = postKept$chain[i]
      mod = models[[ich]][as.character(it)][[1]]
      nn = rowSums(mod$n[,year,])
      nn_ = rowSums(mod$n[,year+1,])
      phi = exp(mod$lphi)
      tmp$isAbs[tmp$model==i] = nn==0
      tmp$growth[tmp$model==i] = nn_>nn 
      Disp = mod$disp_to_from
      seeds_local = diag(Disp) * mod$s_prod[,year]
      survive = round((1-mod$rho)*mod$n[,year,1:(maxAge-1)])
      COMP = sapply(1:nCell,function(j){
        nAvail = max(mod$p[j,year+1]*area[j]*phi-sum(survive[j,]),0)
        nSlots = max(mod$p[j,year+1]*area[j]*phi,mod$s[j,year])
        return(nAvail/nSlots)
      })
      nnew = round(seeds_local * COMP * mod$p[,year+1])
      ndie = nn[maxAge] + sum(round(mod$rho*mod$n[,year,1:(maxAge-1)]))
      tmp$locGrowth[tmp$model==i] = nnew>ndie
      diag(Disp) = 0
      to_from = round(mod$p[,year+1]*COMP*t(mod$s_prod[,year] * t(Disp)))
      tmp$spread[tmp$model==i] = colSums(to_from)>0
    }
    toPlotGro = aggregate(list(
      isAbs=tmp$isAbs,
      growth=tmp$growth,
      locGrowth=tmp$locGrowth,
      spread=tmp$spread),by=list(cell=tmp$cell),
      FUN=function(bol)sum(bol)/dim(postKept)[1])
    toPlotGro = merge(toPlotGro,ctab,by="cell",all.x=T)
    toPlotGro$status = NA
    toPlotGro$status[toPlotGro$growth>thre] = "growth"
    toPlotGro$status[toPlotGro$locGrowth>thre] = "intrinsic growth"
    toPlotGro$status[toPlotGro$spread>thre] = "spread"
    toPlotGro$status[is.na(toPlotGro$status) & toPlotGro$isAbs>thre] = "absence"
    toPlotGro$status[is.na(toPlotGro$status)]="uncertain"
    cd = statusLevels%in%toPlotGro$status
    toPlotGro$status = factor(toPlotGro$status,levels=statusLevels[cd])
    
    pGro= ggmap(mapo)+geom_tile(data=toPlotGro,aes(x=x,y=y,fill=status),alpha=.6)+
      scale_fill_manual('growth syndrome',values=colorsGro[cd])+
      ggtitle(paste('Map of population growth syndrome in',year+1979))+
      xlab('Longitude')+ylab('Latitude')+
      theme(text=element_text(size=18))
  }else{
    pGro=list()
  }
  tmp=as.data.frame(
    expand.grid(cell=ctab$cell[order(ctab$arrayID)],
                model=1:dim(postKept)[1]))
  # Plot pop
  tmp$qt = NA
  for(i in 1:dim(postKept)[1]){
    it = postKept$iteration[i]
    ich = postKept$chain[i]
    mod = models[[ich]][as.character(it)][[1]]
    nn = sapply(1:nYear,function(yea) rowSums(mod$n[,yea,]) )
    nn = nn/max(nn,na.rm=T)
    breakos = c(-.1,0,quantile(as.vector(nn),prob=c(.4,.7,.9)),1)
    nf = data.frame(num=cut(nn[,year],breaks=breakos))
    leg = data.frame(num=levels(nf$num),qt=as.character(c('{0}',bLabels)))
    nf$qt = sapply(1:dim(nf)[1],function(ii)leg$qt[leg$num==nf$num[ii]])
    nf$qt=factor(nf$qt,levels=leg$qt[leg$qt%in%unique(nf$qt)])
    tmp$qt[tmp$model==i] = as.character(nf$qt)
  }
  toPlotPop = aggregate(list(qt=tmp$qt),
                        by=list(cell=tmp$cell),
                        FUN=function(fac){
                          prop=table(fac)/dim(postKept)[1]
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
}

#####
# Figure 6. Relative pop and cell decrease under dispersal ablation
#####

metricos = as.data.frame(expand.grid(
  year=1:nYear,
  label=postKept$label,
  dispMode=c('all','onlyShort','onlyLong'),
  nCell=NA,pop=NA))
metricos = merge(metricos,postKept[,c('label','iteration','chain')],by="label",all.x=T)
metricos$chain=factor(metricos$chain)
for(lab in unique(metricos$label)){
  it = unique(metricos$iteration[metricos$label==lab])
  ich = unique(metricos$chain[metricos$label==lab])
  mod = models[[ich]][as.character(it)][[1]]
  for(dispersal in c('all','onlyShort','onlyLong')){
    if(dispersal=='onlyLong'){
      mod$d_s=0
      mod$d_l = models[[ich]][as.character(it)][[1]]$d_l
      mod = calculate.model(mod)
    }else if(dispersal=='onlyShort'){
      mod$d_s = models[[ich]][as.character(it)][[1]]$d_s
      mod$d_l=0
      mod = calculate.model(mod)
    }
    for(year in 1:nYear){
      cd = metricos$year==year & metricos$label==lab & metricos$dispMode==dispersal
      metricos$pop[cd] = sum(mod$n[,year,])
      metricos$nCell[cd] = sum(rowSums(mod$n[,year,])>0)
    }
  }
}

popRatio = as.data.frame(expand.grid(year=1:nYear,label=unique(metricos$label),dispMode=c('onlyShort','onlyLong')))
popRatio$popDiff = NA
popRatio$cellDiff = NA
for(i in 1:dim(popRatio)[1]){
  yearo = popRatio$year[i]
  labo = popRatio$label[i]
  dispo = as.character(popRatio$dispMode[i])
  cdNow = metricos$year==yearo & metricos$label==labo 
  #cdBef = metricos$year==(yearo-1) & metricos$label==labo
  popRatio$popDiff[i] = metricos$pop[cdNow & metricos$dispMode==dispo] - metricos$pop[cdNow & metricos$dispMode=="all"]
  popRatio$popDiff[i] = popRatio$popDiff[i] / metricos$pop[cdNow & metricos$dispMode=="all"]
  
  popRatio$cellDiff[i] = metricos$nCell[cdNow & metricos$dispMode==dispo] - metricos$nCell[cdNow & metricos$dispMode=="all"]
  popRatio$cellDiff[i] = popRatio$cellDiff[i] / metricos$nCell[cdNow & metricos$dispMode=="all"]
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
  
  vals = popRatio$cellDiff[popRatio$year==tmp$year[i] & popRatio$dispMode==tmp$dispMode[i]]
  tmp$cellDiff[i] = mean(vals)
  tmp$lowerCellDiff[i] = quantile(vals,probs=Alpha/2)
  tmp$upperCellDiff[i] = quantile(vals,probs=1-Alpha/2)
}
tmp$year = tmp$year+1979

pPopDiff = ggplot(tmp,aes(x=year,y=popDiff,group=dispMode,colour=dispMode,fill=dispMode),size=1)+
  geom_line(size=1)+
  geom_ribbon(aes(ymin=lowerPopDiff,ymax=upperPopDiff),alpha=0.3,size=.3)+
  geom_point(size=2)+
  scale_y_continuous(limits=c(-1,max(tmp$upperPopDiff)))+ylab('Relative population difference (vs full model)')+
  theme_bw()+theme(text=element_text(size=20))

pCellDiff = ggplot(tmp,aes(x=year,y=cellDiff,colour=dispMode,fill=dispMode))+
  geom_line(size=1)+
  geom_ribbon(aes(ymin=lowerCellDiff,ymax=upperCellDiff),alpha=0.3,size=.3)+
  geom_point(size=2)+
  scale_y_continuous(limits=c(-1,max(tmp$upperCellDiff)))+ylab('Relative n° of cells difference (vs full model)')+
  theme_bw()+theme(text=element_text(size=20))


png('relative_decrease_ablation.png',width=666,height=1000)
multiplot(list(pPopDiff,pCellDiff),cols=1)
dev.off()

#####
# Figure S4.4 Environmental suitability
#####

load(file = "preliminary_data")
# Parameters
rRef = crop(rRef,spatial_ext)
cellsCov1 = rRef[]
coos = as.data.frame(rasterToPoints(rRef))
area = areas[areas$cell%in%cellsCov1,]
cellsCov2 = area$cell[area$propLand>0]
x1 = tab[tab$cell%in%cellsCov2,]
cd = !colnames(x1)%in%c('cell','year')
x1[,cd]=scale(x1[,cd])
SVD=svd(x1[,cd])
components = t(SVD$v)[1:2,]
colnames(components) = colnames(x1[,cd])
rownames(components) = paste('svd',1:2,sep="")
print(components,digits=3)
ctab = data.frame(arrayID=1:length(cellsCov2),cell=cellsCov2[order(cellsCov2)])
coos = coos[coos[,3]%in%ctab$cell,]
ctab = merge(ctab,coos,by.x="cell",by.y="layer",all.x=T)
ctab = merge(ctab,area[area$propLand>0,],by="cell",all.x=T)
load(file="data_for_model")
#print(sum(SVD$d[1:2])/sum(SVD$d))
#SVD=SVD$u[,1:2,drop=F]

toPlot=as.data.frame(expand.grid(id=1:dim(postKept)[1],
    parameter=c("I(svd1)","I(svd2)","I(svd1^2)",
    "I(svd2^2)","forest","crop","urb"),
    value=NA,iteration=NA,chain=NA,label=NA))

for(id in 1:dim(postKept)[1]){
  cd1 = toPlot$id==id
  it = postKept$iteration[id]
  ich = postKept$chain[id]
  mod=models[[ich]][as.character(it)][[1]]
  Beta = mod$Beta
  names(Beta) = dimnames(x)[[3]]
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

png('Beta.png',width=650,height=450)
print(p)
dev.off()

parSdTab = data.frame(id=1)
rownames(parSdTab)='Std. deviation'
for(par in as.character(unique(toPlot$parameter))){
  eval(parse(text=paste('parSdTab[,"',par,'"]=',sd(x[,,par]),sep="")))
}
print(parSdTab[,-1],digits=3)

if(F){
  xScale=seq(min(x[,,'I(svd1)']),max(x[,,'I(svd1)']),length.out=50)
  yScale=seq(min(x[,,'I(svd2)']),max(x[,,'I(svd2)']),length.out=50)
  toPlot = data.frame(svd1=1,svd2=1,val=1,label=NA);toPlot=toPlot[-1,,drop=F]
  enviVar = c('I(svd1)','I(svd2)','I(svd1^2)','I(svd2^2)')
  for(id in 1:dim(postKept)[1]){
    cd1 = toPlot$id==id
    it = postKept$iteration[id]
    ich = postKept$chain[id]
    mod=models[[ich]][as.character(it)][[1]]
    Beta = mod$Beta
    names(Beta) = dimnames(x)[[3]]
    Beta = matrix(Beta[enviVar],length(enviVar),1)
    tmp1 = as.data.frame(expand.grid(
      svd1 = xScale,svd2 = yScale))
    tmp = model.matrix(object = as.formula(paste('~-1+',paste(enviVar,collapse="+"))),data=tmp1)
    tmpVal = as.vector(tmp %*% Beta)
    tmp1$val = exp(tmpVal)/(1+exp(tmpVal))
    tmp1$label = postKept$label[id]
    toPlot=rbind(toPlot,tmp1)
  }
  toPlot=merge(toPlot,postKept[,c('label','chain')],by="label",all.x=T)
  toPlotTmp = toPlot[toPlot$chain%in%c(1,2,4,3),,drop=F]
  toPloto = aggregate(list(val=toPlotTmp$val),
                      by=list(svd1=toPlotTmp$svd1,svd2=toPlotTmp$svd2),
                      mean)
  plotM = matrix(toPloto$val,
                 length(unique(toPlot$svd1)),
                 length(unique(toPlot$svd2)))
  contour(unique(toPlot$svd1),unique(toPlot$svd2),plotM,
          nlevels=8,lwd=2,labcex=1,
          xlab="svd1",ylab="svd2")
}

#####
# Figure S6.6 Validation
#####

load(file="validation_kit")

valid$legend = 'data deficient'
valid$legend[valid$nNoDetec>0] = 'training cell'
valid$legend[valid$valid] = "validation cell"
valid$legend = factor(valid$legend,levels=c('data deficient','training cell',"validation cell"))
p = ggmap(mapo)+geom_tile(data=valid,aes(x=x,y=y,fill=legend),alpha=.6)+scale_fill_manual('Cell data type',values=c("gray50",'firebrick1','chartreuse2'))+xlab('Longitude')+ylab('Latitude')+theme(text=element_text(size=20))

png('map_validation_cells.png',width=1000,height=400)
print(p)
dev.off()

# Compute correlation prediction vs nEst (ground truth) per model

yearBreaks= c(1999,2015,2021)

volumes$period = cut(volumes$year+1979,breaks=yearBreaks)

postKept$label = paste('ch',postKept$chain,'_it',postKept$iteration,sep="")

validCells = valid$cell[valid$valid]

metrics = as.data.frame(expand.grid(
  period=levels(volumes$period),
  label=postKept$label,
  validation=c(F,T),
  cor=NA,
  auc=NA))
for(i in 1:dim(metrics)[1]){
  per = metrics$period[i]
  lab = metrics$label[i]
  it = postKept$iteration[postKept$label==lab]
  ich = postKept$chain[postKept$label==lab]
  mod = models[[ich]][as.character(it)][[1]]
  
  if(metrics$validation[i]){
    tmp = volumes[!is.na(volumes$period) & volumes$period==per & volumes$cell%in%valid$cell[valid$valid] & volumes$tg>0,]
  }else{
    tmp = volumes[!is.na(volumes$period) & volumes$period==per & !volumes$cell%in%valid$cell[valid$valid] & volumes$tg>0,]
  }
  tmp=refineSample(tmp)
  tmp = merge(tmp,valid[,c('cell','arrayID')],by="cell",all.x=T)
  tmp$nPred = sapply(1:dim(tmp)[1],function(j)sum(mod$n[tmp$arrayID[j],tmp$year[j],]))
  metrics$cor[i] = cor(tmp$nPred,tmp$y/tmp$tg)
  metrics$auc[i] = auc_wmw(tmp$y>0,tmp$nPred)
  #plot(tmp$nPred,tmp$y/TestTmp$tg,ylab="n detections / n TG records",xlab="Predicted population")
  if(i/5==round(i/5)){
    flush.console()
    cat('\r Processed ',round(1000*i/dim(metrics)[1])/10,'%')
  }
}

metrics = merge(metrics,postKept[,c('label','chain')],by="label",all.x=T)

metrics$grouping = paste(metrics$period,metrics$validation)
metrics$grouping =sapply(1:dim(metrics)[1],function(i)if(metrics$validation[i]){paste(metrics$period[i],'valid.')}else{paste(metrics$period[i],'train')})
pValid = ggplot(metrics,aes(x=grouping,y=auc))+geom_violin(aes(fill=validation),alpha=.4)+
  geom_point(aes(group=chain,color=chain), position=position_dodge(width=.15))+
  xlab('Time period - validation/training')+
  ylab('AUC')+
  theme_bw()+theme(text=element_text(size=23))

png('auc_valid_vs_train_periods.png',height=600,width=800)
print(pValid)
dev.off()

#?wilcox.test()$statistic

