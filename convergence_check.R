require(coda)
require(data.table)
require(ggplot2)

chDir = getwd()
prefix = 'samples_pt33_seed'

setwd(chDir)
seeds = c(1:9)
chs = list()
gc(reset=T)
for(seed in seeds){
  print(seed)
  chs[[which(seeds==seed)]]=try(readRDS(paste0(prefix,seed)))
  if(is.character(chs[[which(seeds==seed)]])){
    chs[[which(seeds==seed)]]=NULL
  }
}

# Get likelihoods 
lls = list()
for(seed in seeds){
  lls[[which(seeds==seed)]]=readRDS(paste(prefix,seed,'_logLik',sep=""))
}

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

maxAge=50
minRho=1-(1/8.4e6)^(1/maxAge)

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

#####
# Plot posterior Likelihood (Figure 1 of App S3)
#####

allLik= data.frame(seed=NA,iteration=NA,logLik=NA)
allLik=allLik[-1,,drop=F]
for(ich in 1:length(lls)){
  tmp = cbind(seed=seeds[ich],lls[[ich]])
  allLik = rbind(allLik,tmp)
}

allLik$seed=factor(allLik$seed)
p=ggplot(allLik,aes(x=iteration,y=logLik,group=seed,colour=seed))+geom_line(size=1.5)+
  xlab('MCMC Iterations')+ylab('Log-likelihood')+scale_y_continuous(limits=c(-525,-450))+theme_bw()+
  theme(text=element_text(size=50))

if(!'trace plots'%in%list.files()){dir.create('trace plots')}

png('trace plots/logLik_allSeeds.png',width=2000,height=1000)
print(p)
dev.off()

p=ggplot(allLik,aes(x=iteration,y=logLik,group=seed,colour=seed))+geom_line(size=1.5)+
  xlab('MCMC Iterations')+ylab('Log-likelihood')+scale_y_continuous(limits=c(-650,-590))+theme_bw()+
  theme(text=element_text(size=50))

png('trace plots/logLik_allSeeds_ZOOM.png',width=2000,height=1000)
print(p,cols=3)
dev.off()

# overall maximum loglikelihood
maxLik = max(allLik$logLik)
# Chain where maxLik obtained
maxCh = allLik$seed[allLik$logLik==max(allLik$logLik)][1]
# First iteration of chain maxCh where maxLik obtained
firstItMaxLik = min(allLik$iteration[allLik$seed==maxCh & allLik$logLik==maxLik])
# Minimum logLikelihood obtained in maxCh after passing by maxLik
minLikAfterMax = min(allLik$logLik[allLik$seed==maxCh & allLik$iteration>firstItMaxLik])
chToKeep=data.frame(seed=unique(allLik$seed),
                    maxLik=NA,
                    nIt=NA)
for(seed in unique(allLik$seed)){
  cd = chToKeep$seed==seed
  chToKeep$nIt[cd]=max(allLik$iteration[allLik$seed==seed])
  maxLikSeed = max(allLik$logLik[allLik$seed==seed])
  chToKeep$maxLik[cd]=maxLikSeed
}
print(chToKeep)

######
# Parameters vs iterations
######

for(ich in 1:length(chs)){
  its = as.numeric(sapply(strsplit(names(chs[[ich]]),split = "_"),function(el)el[2]))
  tmp = data.frame(iteration=its[1:length(its)],chain=ich,seed=seeds[ich])
  tmp$lM = sapply(1:length(chs[[ich]]),function(i)chs[[ich]][[i]]$lM)
  tmp$M = exp(tmp$lM)
  tmp$lphi = sapply(1:length(chs[[ich]]),function(i)chs[[ich]][[i]]$lphi)
  tmp$phi = exp(tmp$lphi)
  for(j in 1:length(chs[[1]][[2]]$pdetec)){
    tmp[,paste('pdetec_',j,sep="")] = 
      sapply(1:length(chs[[ich]]),
             function(i)chs[[ich]][[i]]$pdetec[j])}
  for(j in 1:length(chs[[1]][[2]]$Beta)){
    tmp[,paste('Beta_',j,sep="")] =
      sapply(1:length(chs[[ich]]),
             function(i)chs[[ich]][[i]]$Beta[j])
  }
  for(col in c('matur','theta','d_s','d_l','rho')){
    eval(parse(text=
                 paste('tmp$',col,'=sapply(1:length(chs[[ich]]),function(i)chs[[ich]][[i]]$',col,')',sep="")))
  }
  tmp$logLikelihood = lls[[ich]]$logLik[1:dim(tmp)[1]]
  if(ich==1){post=tmp}else{post=rbind(post,tmp)}
}
head(post,digits=2)

for(ich in 1:length(chs)){
  cat('seed ',seedsKept[ich],'\n')
  llVals= unique(post$logLikelihood[post$chain==ich])
  llVals=llVals[!is.na(llVals)]
  jumps  = sapply(llVals,function(ll){
    if(ll==llVals[1]){
      0
    }else{
      min(post$iteration[!is.na(post$logLikelihood) & post$logLikelihood==ll & post$chain==ich])}})
  itsBeforeJump = jumps-c(0,jumps[1:(length(jumps)-1)])
  print(itsBeforeJump)
}

length(unique(paste(post$logLikelihood,post$seed,sep="_")))

thini=2000
burnin = 30000
postT = post[round((post$iteration-1)/thini)==(post$iteration-1)/thini & post$iteration>burnin,]
cat('n samples ',dim(postT)[1],'\n')
cat('n ',length(unique(paste(postT$logLikelihood,postT$seed,sep="_"))))

#####
# Trace plots
#####

cd = postT$iteration<=max(postT$iteration)
sizeText = 30
prefix="allChains_"
postT$seed=factor(postT$seed)

cols = c('logLikelihood','lM','theta','matur',
         'd_s','d_l','rho','lphi',
         paste('Beta_',1:length(chs[[1]][[2]]$Beta),sep=""),
         paste('pdetec_',1:length(chs[[1]][[2]]$pdetec),sep=""))
  
pList=lapply(cols,function(paro){
  if(paro%in%c('d_s','d_l') | regexpr('pdetec',paro)>0){
    return(ggplot(postT[cd,c(paro,'seed','iteration')],
                  aes(x=iteration,y=log10(postT[cd,paro]),
                      group=seed,color=seed))+scale_y_continuous(limits=c(c(max(min(log10(postT[cd,paro]),-10)),max(log10(postT[cd,paro])))))+
             geom_line(size=1)+xlab('iteration')+ylab(paste('log10-',paro))+
             theme(text=element_text(size=sizeText)))
  }else if(paro=='logLikelihood'){
    return(ggplot(postT[cd,c(paro,'seed','iteration')],
                                      aes(x=iteration,y=postT[cd,paro],
                                          group=seed,color=seed))+
             scale_y_continuous(limits=c(c(-650,-600)))+
      geom_line(size=1)+xlab('iteration')+ylab(paste('log10-',paro))+
      theme(text=element_text(size=sizeText)))
  }else{
    return(ggplot(postT[cd,c(paro,'seed','iteration')],
                                      aes(x=iteration,y=postT[cd,paro],group=seed,color=seed))+
      geom_line(size=1)+xlab('iteration')+ylab(paro)+
      theme(text=element_text(size=sizeText)))
  }
})

png(paste('trace plots/',prefix,'allParams.png',sep=""),width=2000,height=3000)
multiplot(pList,cols=3)
dev.off()


popInis = as.data.frame(matrix(NA, dim(postT)[1] , sum(chs[[1]][[2]]$popIni>0) ))
colnames(popInis)=as.character(which(chs[[1]][[2]]$popIni>0))
for(j in 1:dim(popInis)[2]){
  popInis[,j]= sapply(1:dim(popInis)[1],function(i)chs[[postT$chain[i]]][[postT$iteration[i]]]$popIni[as.numeric(colnames(popInis)[j])])
}
popInis$iteration = postT$iteration
popInis$seed = postT$seed
maxPop = max(as.vector(popInis[,1:24]))
pList=lapply(colnames(popInis)[1:(length(colnames(popInis))-3)],function(paro){
  ggplot(popInis[,c(paro,'seed','iteration')],
                  aes(x=iteration,y=popInis[,paro],
                      group=seed,color=seed))+
    geom_line(size=1)+xlab('iteration')+
    ylab(paste('Initial pop of cell',paro))+
    scale_y_continuous(limits=c(0,maxPop))+
    theme(text=element_text(size=sizeText))
  })

png(paste('trace plots/',prefix,'_popIni.png',sep=""),width=2000,height=3000)
multiplot(pList,cols=4)
dev.off()


avgRat = as.data.frame(matrix(NA, dim(postT)[1] , sum(chs[[1]][[2]]$popIni>0) ))
colnames(avgRat)=as.character(which(chs[[1]][[2]]$popIni>0))
for(j in 1:dim(avgRat)[2]){
  avgRat[,j]= sapply(1:dim(avgRat)[1],function(i)chs[[postT$chain[i]]][[postT$iteration[i]]]$avgAgeRatioIni[as.numeric(colnames(avgRat)[j])])
}
avgRat$iteration = postT$iteration
avgRat$seed = postT$seed
pList=lapply(colnames(avgRat)[1:(length(colnames(avgRat))-3)],function(paro){
  ggplot(avgRat[,c(paro,'seed','iteration')],
         aes(x=iteration,y=avgRat[,paro],
             group=seed,color=seed))+
    geom_line(size=1)+xlab('iteration')+
    ylab(paste('Initial pop of cell',paro))+
    scale_y_continuous(limits=c(0,1))+
    theme(text=element_text(size=sizeText))
})

png(paste('trace plots/',prefix,'_avgRatioIni.png',sep=""),width=2000,height=3000)
multiplot(pList,cols=4)
dev.off()


######
# MCMC convergence diagnosis
######

cols = c('theta','matur','d_s','d_l','rho',
         'lphi','lM',
         paste('Beta_',1:3,sep=""),
         paste('Beta_',6:length(chs[[1]][[2]]$Beta),sep=""),
         paste('pdetec_',1:length(chs[[1]][[2]]$pdetec),sep=""))

maxIt=min(sapply(1:length(chs),function(ich)max(postT$iteration[postT$chain==ich])))
# List of mcmc objects
mcmcs=list()
for(ich in 1:length(chs)){
  itos = unique(postT$iteration[postT$chain==ich & postT$iteration<=maxIt])
  mat = matrix(NA,length(itos),length(cols))
  colnames(mat)=cols 
  it = 1
  for(it in 1:length(itos)){
    for(paro in cols){
      if(paro%in%c('d_s','d_l') | regexpr('pdetec',paro)>0){
        mat[it,paro] = log(postT[postT$iteration==itos[it] & postT$chain==ich,paro])
      }else if(paro=="matur"){
        mat[it,paro] = postT[postT$iteration==itos[it] & postT$chain==ich,paro]+rnorm(1,0,.05)
      }else{
        mat[it,paro] = postT[postT$iteration==itos[it] & postT$chain==ich,paro]
      }
    }
    it=it+1
  }
  mcmcs[[ich]] = mcmc(data = mat)
}
mcmcs=as.mcmc.list(mcmcs)

### Gelman and Rubin's univariate + multivariate convergence tests
# Code taken directly from gelman.diag function (coda) which sends error 
x=mcmcs
confidence=.95
Niter <- niter(x)
Nchain <- nchain(x)
Nvar <- nvar(x)
xnames <- varnames(x)
x <- lapply(x, as.matrix)
S2 <- array(sapply(x, var, simplify = TRUE), dim = c(Nvar, 
                                                     Nvar, Nchain))
W <- apply(S2, c(1, 2), mean)
xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE), 
               nrow = Nvar, ncol = Nchain)
B <- Niter * var(t(xbar))
CW <- chol(W)
emax <- eigen(backsolve(CW, t(backsolve(CW, B, transpose = TRUE)), 
                            transpose = TRUE), symmetric = TRUE, only.values = TRUE)$values[1]    
mpsrf <- sqrt((1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter)
w <- diag(W)
b <- diag(B)
s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain)
muhat <- apply(xbar, 1, mean)
var.w <- apply(s2, 1, var)/Nchain
var.b <- (2 * b^2)/(Nchain - 1)
cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) - 2 * 
                                  muhat * var(t(s2), t(xbar)))
V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter
var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b + 
            2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2
df.V <- (2 * V^2)/var.V
df.adj <- (df.V + 3)/(df.V + 1)
B.df <- Nchain - 1
W.df <- (2 * w^2)/var.w
R2.fixed <- (Niter - 1)/Niter
R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w)
R2.estimate <- R2.fixed + R2.random
R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) * 
  R2.random
psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
out <- list(psrf = psrf, mpsrf = mpsrf)
print('MPSRF (multivariate bound)')
print(out$mpsrf)
print('PSRFs (univariate bounds)')
print(out$psrf)


### thining based on ESS
extr = as.mcmc(mcmcs[[1]][1:sum(!is.na(mcmcs[[1]][,1])),])
totEffSizes =effectiveSize(extr)
for(ich in 2:length(chs)){
  extr = as.mcmc(mcmcs[[ich]][1:sum(!is.na(mcmcs[[ich]][,1])),])
  totEffSizes = totEffSizes+ effectiveSize(extr)
}
print(data.frame(param=names(totEffSizes),eff_size=as.numeric(totEffSizes)))
setwd(paste0(chDir,'trace plots/'))
save(effS,file = 'effect_size')

######
# Parameters vs iterations
######

burnin = 15000
thin = 450

colos=c(setdiff(allParams,c('pdetec','Beta','popIni','avgAgeRatioIni')))
added = c('M','phi','totIniPop','maxAgeIni','mean_p','nTot',paste0('pdetec_',1:nDS),paste0('Beta_',1:7))
for(ich in 1:length(chs)){
  cat('\n chain ',ich,' \n  \n')
  its = seq(length(chs[[ich]])-thin*((length(chs[[ich]])-burnin)%/%thin),length(chs[[ich]]),by=thin)
  tmp = data.frame(iteration=its,chain=ich)
  for(col in c(colos,added)){
    eval(parse(text=paste('tmp$',col,'=NA',sep="")))}
  for(i in 1:dim(tmp)[1]){
    flush.console()
    cat('\r it',tmp$iteration[i])
    model = get.model(chs[[ich]][[tmp$iteration[i]]],params = allParams)
    tmp$totIniPop[i] = sum(model$popIni)
    tmp$maxAgeIni[i] = max(which(colSums(getIniPop(model$popIni,model$avgAgeRatioIni,maxAge))>0))   
    tmp$M[i] = exp(model$lM)
    tmp$phi[i] = exp(model$lphi)
    for(j in 1:nDS){tmp[i,paste('pdetec_',j,sep="")] = model$pdetec[j]}
    for(j in 1:length(chs[[1]][[2]]$Beta)){
      tmp[i,paste('Beta_',j,sep="")] = model$Beta[j]}
    
    for(col in setdiff(colos,added)){
      eval(parse(text=
                   paste('tmp$',col,'[i]=model$',col,sep="")))
    }
    # Summary statistics
    tmp$mean_p[i] = mean(model$p,na.rm=T)
    tmp$nTot[i] = sum(as.vector(model$n[,,1]))+sum(as.vector(model$iniPop))
  }
  if(ich==1){post=tmp}else{post=rbind(post,tmp)}
}

head(post,digits=2)

post$label= paste('ch',post$chain,'_it',post$iteration,sep="")

# Remove thinned samples
chsNew = list()
for(ich in 1:length(chs)){
  its = seq(length(chs[[ich]])-thin*((length(chs[[ich]])-burnin)%/%thin),length(chs[[ich]]),by=thin)
  tmp=lapply(its,function(it)chs[[ich]][paste0('iteraction_',it)][[1]])
  names(tmp)=names(chs[[ich]])[its]
  chsNew[[ich]] = tmp
  gc(reset=T)
}
chs=chsNew;rm(chsNew);gc(reset=T)

setwd(chDir)
save(post,chs,file = 'toAnalyse.Rdata')

