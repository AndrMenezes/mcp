rm(list = ls())
options(digits=22)
# Bibliotecas -------------------------------------------------------------
library(plyr)
library(doParallel)
library(foreach)
library(scmamp)

# Functions-----------------------------------------------------------------
ngroup<-3
nteste<-11

p<-matrix(integer(ngroup*nteste), ncol = nteste) 
valor.p<-function(n,tt){
  grp<-factor(rep(1:ngroup, each = n))
  pv<-as.vector(na.exclude(as.vector(pairwise.t.test(tt,grp,p.adjust.method = "none")$p.value))) 
  p[,1]<-pv
  p[,2]<-p.adjust(pv, method = "bonferroni")
  p[,3]<-p.adjust(pv, method = "holm")
  p[,4]<-p.adjust(pv, method = "hochberg")
  p[,5]<-p.adjust(pv, method = "hommel")
  p[,6]<-p.adjust(pv, method = "BH")
  p[,7]<-p.adjust(pv, method = "BY")
  p[,8]<-adjustFinner(pv)
  p[,9]<-adjustHolland(pv)
  p[,10]<-adjustLi(pv)
  p[,11]<-adjustRom(pv)
  return(p)
}
names_test<-c("None","Bonferroni","Holm","Hochberg","Hommel","BH","BY","Finner", "Holland", "Li", "Rom")

B <- 10000
ni <- c(2,3,5,10,15,20)
mu <- 0
mu.a <- -6:6
sd <- c(1, 2)

x6 <- array(dim = c(ngroup,nteste,length(ni),length(sd),length(mu.a),B))
dimnames(x6) <- list(c("2-1", "3-1", "3-2"),names_test, 
                     paste("n = ", ni),paste("sd = ", sd),
                     paste("mu = ", mu.a), paste("sim",1:B))

# Cria cluster e registra conforme o número of CPU cores.
nodes <- detectCores()
cl <- makeCluster(nodes)
registerDoParallel(cl)

inicio <- proc.time()
set.seed(1502)
foreach(i=1:length(sd)) %do% {
  dados_xy <- rnorm(max(ni)*2*B,mean = mu, sd=sd[i])
  foreach(j=1:length(mu.a)) %do% {
    dados <- c(dados_xy,rnorm(max(ni)*B, mean = mu.a[j], sd=sd[i]))
    mat <- cbind(as.data.frame(matrix(dados, nrow=3*max(ni), ncol=B, byrow=T)),
                 grupo=rep(1:ngroup, each=max(ni)))
    foreach(k=1:length(ni)) %do% {
      x3 <- ddply(mat, .(grupo), function(u) u[1:ni[k],])[,-(B+1)]
      x6[,,k,i,j,] <- array(sapply(x3, function(u) valor.p(n=ni[k], tt=u), 
                                   simplify = "array"),c(ngroup, nteste, B))
      cat(i, j, k, "\n")
    }
  }
}
(fim <- proc.time() - inicio)
stopCluster(cl) # fecha o cluster

# Poder do teste
poder<-apply(apply(x6, 1:6, function(j) sum(j < 0.05))[-1,,,,,], 2:5, function(k) mean(k!=0))
# # Erro Tipo I por experimento (familywise)
alpha<-apply(apply(x6, 2:6, function(j) sum(j < 0.05)), 1:4, function(k) mean(k!=0))

setwd("C:/Users/User/Dropbox/Projetos/Multiple Comparison Procedure/Paper - Correções/Dados")

# Salvando valores-p ------------------------------------------------------
saveRDS(x6, "valorp-normal.rds")

# Salvando poder ----------------------------------------------------------
dados1<-adply(poder,1:4)
colnames(dados1)<-c("Teste","n","sd","mu" ,"Poder")
levels(dados1$Teste)<-names_test
levels(dados1$mu)<-mu.a
levels(dados1$n)<-ni
levels(dados1$sd)<-sd
write.table(dados1,file="poder-normal-10-corrections.csv",sep = ";",row.names = F)

# Salvando alpha (tamanho estimado) ---------------------------------------
dados2<-adply(alpha,1:4)
colnames(dados2)<-c("Teste","n","sd","mu" ,"Tamanho")
levels(dados2$Teste)<-names_test
levels(dados2$mu)<-mu.a
levels(dados2$n)<-ni
levels(dados2$sd)<-sd
dados2<-subset(dados2, mu==0)
write.table(dados2,file="errotipoI-normal-10-corrections.csv",sep = ";",row.names = F)
