rm(list = ls())
options(digits=22)
# Bibliotecas -------------------------------------------------------------
pacotes = c("PMCMR", "plyr", "agricolae", "mutoss" ,"doParallel", "foreach")
sapply(pacotes, require, character.only=T)

# Functions-----------------------------------------------------------------
ngroup<-3
nteste<-10

p<-matrix(integer(ngroup*nteste), ncol = nteste) 
valor.p<-function(n,tt){
  grp<-factor(rep(1:ngroup, each = n))
  mod<-aov(tt ~ grp)
  p[,1]<-LSD.test(mod,"grp",p.adj="none",group=F,console = F)$comparison[,2]
  p[,2]<-as.vector(na.exclude(as.vector(pairwise.t.test(tt,grp,p.adjust.method = "bonferroni")$p.value)))
  p[,3]<-TukeyHSD(mod)$grp[,4]
  p[,4]<-SNK.test(mod,"grp",group=F,console = F)$comparison[,2]
  p[,5]<-duncan.test(mod,"grp",group=F,console = F)$comparison[,2]
  p[,6]<-scheffe.test(mod,"grp",group=F,console = F)$comparison[,2]
  # p[,7]<-as.numeric(as.vector(regwq(tt~grp, data=data.frame(tt,grp), alpha=0.05,MSE=NULL, df=NULL, silent = TRUE)$adjPValue))[c(3,1,2)]
  p[,7]<-as.vector(na.exclude(as.vector(posthoc.kruskal.nemenyi.test(x=tt, g=grp, dist="Chisquare")$p.value)))
  p[,8]<-as.vector(na.exclude(as.vector(posthoc.kruskal.dunn.test(x=tt, g=grp, p.adjust.method="bonferroni")$p.value)))
  p[,9]<-as.vector(na.exclude(as.vector(posthoc.kruskal.conover.test(x=tt, g=grp, p.adjust.method="bonferroni")$p.value)))
  p[,10]<-as.vector(na.exclude(as.vector(posthoc.vanWaerden.test(x=tt, g=grp, p.adjust.method="bonferroni")$p.value)))
  return(p)
}

names_test<-c("LSD","t-Bonferroni","Tukey","SNK","Duncan","Scheffé","Nemenyi","Dunn","Conover", "vanWaerden")

B  <- 5000
ni <- c(2, 3, 5, 10, 20)
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
      # exporta para os cluster a variável a ser trabalhada
      # sem esse comando os clusters não reconhecem a variável direto do global env.
      clusterExport(cl, varlist = c("x3", "ngroup", "k", "ni", "B", "nteste", "LSD.test",
                                    "p", "SNK.test", "duncan.test", "scheffe.test","valor.p",
                                    "posthoc.kruskal.nemenyi.test", "posthoc.kruskal.dunn.test",
                                    "posthoc.vanWaerden.test","posthoc.kruskal.conover.test"))      
      system.time(x6[,,k,i,j,] <- array(parSapply(cl, x3, function(u) valor.p(n=ni[k], tt=u), simplify = "array"),c(ngroup, nteste, B)))[3]

      cat(i, j, k, "\n")
    }
  }
}
(fim <- proc.time() - inicio)
stopCluster(cl) # fecha o cluster

# Poder do teste
poder<-apply(apply(x6, 1:6, function(j) sum(j < 0.05))[-1,,,,,], 2:5, function(k) mean(k!=0))
# Erro Tipo I por cexperimento (familywise)
alpha<-apply(apply(x6, 2:6, function(j) sum(j < 0.05)), 1:4, function(k) mean(k!=0))

# Salvando valores-p ------------------------------------------------------
saveRDS(x6, "Simulations_pvalues.rds")

# Salvando poder ----------------------------------------------------------
dados1<-adply(poder,1:4)
colnames(dados1)<-c("Teste","n","sd","mu" ,"Poder")
levels(dados1$Teste)<-names_test
levels(dados1$mu)<-mu.a
levels(dados1$n)<-ni
levels(dados1$sd)<-sd
write.table(dados1,file="poder-TCMM.csv",sep = ";",row.names = F)

# Salvando alpha (tamanho estimado) ---------------------------------------
dados2<-adply(alpha,1:4)
colnames(dados2)<-c("Teste","n","sd","mu" ,"Tamanho")
levels(dados2$Teste)<-names_test
levels(dados2$mu)<-mu.a
levels(dados2$n)<-ni
levels(dados2$sd)<-sd
dados2<-subset(dados2, mu==0)
write.table(dados2,file="errotipoI-TCMM.csv",sep = ";",row.names = F)
