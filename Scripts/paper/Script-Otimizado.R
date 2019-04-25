rm(list = ls())
options(digits=22)
# Bibliotecas -------------------------------------------------------------
library(PMCMR)
library(plyr)
library(agricolae)

# Functions-----------------------------------------------------------------
ngroup<-3
nteste<-10
# sd.y<-sd.x<-1  
# sd.z<-1


p<-matrix(integer(ngroup*nteste), ncol = nteste) 

valor.p <- function(n, mu.x, mu.z, sd){

mu.y<-mu.x 
sd.x <- sd.y <- sd.z <- sd

tt<-c(rnorm(n, mu.x, sd.x),rnorm(n, mu.y, sd.y),rnorm(n, mu.z, sd.z))
grp<-factor(rep(1:ngroup, each = n))

mod<-aov(tt ~ grp)
  
p[,1]<-TukeyHSD(mod)$grp[,4]
p[,2]<-as.vector(na.exclude(as.vector(pairwise.t.test(tt,grp,p.adjust.method = "none")$p.value)))
p[,3]<-as.vector(na.exclude(as.vector(posthoc.kruskal.nemenyi.test(x=tt, g=grp, dist="Chisquare")$p.value)))
p[,4]<-as.vector(na.exclude(as.vector(posthoc.kruskal.dunn.test(x=tt, g=grp, p.adjust.method="bonferroni")$p.value)))
p[,5]<-as.vector(na.exclude(as.vector(posthoc.kruskal.conover.test(x=tt, g=grp, p.adjust.method="bonferroni")$p.value)))
p[,6]<-as.vector(na.exclude(as.vector(pairwise.t.test(tt,grp,p.adjust.method = "bonferroni")$p.value)))
# p[,7]<-as.vector(na.exclude(as.vector(pairwise.t.test(tt,grp,p.adjust.method = "holm")$p.value)))
# p[,8]<-as.vector(na.exclude(as.vector(pairwise.t.test(tt,grp,p.adjust.method = "hochberg")$p.value)))
# p[,9]<-as.vector(na.exclude(as.vector(pairwise.t.test(tt,grp,p.adjust.method = "hommel")$p.value)))
# p[,10]<-as.vector(na.exclude(as.vector(pairwise.t.test(tt,grp,p.adjust.method = "BH")$p.value)))
# p[,11]<-as.vector(na.exclude(as.vector(pairwise.t.test(tt,grp,p.adjust.method = "BY")$p.value)))
# p[,12]<-as.vector(na.exclude(as.vector(pairwise.t.test(tt,grp,p.adjust.method = "fdr")$p.value)))
p[,7]<-LSD.test(mod,"grp",p.adj="bon",group=F,console = F)$comparison[,2]
p[,8]<-SNK.test(mod,"grp",group=F,console = F)$comparison[,2]
p[,9]<-duncan.test(mod,"grp",group=F,console = F)$comparison[,2]
p[,10]<-HSD.test(mod,"grp",group=F,console = F)$comparison[,2]
return(p)
}

# names_test<-c("Tukey","t","Nemenyi","Dunn","Conover","t-Bonferroni","t-Holm","t-Hochberg","t-Hommel",
#               "t-BH","t-BY","t-fdr","LSD","SNK","Duncan",'HSD')
names_test<-c("Tukey","t","Nemenyi","Dunn","Conover","t-Bonferroni",
              "LSD","SNK","Duncan",'HSD')

cmapply <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, 
                    USE.NAMES = TRUE)
{ return(do.call(mapply, c(
    list(FUN=FUN, MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES), 
    expand.grid(..., stringsAsFactors=FALSE)
  )))
}


#alfa <- function(arrray){
#  apply(apply(arrray, c(2,3,4,5), function(j) sum(j < 0.05)), c(1,2,3), function(k) mean(k!=0))
#}

#poder <- function(arrray){
#  apply(apply(array, c(1,2,3,4,5), function(j) sum(j < 0.05))[-1,,,,], c(2,3,4), function(k) mean(k!=0))
#}


# Poder do teste ----------------------------------------------------------
B <- 500
ni = c(2,3,5,10,20)
mu <- 0
mu.a <- 0
sd.y <- c(1,2,3)

inicio <- proc.time()
x2<-array(replicate(B, cmapply(valor.p, n = ni, mu.x = mu, 
                               mu.z = mu.a, sd = sd.y,
                               SIMPLIFY = "array")),
          c(ngroup,nteste,length(ni),length(mu.a),length(sd.y),B))
fim <- proc.time() - inicio

dimnames(x2) <- list(c("2-1", "3-1", "3-2"),names_test, 
                    paste("n = ", ni),paste("mu = ", mu.a),
                    paste("sd = ", sd.y),
                    paste("sim",1:B))

alpha <- apply(apply(x2, 2:6, function(j) sum(j <= 0.05)),1:4 , function(k) mean(k!=0))


dados1 <- adply(alpha, 1:4)[,-3]
colnames(dados1)<-c("Teste","n","sd","Tamanho")
levels(dados1$Teste)<-names_test
levels(dados1$n)<-ni
levels(dados1$sd)<-sd.y
head(dados1)

options(digits = 4)
dados1 %>% arrange(Teste) %>% filter(sd==1) %>%
  ggplot(aes(n,Tamanho,group=Teste,col=Teste))+
  geom_line()+geom_point()+
  scale_colour_brewer(palette = "Paired")+
  scale_y_continuous(breaks = seq(0, 0.25, by=0.05))+
  geom_hline(yintercept = 0.05,linetype="dashed")+
  theme(axis.title.x = element_blank()) -> p1

dados1 %>% arrange(Teste) %>% filter(sd==2) %>%
  ggplot(aes(n,Tamanho,group=Teste,col=Teste))+
  geom_line()+geom_point()+
  scale_colour_brewer(palette = "Paired")+
  scale_y_continuous(breaks = seq(0, 0.25, by=0.05))+
  geom_hline(yintercept = 0.05,linetype="dashed")+
  theme(axis.title.x = element_blank()) -> p2

dados1 %>% arrange(Teste) %>% filter(sd==3) %>%
  ggplot(aes(n,Tamanho,group=Teste,col=Teste))+
  geom_line()+geom_point()+
  scale_colour_brewer(palette = "Paired")+
  scale_y_continuous(breaks = seq(0, 0.25, by=0.05))+
  geom_hline(yintercept = 0.05,linetype="dashed")+
  theme(axis.title.x = element_blank()) -> p3

dados1 %>% arrange(Teste) %>%
  ggplot(aes(n,Tamanho,group=Teste,col=Teste))+
  facet_grid(.~sd,scales="free")+
  geom_line()+geom_point()+
  scale_colour_brewer(palette = "Paired")+
  scale_y_continuous(breaks = seq(0, 0.25, by=0.05))+
  geom_hline(yintercept = 0.05,linetype="dashed")+
  theme(axis.title.x = element_blank()) -> p4

x2
#Salva valores p
saveRDS(x2, "Simulations_pvalues.rds")

poder<-apply(apply(x2, c(1,2,3,4,5), function(j) sum(j < 0.05))[-1,,,,], c(2,3,4), function(k) mean(k!=0))

dados<-adply(poder,1:3)
colnames(dados)<-c("Teste","n","mu" ,"Poder")
levels(dados$Teste)<-names_test
levels(dados$mu)<-mu.a

#Salva o poder
write.table(dados,file="Simulations_poder.csv",sep = ";",row.names = F)



# Plots -----------------------------------------------------------------
#HeatMap
library(dplyr)
library(ggplot2)

cor<-c('green','yellow','red')

p1<-dados%>%arrange(Teste)%>%ggplot( aes(mu, Teste, fill = Poder))+ facet_grid(.~n,scales="free")
p1<- p1+  geom_tile(color = "white")+ labs(list(x="Média",y=""))
p1<- p1+   scale_fill_gradient2(low = cor[1], high = cor[3], mid = cor[2], 
                                midpoint = .5, limit = c(0,1), space = "Lab",name="Poder")  
p1<- p1+ theme_minimal()+   coord_fixed()
#p1<- p1+ geom_text(aes(x=mu, y=Teste,fill = Poder, label = round(Poder, 2)))
p1<- p1+ theme(legend.position="left",legend.key.size = unit(1.5, "cm"))
#p1<- p1+ theme(axis.text.x=element_text())
p1<- p1+ scale_x_discrete(breaks = seq(-8,8,by=2),labels=seq(-8,8,by=2))
p1
#Alpha
dados%>%filter(mu==0)%>%ggplot(aes(Teste,Poder))+geom_point()+facet_grid(~n)+
  theme(axis.text.x=element_text(angle=90)) +geom_hline(yintercept = 0.05,linetype="dashed")+
  theme(axis.title.x = element_blank())

dados%>%filter(mu==0)%>%ggplot(aes(n,Poder,fill=Teste))+geom_point()+
  theme(axis.text.x=element_text(angle=90)) +geom_hline(yintercept = 0.05,linetype="dashed")+
  theme(axis.title.x = element_blank())
