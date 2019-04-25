# Libraries ---------------------------------------------------------------
library(data.table)
library(tidyverse)
library(ggplot2)
library(forcats)

# Import ------------------------------------------------------------------
bd<- "C:\\Users\\User\\Dropbox\\Projetos\\Multiple Comparison Procedure\\Evento - 5WPSM\\Dados"
setwd(bd)


files <- list.files(path = bd,pattern = "errotipoI")
temp <- lapply(files,  function(x)read.csv(x,sep=";"))
tipoI<-rbindlist(temp)
tipoI$dist<-rep(c("Gumbel","Logistic","Normal"),each=84)
tipoI$measure<-rep("Type I",84*3)

files <- list.files(path = bd,pattern = "poder")
temp <- lapply(files, function(x)read.csv(x,sep=";"))
poder<-rbindlist(temp)
poder$dist<-rep(c("Gumbel","Logistic","Normal"),each=1092)
poder$measure<-rep("Power",1092*3)
poder%>%rename(Tamanho=Poder)->poder

dados<-full_join(tipoI,poder)
str(dados)

# Manipulation ------------------------------------------------------------
dados$Teste<-as.factor(dados$Teste)
levels(dados$Teste)[4:7]<-c("Hochberg","Holm","Hommel","None")
dados$Teste<-fct_relevel(dados$Teste,"None","Bonferroni","Holm","Hochberg","Hommel",
                         "BH",'BY')
dados$Teste<-fct_rev(dados$Teste)

#dados%>%select(-Poder)->dados
dados$mu<-as.factor(dados$mu)
dados$sd<-as.factor(dados$sd)
dados$n<-as.factor(dados$n)

levels(dados$sd)<-c("sigma == 1","sigma == 2")

# Type I error ------------------------------------------------------------
head(dados)

dados%>%
  filter(measure=="Type I")%>%
  group_by(dist,sd,n)%>%
  mutate(erro=abs(Tamanho-0.05),
         rank=rank(erro),
         label=paste(round(Tamanho,2),"^(",rank,")",sep = ""))->dados1


  dados1%>%
  ggplot(aes(y=Teste,x=n))+
  facet_grid(sd~dist,labeller = label_parsed)+
  geom_tile(aes(fill=Tamanho),col="black")+
  scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",mid = "white",midpoint = .05,
                       limits=c(min(dados1$Tamanho), max(dados1$Tamanho)))+
  labs(x="Group size",y="Correction",fill="Experimentwise Type I Error Rate")+
  geom_text(aes(label=label),parse=T)+
  theme_bw()+
  theme(legend.position = "top")+
  theme(legend.key.width = unit(2.5, "cm"))+
    theme(text = element_text(size=20))
  
  
  
  dados1%>%filter(dist=="Normal")%>%
    ggplot(aes(y=Teste,x=n))+
    facet_grid(.~sd,labeller = label_parsed)+
    geom_tile(aes(fill=Tamanho),col="black")+
    scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",mid = "white",midpoint = .05,
                         limits=c(min(dados1$Tamanho), max(dados1$Tamanho)))+
    labs(x="Group size",y="Correction",fill="Experimentwise Type I Error Rate")+
    geom_text(aes(label=label),parse=T)+
    theme_bw()+
    theme(legend.position = "top")+
    theme(legend.key.width = unit(2.5, "cm"))+
    theme(text = element_text(size=20))
  
  dados1%>%filter(dist=="Logistic")%>%
    ggplot(aes(y=Teste,x=n))+
    facet_grid(.~sd,labeller = label_parsed)+
    geom_tile(aes(fill=Tamanho),col="black")+
    scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",mid = "white",midpoint = .05,
                         limits=c(min(dados1$Tamanho), max(dados1$Tamanho)))+
    labs(x="Group size",y="Correction",fill="Experimentwise Type I Error Rate")+
    geom_text(aes(label=label),parse=T)+
    theme_bw()+
    theme(legend.position = "top")+
    theme(legend.key.width = unit(2.5, "cm"))+
    theme(text = element_text(size=20))
  

# Power ------------------------------------------------------------
rank(1:3)
  
dados%>%
    filter(measure=="Power")%>%
    mutate(n=paste("n == ",n,sep=""))%>%
    group_by(n,mu,sd,dist)%>%
      mutate(label=paste("(",rank(-Tamanho),")",sep=""))%>%
      filter(mu!=0)->dados2

dados2$n<-as.factor(dados2$n)    
dados2$n<-fct_relevel(dados2$n, "n == 2","n == 3","n == 5","n == 10","n == 15")

dados2$mu<-droplevels(dados2$mu)
      
dados2%>%
  filter(dist=="Gumbel")%>%
    ggplot(aes(y=Teste,x=mu))+
    facet_grid(n~sd,labeller = label_parsed)+
    geom_tile(aes(fill=Tamanho),col="black")+
    scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",mid = "white",midpoint = .5)+
    labs(x="Location Parameter",y="Correction",fill="Empirical Power")+
    geom_text(aes(label=label))+
    theme_bw()+
    theme(legend.position = "top")+
    theme(legend.key.width = unit(2.5, "cm"))+
    theme(text = element_text(size=20))

dados2%>%
  filter(dist=="Logistic")%>%
  ggplot(aes(y=Teste,x=mu))+
  facet_grid(n~sd,labeller = label_parsed)+
  geom_tile(aes(fill=Tamanho),col="black")+
  scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",mid = "white",midpoint = .5)+
  labs(x="Location Parameter",y="Correction",fill="Empirical Power")+
  geom_text(aes(label=label))+
  theme_bw()+
  theme(legend.position = "top")+
  theme(legend.key.width = unit(2.5, "cm"))+
  theme(text = element_text(size=20))  

dados2%>%
  filter(dist=="Normal")%>%
  ggplot(aes(y=Teste,x=mu))+
  facet_grid(n~sd,labeller = label_parsed)+
  geom_tile(aes(fill=Tamanho),col="black")+
  scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",mid = "white",midpoint = .5)+
  labs(x="Location Parameter",y="Correction",fill="Empirical Power")+
  geom_text(aes(label=label))+
  theme_bw()+
  theme(legend.position = "top")+
  theme(legend.key.width = unit(2.5, "cm"))+
  theme(text = element_text(size=20))  

seq<-c(2,3,5,10,15,20)
aux<-paste("n == ",seq,sep="")
lista<-list()



pdf(file = 'Gumbel_power.pdf',width = 9,height = 6)
for( i in 1:6){
dados2%>%
  filter(dist=="Gumbel" & n==aux[i])%>%
  ggplot(aes(y=Teste,x=mu))+
  facet_grid(.~sd,labeller = label_parsed)+
  geom_tile(aes(fill=Tamanho),col="black")+
  scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",
                       mid = "white",midpoint = .5,limits=c(0,1))+
  labs(x="Location Parameter",y="Correction",fill="Empirical Power")+
  geom_text(aes(label=label))+
  theme_bw()+
  theme(legend.position = "top")+
  theme(legend.key.width = unit(2.5, "cm"))+
  theme(text = element_text(size=20))+
  ggtitle(label = paste("n = ",seq[i],sep=""))->p
  print(p)
}
dev.off()


pdf(file = 'Normal_power.pdf',width = 9,height = 6)
for( i in 1:6){
  dados2%>%
    filter(dist=="Normal" & n==aux[i])%>%
    ggplot(aes(y=Teste,x=mu))+
    facet_grid(.~sd,labeller = label_parsed)+
    geom_tile(aes(fill=Tamanho),col="black")+
    scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",
                         mid = "white",midpoint = .5,limits=c(0,1))+
    labs(x="Location Parameter",y="Correction",fill="Empirical Power")+
    geom_text(aes(label=label))+
    theme_bw()+
    theme(legend.position = "top")+
    theme(legend.key.width = unit(2.5, "cm"))+
    theme(text = element_text(size=20))+
    ggtitle(label = paste("n = ",seq[i],sep=""))->p
  print(p)
}
dev.off()

pdf(file = 'Logistic_power.pdf',width = 9,height = 6)
for( i in 1:6){
  dados2%>%
    filter(dist=="Logistic" & n==aux[i])%>%
    ggplot(aes(y=Teste,x=mu))+
    facet_grid(.~sd,labeller = label_parsed)+
    geom_tile(aes(fill=Tamanho),col="black")+
    scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",
                         mid = "white",midpoint = .5,limits=c(0,1))+
    labs(x="Location Parameter",y="Correction",fill="Empirical Power")+
    geom_text(aes(label=label))+
    theme_bw()+
    theme(legend.position = "top")+
    theme(legend.key.width = unit(2.5, "cm"))+
    theme(text = element_text(size=20))+
    ggtitle(label = paste("n = ",seq[i],sep=""))->p
  print(p)
}
dev.off()


###########################################
setwd("C:\\Users\\User\\Dropbox\\UEM\\Projeto-Simulação\\Evento - 5WPSM\\Apresentação\\Figuras")
seq<-c(2,3,5,10,15,20)
aux<-paste("n == ",seq,sep="")
lista<-list()
for( i in 1:6){
  pdf(file = paste('Gumbel_power', i, '.pdf', sep = ''),width = 9,height = 6)
  dados2%>%
    filter(dist=="Gumbel" & n==aux[i])%>%
    ggplot(aes(y=Teste,x=mu))+
    facet_grid(.~sd,labeller = label_parsed)+
    geom_tile(aes(fill=Tamanho),col="black")+
    scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",
                         mid = "white",midpoint = .5,limits=c(0,1))+
    labs(x="Location Parameter",y="Correction",fill="Empirical Power")+
    geom_text(aes(label=label))+
    theme_bw()+
    theme(legend.position = "top")+
    theme(legend.key.width = unit(2.5, "cm"))+
    theme(text = element_text(size=20))+
    ggtitle(label = paste("n = ",seq[i],sep=""))->p
  print(p);dev.off()
}
dev.off()

for( i in 1:6){
  pdf(file = paste('Normal_power', i, '.pdf', sep = ''),width = 9,height = 6)
  dados2%>%
    filter(dist=="Normal" & n==aux[i])%>%
    ggplot(aes(y=Teste,x=mu))+
    facet_grid(.~sd,labeller = label_parsed)+
    geom_tile(aes(fill=Tamanho),col="black")+
    scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",
                         mid = "white",midpoint = .5,limits=c(0,1))+
    labs(x="Location Parameter",y="Correction",fill="Empirical Power")+
    geom_text(aes(label=label))+
    theme_bw()+
    theme(legend.position = "top")+
    theme(legend.key.width = unit(2.5, "cm"))+
    theme(text = element_text(size=20))+
    ggtitle(label = paste("n = ",seq[i],sep=""))->p
  print(p);dev.off()
}
dev.off()



for( i in 1:6){
  pdf(file = paste('Logistic_power', i, '.pdf', sep = ''),width = 9,height = 6)
  dados2%>%
    filter(dist=="Logistic" & n==aux[i])%>%
    ggplot(aes(y=Teste,x=mu))+
    facet_grid(.~sd,labeller = label_parsed)+
    geom_tile(aes(fill=Tamanho),col="black")+
    scale_fill_gradient2(low ="dodgerblue2",high ="firebrick2",
                         mid = "white",midpoint = .5,limits=c(0,1))+
    labs(x="Location Parameter",y="Correction",fill="Empirical Power")+
    geom_text(aes(label=label))+
    theme_bw()+
    theme(legend.position = "top")+
    theme(legend.key.width = unit(2.5, "cm"))+
    theme(text = element_text(size=20))+
    ggtitle(label = paste("n = ",seq[i],sep=""))->p
  print(p);dev.off()
}
dev.off()

























