homedir <- getwd()
setwd(homedir)

### Load libraries
library(survminer)
library(survival)
library(stringr)
library(xlsx)
dir.create("results")
dir.create("plots")

#### Create Survival object for OV patients (PFS)
df <- read.xlsx2("data/Patients/PFS.xlsx", sheetIndex = 1)
rownames(df) <- df$OC.AXIOMarray.code
df$time <- round(as.numeric(df$PFS..mesi.), digits = 1)
df$status <- as.numeric(df$RECIDIVA.1.yes)
df <- df[,c(4,5)]
df <- na.omit(df)
OSurv <- Surv(df$time,df$status)
names(OSurv) <- rownames(df)

# load table of patients genotypes
tab <- read.delim("results/Genotyping_Table.txt",skip = 5)
colnames(tab) <- sub("_\\..*", "", colnames(tab))
tab[1:3,1:3]
rownames(tab) <- tab[,1]
tab <- tab[,-1]
common <- intersect(names(OSurv), colnames(tab))
OSurv <- OSurv[common]
tab <- tab[, common]

### add filtered SNPs list (coming from Array Analysis with the Axiom Analysis Suite)
ref <- read.delim("results/Genotyping_table_filtered.txt")
rownames(ref) <- ref$probeset_id
#https://erc.bioscientifica.com/view/journals/erc/23/3/171.xml
bh <- 0.05/nrow(ref)
### Apply a log-rank test
# https://github.com/kassambara/survminer/issues/97

source("D:/Archive/stepsurvival_cox.R")
# Create a Vector with Column names
columns = c("Affy-ID","Gene","SNP rs ID", "p-val") 

#Create a Empty DataFrame with 0 rows and n columns
tabella = data.frame(matrix(nrow = 0, ncol = length(columns))) 

# Assign column names
#colnames(tabella) = columns

# print DataFrame
print(tabella)
values <- rep(NA,4)
tabella<-rbind(tabella, values)
colnames(tabella) = columns

for (mysnp in rownames(ref)){
  message("Doing: ",mysnp)
  oritrack <- tab[mysnp,]
  track <- oritrack
  track <- as.factor(track)
  if (nlevels(track)<=1) {
    print("less than 1")
  } else if (nlevels(track)==4){
    print("4")
    comm <- intersect(names(OSurv),names(track))
    survival <- OSurv[comm]
    genus <- track[comm]
    survival <- cbind(survival)
    datafr <- data.frame(survival,genus)
    colnames(datafr)[3]<-"genotype"
    datafr <- as.data.frame(datafr)
    datafr$time <- as.numeric(datafr$time)
    datafr$status <- as.numeric(datafr$status)
    datafr <- datafr[datafr$genotype!="NoCall",]
    # comm <- intersect(rownames(survival),rownames(datafr))
    # survival <- survival[comm]
    # datafr <- datafr[comm,]
    cox <- survdiff(Surv(datafr$time,datafr$status) ~ genotype, data = datafr)
    l <- cox$obs
    m <- (l[1]+0.1)/(l[2]+0.1)
    p <- 1-pchisq(cox$chisq, df=length(cox$n) - 1)
    message("p-value: ",p)
    #p <- p*3
    events_difference <- cox$obs-cox$exp
      if (p <= bh & min(table(track))>=5){
        km_fit <- survfit(Surv(datafr$time,datafr$status) ~ genotype, data = datafr)
        rs <- ref[mysnp,c(29:30)]
        gene <- strsplit(rs$Associated.Gene,"//")
        gene <- unlist(gene)[c(5,12)]
        gene <- gsub(" ","",gene)
        ifelse(gene[1]!="---",gene <- gene[1],gene <- paste0("dn-",gene[2]))
        TIT <- paste(mysnp,",", gene[1],",",rs$dbSNP.RS.ID, "adj-p = ",signif(p,3))
        gp <- ggsurvplot(km_fit,data = datafr, title = TIT,
                         risk.table = T)+xlab("Months")
        png(paste0("plots/000_",mysnp,"_PFS.png"),w=2500,h=2500,res=300)
        print(gp)
        dev.off()
        values<-c(mysnp,gene,rs$dbSNP.RS.ID,signif(p,3))
        tabella <- rbind(tabella, values)
      }
  } else {
    print("good")
    comm <- intersect(names(OSurv),names(track))
    survival <- OSurv[comm]
    genus <- track[comm]
    survival <- cbind(survival)
    datafr <- data.frame(survival,genus)
    colnames(datafr)[3]<-"genotype"
    datafr <- as.data.frame(datafr)
    datafr$time <- as.numeric(datafr$time)
    datafr$status <- as.numeric(datafr$status)
    cox <- survdiff(Surv(datafr$time,datafr$status) ~ genotype, data = datafr) #pairwise_survdiff would be better...
    # datafr <- cbind(OSurv,genotype)
    # datafr <-as.data.frame(datafr)
    # cox <- pairwise_survdiff(Surv(time, status) ~ genotype, data = datafr)
    l <- cox$obs
    m <- (l[1]+0.1)/(l[2]+0.1)
    p <- 1-pchisq(cox$chisq, df=length(cox$n) - 1)
    message("p-value: ",p)
    #p <- p*3
    events_difference <- cox$obs-cox$exp
      if (p <= bh & min(table(track))>=5){
        km_fit <- survfit(Surv(datafr$time,datafr$status) ~ genotype, data = datafr)
        rs <- ref[mysnp,c(29:30)]
        gene <- strsplit(rs$Associated.Gene,"//")
        gene <- unlist(gene)[c(5,12)]
        gene <- gsub(" ","",gene)
        ifelse(gene[1]!="---",gene <- gene[1],gene <- paste0("dn-",gene[2]))
        TIT <- paste(mysnp,",", gene,",", rs$dbSNP.RS.ID, "adj-p = ",signif(p,3))
        gp <- ggsurvplot(km_fit,data = datafr, title = TIT,
                         risk.table = T)+xlab("Months")
        png(paste0("plots/000_",mysnp,"_PFS.png"),w=2500,h=2500,res=300)
        print(gp)
        dev.off()
        values<-c(mysnp,gene,rs$dbSNP.RS.ID,signif(p,3))
        tabella <- rbind(tabella, values)
      }
  }
}

save(tabella,file="results/tabella_PFS.rda")
write.xlsx(tabella, file = "results/tabella_PFS.xlsx",row.names = F)
nrow(tabella[tabella$Gene!="---",])
