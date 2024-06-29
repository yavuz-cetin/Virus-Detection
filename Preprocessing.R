if (!require("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") 

BiocManager::install("GEOquery")

library(GEOquery)

if(!file.exists("geo_downloads")) 
  dir.create("geo_downloads")

my.gse <- "GSE71766"  

if(!file.exists(paste0("./geo_downloads/",my.gse)))
  getGEOSuppFiles(my.gse, makeDirectory=T, baseDir="geo_downloads") #makeDirectory T olduğu için subdirectorye ekledi. bu satırın çalışması için options(download.file.method.GEOquery = "curl"). normalde sanırım 'auto'

my.geo.gse <- getGEO(GEO=my.gse, filename=NULL, destdir="./geo_downloads", GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)
my.geo.gse
my.geo.gse <- my.geo.gse[[1]]

untar(paste0("geo_downloads/",my.gse,"/",my.gse,"_RAW.tar"), exdir=paste0("geo_downloads/",my.gse,"/CEL"))

library(stringr)

# Dosyaların bulunduğu dizin
dosya_yolu <- paste0("geo_downloads/", my.gse, "/CEL/")

# Dosyaları listele
dosyalar <- list.files(dosya_yolu, pattern = "\\.gz$")

# Dosya adlarını değiştir
for (dosya in dosyalar) {
  yeni_isim <- str_replace(dosya, "^(GSM\\d+)_.*$", "\\1.CEL.gz")
  file.rename(file.path(dosya_yolu, dosya), file.path(dosya_yolu, yeni_isim))
}

my.cels <- list.files(dosya_yolu)

#Preparing the Phenodata
my.pdata <- as.data.frame(pData(my.geo.gse), stringsAsFactors=F)
head(my.pdata)
dim(my.pdata)
colnames(my.pdata)

table(rownames(my.pdata) == my.cels) #kontrol yapıyoruz, false dönecek

head(my.cels)
head(rownames(my.pdata))

temp.rownames <- paste(rownames(my.pdata), ".CEL.gz", sep="")
table(temp.rownames == my.cels)

rownames(my.pdata) <- temp.rownames
rm(temp.rownames)
rm(dosya)
rm(dosya_yolu)
rm(dosyalar)
rm(yeni_isim)
table(rownames(my.pdata) == my.cels)

write.table(my.pdata, file=paste0("geo_downloads/",my.gse,"/CEL/",my.gse,"_PhenoData.txt"), sep="\t", quote=F)

#Reading the CEL Files
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")

library(affy)
cel.path <- paste0("geo_downloads/",my.gse,"/CEL")
my.affy <- ReadAffy(celfile.path=cel.path, phenoData=paste(cel.path, paste0(my.gse,"_PhenoData.txt"), sep="/"))
show(my.affy)
exprs(my.affy)
head(exprs(my.affy))

colnames(pData(my.affy))
pData(my.affy)
pData(my.affy)$title
pData(my.affy)$description

#Affy nesnesini rma nesnesine dönüştürüp normalize eden fonksiyon:
my.rma <- rma(my.affy, normalize=T, background=T)  #quantile normalizasyonu.
head(exprs(my.rma))

#Annotation#
my.rma@annotation

if (!require("hgu219.db", quietly = TRUE))
  BiocManager::install("hgu219.db")

library(hgu219.db)
library(annotate)
library(R2HTML)

ID <- featureNames(my.rma)
Symbol<-getSYMBOL(ID,"hgu219.db")
sym <- as.data.frame(Symbol)

data <- as.data.frame(exprs(my.rma)) 
data <- cbind(sym,data)  #sütun ekleme

i <- which(is.na(data$Symbol) == TRUE)
data<-data[-c(i),]

rownames(data) <- data[,1] #bu satırda amac satir isimlerini sembol sutunundaki degerlere atamak ancak hata verecek, cunku tekrarli sembol isimleri var.

X <- data.table::as.data.table(data)  #library(data.table), data.table, data.frame'e gore daha fazla islem yapabilmeye olanak saglar.
final_data <- X[,lapply(.SD,mean),"Symbol"] #Symbol sütununda tekrar eden genler var, bu satırların ortalaması alınıp tek bir satır olarak yazılacak.
final_data <- as.data.frame(final_data)
rownames(final_data) <- final_data[,1] 
final_data <- final_data[,-c(1)]

saveRDS(final_data,"geo_downloads/GSE71766/GSE71766_raw.RDS") #okumak için final_data = readRDS("geo_downloads/GSE71766/GSE71766_raw.RDS")

final_data  = t(final_data) 
metadata = pData(my.affy)
table (rownames(final_data) == rownames(metadata)) #her iki datasetteki samplelarin ayni sirada oldugunu kontrol ettik.

final_data = as.data.frame(final_data)
final_data$stage = metadata$virus.infection.ch1
final_data[,"stage"]
final_data$stage = as.numeric(factor(final_data$stage, levels = c("uninfected", "Influenza virus infected", "Rhino virus infected", "Influenza & Rhino virus infected"))) #stage degerlerinin char degil, sayi olmasini istersek bu sekilde degistiririz.
# Uninfected 1  -  Influenza virus infected 2  -  Rhino virus infected 3  -  Influenza & Rhino virus infected 4

write.csv(final_data,file="geo_downloads/GSE71766/GSE71766.csv")