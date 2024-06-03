# Prateća skripta za udžbenik za TOBS predmet pod nazivom
# "Analiza biosignala sa praktičnim primerima u programskom jeziku R"
# autorke Nadice Miljković

################################# POGLAVLJE 1. #################################
## 1. Osnovni pojmovi i korisne informacije za rad u programskom jeziku R

# funkcija za instalaciju paketa
install.packages("dplyr")

# funkcija za aktivaciju paketa (neophodna za korišćenje funkcija i podataka 
# iz paketa)
library(dplyr)

# računarska reproducibilnost za dplyr paket
if (!require("groundhog")) install.packages("groundhog")
library(groundhog)
groundhog.day = "2023-07-30"
pkgs = c('dplyr')

# za računarsku reproducibilnost veceg broja paketa
# potrebno je promeniti samo prethodnu liniju koda:
# 
groundhog.library(pkgs, groundhog.day, tolerate.R.version = '4.1.2')

# deljenje sa nulom je moguće u R-u
p <- 4
p / 0

## 1.1 Tipovi podataka u R-u

0 / p # operacija je validna 
0 / 0 # rezultat je NaN
TRUE # logicka promenljiva
T # skraćen zapis logicke promenljive
FALSE
F

# definisanje nizova i implicitna tj. automatska konverzija
niz1 <- c(TRUE, FALSE, 15, 13, "c", "tobs")
niz1
class(niz1)

niz2 <- c(TRUE, FALSE, 15, 13)
niz2
class(niz2)

niz3 <- c(TRUE, FALSE, "c", "tobs")
niz3
class(niz3)

# primer eksplicitne konverzije
as.logical(niz3)

# boje u R-u
# nazivi boja su preuzeti iz kataloga:
# http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf, pristupljeno 17.05.2024.
c("aquamarine2", "red", "blue", "chocolate",
  "goldenrod", "deeppink")

# kreiranje listi i tabelarnih podataka
dat <- list()
dat$niz2 <- niz2
dat$niz3 <- niz3
dat

dat <- data.frame(dat)
dat

dat2 <- data.frame()
dat2$niz2 <- niz2
dat2

class(dat)
str(dat) # str funkcija se koristi za proveru strukture objekta

# kreiranje matrica 
m1 <- matrix(37:44, ncol = 4, nrow = 2)
m1

dimnames(m1)
dimnames(m1) <- list(c("vrsta1", "vrsta"),
                     c("kol1", "kol2", "kol3", "kol4"))

m1
attributes(m1)

# kreiranje faktorske promenljive
# ilustracija za radne dane u sedmici
f1 <- as.factor(c("pon", "ut", "sr", "cet", "pet"))
f1
f1 <- factor(c("pon", "ut", "sr", "cet", "pet"),
             levels = c("pon", "ut", "sr", "cet", "pet"),
             ordered = TRUE)
f1

# korišćenje operatora za indeksiranje "$"
dat$niz2
naziv_kolone = "niz2"
dat$naziv_kolone

# logičko indeksiranje (podaci iz ISwR paketa)
if (!require("ISwR")) install.packages("ISwR")
library(ISwR)
pod <- stroke
head(pod, 2)
class(pod$sex)
sum(pod$sex == "Female")
length( pod$sex[pod$sex == "Female"] )

# broj dana od poslednje stabilne verzije R-a do danasnjeg dana
Sys.time()
datum1 <- as.Date("2023-06-16")
datum2 <- as.Date("2023-07-30")
datum2 - datum1

# primer upotrebe funkcije za pseudoslučajno odabiranje
# obratiti paznju da ponovljeno pokretanje sample() funkcije
# daje različite odnosno pseudoslučajne rezultate
d <- 1:20
sample(d, 5)
sample(d[d < 5], 5, replace = TRUE)

# primer upotrebe funkcija za rad sa NA vrednostima
x <- airquality[1:6,]
x
is.na(x$Solar.R)
sum(is.na(x$Solar.R))
x[complete.cases(x),]

# 1.1.2 Uvoz podataka u R
?read.table
?read.csv
?write.table

# primeri korišćenja url() funkcije
# primer iz vežbe 3 TOBS predmeta
# sajtu je pristupljeno 17.05.2024.
con <- url("https://automatika.etf.bg.ac.rs/sr/13m051tobs")
x <- readLines(con)
head(x, 10)
close(con)

# primer iz vežbe 5 TOBS predmeta
dat <- read.table(url("http://www.statsci.org/data/oz/ms212.txt"),
                  header = TRUE)
head(dat)

# učitavanje podataka i postavljanje radnog direktorijum
dat <- read.table("C:\Users\BioEng\Downloads\ms212.txt",
                  header = TRUE) # neispravno
dat <- read.table("C:/Users/BioEng/Downloads/ms212.txt",
                  header = TRUE) # ispravno
head(dat, 2)
setwd("C:/Users/BioEng/Downloads")
getwd()
if (!require("rstudioapi")) install.packages("rstudioapi")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# merenje signala u R-u primenom serial paketa
# deo koda je preuzet/promenjen sa:
# https://stackoverflow.com/questions/34522535/how-to-read-data-from-serial-port-in-r
# (pristupljeno 17.05.2024.)

library(serial)
con <- serialConnection(name = "test_con",
                        port = "COM3",
                        mode = "9600,n,8,1",
                        newline = 1,
                        translation = "cr")

# otvaranje konekcije (uvek je prati close() komanda)
open(con)

n <- 1
foo <- ""
textSize <- 0

stopT <- Sys.time() + 7
foo <- ""
textSize <- 0
while(Sys.time() < stopT){
  # čitanje tekstualnog podatka sa serijskog porta
  newText <- read.serialConnection(con)
  if(0 < nchar(newText))
  {
    # smeštanje u promenljivu foo
    foo <- paste(foo, newText)
  }
}

# čuvanje podataka u datoteku
cat("\r\n", foo, "\r\n", file = "podaci.txt")

# zatvaranje konekcije
close(con)

# učitavanje sačuvanih podataka i prikaz na grafiku
dat <- read.table("podaci.txt")
plot(dat$V1, type = "l")

# uputstvo za plot funkciju
?plot

# 1.2 Operatori u R-u i kontrolne strukture
# vektorske operacije nad nizovima
x <- 1:4; y <- 5:8
x
y
x + y
x > 2
x > 100
x >= 3
y == 6
x * y
y / x

# vektorske operacije nad nizovima različitih dužina
x <- 1:4; y <- 1:3
x - y
x * y

# dodela vrednosti
x <- 1:50
x
x <- 1
x
x <- "gold"
print(x)
x <- 'gold'
x

# množenje matrica
x <- matrix(1:4, 2, 2)
y <- matrix(rep(10, 4), 2, 2)
x
y
x * y
x %*% y

# kontrolne strukture
# uslovna if/else struktura
x <- 3
if (x > 3) {
  y <- "veceOd3"
} else {
  y <- "manjeIliJednako3"
}
y

x <- 2
y <- if (x > 3) {
  "veceOd3"
} else {
  "manjeIliJednako3"
}
y

x <- 7
y <- if (x > 3) "veceOd3" else "manjeIliJednako3"
y

# for struktura
for (ind in 1:5) {
  print(ind)
}

for (ix in 1:3) print(ix)

x <- c("t", "o", "b", "s")
for (ix in seq_along(x)) print(x[ix])
for (letter in x) print(letter)

# kompleksni brojevi i indeksi for petlje
sqrt(-1)
sqrt(as.complex(-1))
sqrt(-1 + 0i)
2 + 3i
2 + 3j
i**2
1i^2

# primer primene while strukture
iterator <- 0
while (iterator < 5) {
  iterator <- iterator + 1
  print(iterator)
}
iterator

while (iterator < 5 && iterator >=0) {
  iterator <- iterator + 1
  print(iterator)
}
iterator

# kontrolne funkcije
# primer primene lapply() funkcije
x <- list(prvi = 1:10, drugi = c(NA, 2, 5, 9, 10))
lapply(x, mean)
lapply(x, mean, na.rm <- TRUE)
lapply(x, mean, na.rm = TRUE)

# definisanje anonimne funkcije u okviru kontrolne lapply() funkcije
x <- list(prvi = matrix(1:10, 2, 5), drugi = matrix(1:10, 5, 2))
x
lapply(x, function(pom) pom[,1])

# primena gl(), kao i split() funkcije u kombinaciji sa sapply() funkcijom
gl(3,5)
library(help = "datasets")
head(PlantGrowth, 3)
unique(PlantGrowth$group)
dat <- split(PlantGrowth$weight, PlantGrowth$group)
class(dat)
rez <- sapply(dat, mean, na.rm = TRUE)
rez

# primer primene mapply() funkcije
rep(2)
rep(2, 4)
rep(c(3, 5, 4), 12)
rep(1:3, 3)
rep(3, 1:3)
mapply(rep(3, 1:3))
mapply(rep, 3, 1:3)

# 1.3 Funkcije u R-u
class(mean)
class(sum)

# primeri funkcija
f1 <- function() {}
f2 <- function() {
  print("Primer funkcije.")
}

f1()
f2()
f3()

f3 <- function(ar) {
  print("Primer funkcije sa jednim argumentom.")
  print(ar)
}
f3()
f3(ar = 37)

# primer funkcije sa podrazumevanim argumentom
f4 <- function(ar = 5) {
  print("Primer funkcije sa podrazumevanim argumentom.")
  print(ar)
}
f4()
formals(f4)

# argument "..."
?rnorm
noviRnorm <- function(mean = 2.5, ...) {
  rnorm(mean = 2.5, ...)
}
x <- noviRnorm(300, sd = 1, mean = 2.5)
mean(x)
length(x)

# definisanje funkcije koja ima isto ime kao postojeća
mean
mean <- function(x) x**2
mean
mean(3)
search()
mean <- "tobs"
mean

# provera koji su objekti definisani u okruženju
ls()

# leksička pravila pretrage argumenata
f6 <- function(x, y) {
  x + y + z
}
f6(2, 3)
z <- 1
f6(2, 3)
f7 <- function(x, y) {
  x + y + z
  invisible()
}
f7(2, 3)

# dinamička (funkcija1(10) = 250) 
# vs. leksička pretraga (funkcija1(10) = 100)
x <- 2
funkcija1 <- function(arg1) {
  x <- 5
  x * funkcija2(arg1)
}
funkcija2 <- function(arg1) {
  arg1 * x
}
funkcija1(10)

################################# POGLAVLJE 2. #################################
# 2.2.1.1 Primena dplyr paketa za pripremu podataka za dalju obradu
dat <- read.csv("PrimerPodatakaProsireno.csv")
head(dat, 3)
dat <- read.csv("PrimerPodatakaProsireno.csv", sep = ";")
head(dat, 3)

class(dat$pol)
class(dat$zdrav)
dat$pol <- as.factor(dat$pol)
dat$zdrav <- as.factor(dat$zdrav)
class(dat$pol)
class(dat$zdrav)

# filter() funkcija iz dplyr paketa
library(dplyr)
dat1 <- filter(dat, dat$skocniZglob < 65)
class(dat1)
dat1
dat2 <- dplyr::filter(dat, skocniZglob < 65)
dat2 == dat1
identical(dat1, dat2)

# select() funkcija iz dplyr paketa
head(dat, 3)
tail(select(dat, koleno:zdrav), 3)
tail(select(dat, -(koleno:zdrav)), 3)
tail(select(dat, skocniZglob, koleno), 3)

# funkcija arrange() iz dplyr paketa
head(dat, 3)
head(arrange(dat, koleno), 3)
tail(arrange(dat, koleno), 3)
head(arrange(dat, desc(koleno)), 3)

# funkcije rename() i mutate() iz dplyr paketa
head(dat, 3)
head(rename(dat, ugaoKoleno = koleno), 3)
head(mutate(dat, koleno = koleno - mean(koleno)), 3)
head(mutate(dat, koleno = koleno*pi/180), 3)

# linearni model za transformaciju jedinica iz stepena u V
dat$skocniZglobV <- dat$skocniZglob/90 + 2.5
head(dat$skocniZglob, 3)
head(dat$skocniZglobV, 3)

plot(dat$skocniZglob, dat$skocniZglobV)
lines(dat$skocniZglob, dat$skocniZglobV, col = "plum3")

A <- matrix(c(dat$skocniZglob[1], dat$skocniZglob[2], 
              1, 1), nrow = 2, ncol = 2)
p <- matrix(c(dat$skocniZglobV[1], dat$skocniZglobV[2]),
            nrow = 2, ncol = 1)
A
solve(A)
solve(A) %*% p # ili solve(A, p)

# računanje srednje vrednosti za ROM u skočnom zglobu za
# zdrave ispitanike i pacijente
dat$skocniZglob[c(1, 3, 4, 5, 7)] <- NA
mean(dat$skocniZglob[dat$zdrav == "da"])
mean(dat$skocniZglob[dat$zdrav == "da"], na.rm = T)
mean(dat$skocniZglob[dat$zdrav == "ne"], na.rm = T)

sum(is.na(dat$skocniZglob[dat$zdrav == "da"]))
sum(dat$zdrav == "da")
sum(is.na(dat$skocniZglob[dat$zdrav == "ne"]))
sum(dat$zdrav == "ne")

# primena funkcija group_by() i summarise() iz dplyr paketa
datNovo <- group_by(dat, zdrav)
class(datNovo)
head(datNovo, 3)
summarize(datNovo, skocniZglob = median(skocniZglob))
summarize(datNovo, skocniZglob = median(skocniZglob, na.rm = T))

# korišćenje pipeline operatora za prethodni primer
dat %>% group_by(zdrav) %>% summarize(median(skocniZglob, na.rm = T)) 

# priprema podataka za dalju analizu
# provera osnovnih elemenata u podacima
dat <- read.csv("PrimerPodatakaProsireno.csv", sep = ";")
sum(dat$zdrav == "da")
sum(dat$zdrav == "ne")
sum(dat$skocniZglob < 0)
sum(dat$skocniZglob < 65)
sum(dat$koleno < 0)
sum(dat$koleno < 130)
dat$koleno
sum(dat$koleno < 130*0.9)

# manipulacija nedostajućim vrednostima
# najpre se nedostajuće vrednosti postavljaju, jer ih nema u podacima
# ovo je urađeno radi ilustracije kako manipulacija nedostajućim 
# vrednostima utiče na krajnji rezultat
dat <- read.csv("PrimerPodatakaProsireno.csv", sep = ";")
median(dat$koleno)
median(dat$koleno[dat$zdrav == "da"])
median(dat$koleno[dat$zdrav == "ne"])
dat$koleno[c(1, 2, 3, 4, 19, 20, 21, 22)] <- NA
dat$koleno[c(1:4, 19:22)] <- NA
kol <- median(dat$koleno, na.rm = TRUE)
kol
dat$koleno[c(1:4, 19:22)] <- kol
median(dat$koleno)
median(dat$koleno[dat$zdrav == "da"])
median(dat$koleno[dat$zdrav == "ne"])
length(dat$skocniZglob)

# bolje rešenje za zamenu nedostajućih vrednosti
dat <- read.csv("PrimerPodatakaProsireno.csv", sep = ";")
dat$koleno[c(1:4, 19:22)] <- NA
kolZ <- dat %>% filter(zdrav == "da") %>% select(koleno) %>%
  unlist() %>% median(na.rm = TRUE)
kolZ
kolP <- dat %>% filter(zdrav == "ne") %>% select(koleno) %>%
  unlist() %>% median(na.rm = TRUE)
kolP
dat$koleno[is.na(dat$koleno) == TRUE & dat$zdrav == "da"] <- kolZ
dat$koleno[is.na(dat$koleno) == TRUE & dat$zdrav == "ne"] <- kolP
median(dat$koleno)
median(dat$koleno[dat$zdrav == "da"])
median(dat$koleno[dat$zdrav == "ne"])

# 2.2.2 Eliminacija šuma na primeru EMG biosignala
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
datEMG <- read.table("emg-sve.txt")
datEMGkabl <- read.table("emg-kabl.txt")
fs <- 1000
timeAxis <- seq(0, length(datEMG$V1)/fs - 1/fs, by = 1/fs)
timeAxisKabl <- seq(0, length(datEMGkabl$V1)/fs - 1/fs, by = 1/fs)

jpeg(file = "emg.jpg", width = 6, height = 7, 
     units = 'in', res = 300)
par(mar = c(4.2, 4.2, 1.7, 1.7)) 
par(mfrow = c(3, 1))
plot(timeAxis, datEMG$V1, type = "l", xlim = c(64, 70),
     ylim = c(-1.7, 1.7),
     xlab = "vreme [s]", ylab = "napon [mV]",
     main = "sirov EMG signal")
grid()
plot(timeAxis, datEMG$V1, type = "l", xlim = c(154, 158),
     ylim = c(-5.2, 5.2),
     xlab = "vreme [s]", ylab = "napon [mV]",
     main = "EMG signal sa odlepljivanjem elektrode")
grid()
plot(timeAxisKabl, datEMGkabl$V1, type = "l", xlim = c(40, 50),
     ylim = c(-0.45, 0.45),
     xlab = "vreme [s]", ylab = "napon [mV]",
     main = "EMG signal sa pomeranjem kablova")
grid()
dev.off()

# Furijeova transformacija
library(signal)        

# FFT signala bez šuma - emgm
podFFT <- fft(datEMG$V1[(64*fs):(70*fs)])
podMag <- Mod(podFFT)        
podMag <- podMag[1:length(podMag)/2]
fosa <- 1:length(podMag)/
  (length(datEMG$V1[(64*fs):(70*fs)])/fs)

jpeg(file = "emg-fft.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
  plot(fosa, podMag, type = "l", xlim = c(0, 500),
       xlab = "frekvencija [Hz]",
       ylab = "magnituda [a.u.]")
    grid()
dev.off()

# FFT signala sa šumom
podFFT <- fft(datEMGkabl$V1[(40*fs):(50*fs)])
podMag <- Mod(podFFT)        
podMag <- podMag[1:length(podMag)/2]
fosa <- 1:length(podMag)/
  (length(datEMGkabl$V1[(40*fs):(50*fs)])/fs)

jpeg(file = "emg-fft-kabl.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
plot(fosa, podMag, type = "l", xlim = c(0, 500),
     xlab = "frekvencija [Hz]",
     ylab = "magnituda [a.u.]")
grid()
dev.off()

# FFT signala sa šumom/odlepljene elektrode
podFFT <- fft(datEMG$V1[(154*fs):(164*fs)])
podMag <- Mod(podFFT)        
podMag <- podMag[1:length(podMag)/2]
fosa <- 1:length(podMag)/
  (length(datEMG$V1[(154*fs):(164*fs)])/fs)

jpeg(file = "emg-fft-50.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
plot(fosa, podMag, type = "l", xlim = c(0, 500),
     xlab = "frekvencija [Hz]",
     ylab = "magnituda [a.u.]")
grid()
dev.off()

# Primena rev() funkcije
?rev
rev(c(1:7))

# Batervortov filtar reda 2, 4. i 8. reda
library(signal)
fs <- 1000 # frekvencija odabiranja
filt2 <- butter(2, 10/(fs/2), "high")
filt4 <- butter(4, 10/(fs/2), "high")
filt8 <- butter(8, 10/(fs/2), "high")
f2 <- freqz(filt2, Fs = fs)
f4 <- freqz(filt4, Fs = fs)
f8 <- freqz(filt8, Fs = fs)

jpeg(file = "Bater-tri-reda.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 0.1, 0.1)) 
plot(f2$f, abs(f2$h), type="l", lwd=3, lty=1, # kontinualna linija (eng. solid)
     xlab="frekvencija [Hz]", ylab="magnituda [a.u.]",
     ylim=c(0, 1), xlim=c(0, 40)) 
par(new = T)
plot(f4$f, abs(f4$h), type="l", lwd=3, lty=2, # isprekidana linija (eng. dashed)
     xlab="frekvencija [Hz]", ylab="magnituda [a.u.]", 
     ylim=c(0, 1), xlim=c(0, 40)) 
par(new = T)  
plot(f8$f, abs(f8$h), type="l", lwd=3, lty=3, # tačkasta linija (eng. dotted)
     xlab="frekvencija [Hz]", ylab="magnituda [a.u.]", 
     ylim=c(0, 1), xlim=c(0, 40)) 
grid()
legend(0, 1, legend=c('2. red','4. red','8. red'),  
       lty=1:3, lwd=3)
dev.off() 

# Odskočni odziv Batervort filtara na jediničnu funkciju
n = 2^8
jpeg(file = "Bater-tri-reda-step.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 1, 1)) 
par(mfrow = c(1, 2))
plot(filter(filt2, c(0.5, rep(1, n-1))), type="l", lwd=3, lty=1, # kontinualna linija (eng. solid)
     xlab="odbirci", ylab="amplituda [a.u.]",
     ylim = c(-0.5, 1), xlim = c(0, n),
     main = "filter")
par(new = T)
plot(filter(filt4, c(0.5, rep(1,n-1))), type="l", lwd=3, lty=2, # isprekidana linija (eng. dashed)
     xlab="odbirci", ylab="amplituda [a.u.]",
     ylim = c(-0.5, 1), xlim = c(0, n))
par(new = T)  
plot(filter(filt8, c(0.5, rep(1, n-1))), type="l", lwd=3, lty=3, # tačkasta linija (eng. dotted)
     xlab="odbirci", ylab="amplituda [a.u.]", 
     ylim=c(-0.5, 1), xlim=c(0, n)) 
grid()
# filtfilt() funkcija
plot(filtfilt(filt2, c(0.5, rep(1, n-1))), type="l", lwd=3, lty=1, # kontinualna linija (eng. solid)
     xlab="", ylab="",
     ylim = c(-0.5, 1), xlim = c(0, n),
     main = "filtfilt")
par(new = T)
plot(filtfilt(filt4, c(0.5, rep(1,n-1))), type="l", lwd=3, lty=2, # isprekidana linija (eng. dashed)
     xlab="", ylab="",
     ylim = c(-0.5, 1), xlim = c(0, n))
par(new = T)  
plot(filtfilt(filt8, c(0.5, rep(1, n-1))), type="l", lwd=3, lty=3, # tačkasta linija (eng. dotted)
     xlab="", ylab="", 
     ylim=c(-0.5, 1), xlim=c(0, n)) 
grid()
legend(20, 1, legend=c('2. red','4. red','8. red'),  
       lty=1:3, lwd=3)
dev.off() 

# Odskočni odziv na izmenjenu realizaciju Hevisajdove funkcije
jpeg(file = "Bater-tri-reda-step-2.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 1, 1)) 
par(mfrow = c(1, 2))
plot(filter(filt2, rep(1, n)), type="l", lwd=3, lty=1, # kontinualna linija (eng. solid)
     xlab="odbirci", ylab="amplituda [a.u.]",
     ylim = c(-0.5, 1), xlim = c(0, n),
     main = "filter")
par(new = T)
plot(filter(filt4, rep(1,n)), type="l", lwd=3, lty=2, # isprekidana linija (eng. dashed)
     xlab="odbirci", ylab="amplituda [a.u.]",
     ylim = c(-0.5, 1), xlim = c(0, n))
par(new = T)  
plot(filter(filt8, rep(1, n)), type="l", lwd=3, lty=3, # tačkasta linija (eng. dotted)
     xlab="odbirci", ylab="amplituda [a.u.]", 
     ylim=c(-0.5, 1), xlim=c(0, n)) 
grid()
# filtfilt() funkcija
plot(filtfilt(filt2, rep(1, n)), type="l", lwd=3, lty=1, # kontinualna linija (eng. solid)
     xlab="", ylab="",
     ylim = c(-0.5, 1), xlim = c(0, n),
     main = "filtfilt")
par(new = T)
plot(filtfilt(filt4, rep(1,n)), type="l", lwd=3, lty=2, # isprekidana linija (eng. dashed)
     xlab="", ylab="",
     ylim = c(-0.5, 1), xlim = c(0, n))
par(new = T)  
plot(filtfilt(filt8, rep(1, n)), type="l", lwd=3, lty=3, # tačkasta linija (eng. dotted)
     xlab="", ylab="", 
     ylim=c(-0.5, 1), xlim=c(0, n)) 
grid()
legend(50, 1, legend=c('2. red','4. red','8. red'),  
       lty=1:3, lwd=3)
dev.off() 

# Filtriranje EMG signala
# filtriranje šuma pokreta i šuma napajanja
library(signal)
datEMGkablAll <- read.table("emg-sve.txt")
fs <- 1000
datEMGkabl <- datEMGkablAll$V1[(4*fs):(14*fs)]
timeAxisKabl <- seq(0, length(datEMGkabl)/fs - 1/fs, by = 1/fs)

fb <- butter(2, 10/(fs/2), "high")
f50 <- butter(3, c(49/(fs/2), 51/(fs/2)), "stop")
datEMGkablFilt <- filtfilt(f50, filtfilt(fb, datEMGkabl))

emgFFT <- fft(datEMGkabl); emgMag <- Mod(emgFFT)
emgMag <- emgMag[1:(length(emgMag)/2)]
emgFiltFFT <- fft(datEMGkablFilt); emgFiltMag <- Mod(emgFiltFFT)        
emgFiltMag <- emgFiltMag[1:(length(emgFiltMag)/2)]
fosa <- 1:length(emgMag) / (length(datEMGkabl)/fs)

jpeg(file = "EMG-dva-filtra.jpg", width = 8, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 0.1, 0.1)) 
par(mfrow = c(2, 2))
plot(timeAxisKabl, datEMGkabl, type = "l", 
     ylab = "amplituda [mV]", xlab = "vreme [s]")
grid()
plot(timeAxisKabl, datEMGkablFilt, type = "l",
     ylab = "", xlab = "vreme [s]")
grid()
plot(fosa, emgMag, type = "l",
     ylab = "magnituda [a.u.]", xlab = "frekvencija [Hz]")
grid()
plot(fosa, emgFiltMag, type = "l",
     ylab = "", xlab = "frekvencija [Hz]")
grid()
dev.off()

# Računanje obvojnice EMG signala
fs <- 1000
datEMGkablAll <- read.table("emg-sve.txt")
datEMGkabl <- datEMGkablAll$V1[(4*fs):(14*fs)]
timeAxisKabl <- seq(0, length(datEMGkabl)/fs - 1/fs, by = 1/fs)

datEMGkabl <- abs(datEMGkabl) # ispravljanje signala
ma50 <- rep(1/50, 50)
datEMGkablMA <- stats::filter(datEMGkabl, ma50, sides = 1)

jpeg(file = "ma-red-50.jpg", width = 8, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 0.1, 0.1)) 
plot(timeAxisKabl, datEMGkabl, type = "l",
     ylab = "amplituda [mV]", xlab = "vreme [s]",
     col = "cyan", ylim = c(0, 1.3))
grid()
par(new = TRUE)
plot(timeAxisKabl, datEMGkablMA, type = "l",
     ylab = "amplituda [mV]", xlab = "vreme [s]",
     ylim = c(0, 1.3))
legend(0, 1.2, legend = c("ispravljen EMG", "obvojnica"),  
       fill = c("cyan", "black"))
dev.off()

# Drugi način realizacije MA filtra
library(signal)
a50 <- c(1, rep(0, 50 - 1))
b50 <- rep(1/50, 50)

datEMGkablMA2 <- signal::filter(b50, a50, datEMGkabl)

jpeg(file = "ma-red-50-2.jpg", width = 8, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 0.1, 0.1)) 
plot(timeAxisKabl, datEMGkabl, type = "l",
     ylab = "amplituda [mV]", xlab = "vreme [s]",
     col = "cyan", ylim = c(0, 1.3))
grid()
par(new = TRUE)
plot(timeAxisKabl, datEMGkablMA2, type = "l",
     ylab = "amplituda [mV]", xlab = "vreme [s]",
     ylim = c(0, 1.3))
legend(0, 1.2, legend = c("ispravljen EMG", "obvojnica"),  
       fill = c("cyan", "black"))
dev.off()

# frekvencijska karakteristika MA filtra
a50 <- c(1, rep(0, 50 - 1))
b50 <- rep(1/50, 50)
f50 <- freqz(b50, a50, Fs = fs)

a100 <- c(1, rep(0, 100 - 1))
b100 <- rep(1/100, 100)
f100 <- freqz(b100, a100, Fs = fs)

a200 <- c(1, rep(0, 200 - 1))
b200 <- rep(1/200, 200)
f200 <- freqz(b200, a200, Fs = fs)

jpeg(file = "MA-tri-reda.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 0.1, 0.1)) 
plot(f50$f, abs(f50$h), type="l", lwd=3, lty=1,
     xlab="frekvencija [Hz]", ylab="magnituda [a.u.]",
     ylim=c(0, 1), xlim=c(0, 40)) 
par(new = T)
plot(f100$f, abs(f100$h), type="l", lwd=3, lty=2,
     xlab="frekvencija [Hz]", ylab="magnituda [a.u.]", 
     ylim=c(0, 1), xlim=c(0, 40)) 
par(new = T)  
plot(f200$f, abs(f200$h), type="l", lwd=3, lty=3,
     xlab="frekvencija [Hz]", ylab="magnituda [a.u.]", 
     ylim=c(0, 1), xlim=c(0, 40)) 
grid()
legend(25, 1, legend=c('50 odbiraka','100 odbiraka','200 odbiraka'),  
       lty=1:3, lwd=3)
dev.off() 

# sinc funkcija
x <- seq(-50, 50, by = 0.01)
sinc <- sin(x)/x; sinc[sinc == Inf] <- 1;

jpeg(file = "sinc.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 0.1, 0.1)) 
plot(x, sinc, type="l", lwd=3, lty=1,
     ylab="sinc(x)", xlab="x") 
grid()
dev.off()

# Hanov prozor i faza linearnih filtara u vremenskom domenu
library(signal)
ma5 <- rep(1/5, 5)
hn5 <- hanning(5)
fs <- 1000

fma <- freqz(ma5, Fs = fs)
fhn <- freqz(hn5, Fs = fs)

jpeg(file = "faza-frekvencija-ma-hn.jpg", width = 6.5, height = 5, 
     units = 'in', res = 300)
par(mar = c(4.2, 4.2, 0.3, 0.3)) 
par(mfrow = c(2, 1))
plot(fma$f, abs(fma$h), type="l", lwd=3, lty=1,
     xlab="", ylab="magnituda [a.u.]",
     ylim=c(0, 2), xlim=c(0, 500)) 
grid()
par(new = T)
plot(fhn$f, abs(fhn$h), type="l", lwd=3, lty=2,
     xlab="", ylab="magnituda [a.u.]", 
     ylim=c(0, 2), xlim=c(0, 500)) 
legend(300, 2, legend=c('MA filtar','Hanov filtar'),  
       lty=1:2, lwd=3)
plot(fma$f, atan2(Im(fma$h), Re(fma$h)), type="l", lwd=3, lty=1,
     xlab="frekvencija [Hz]", ylab="faza [rad]", 
     ylim = c(-3, 3), xlim=c(0, 500)) 
grid()
par(new = T)
plot(fhn$f, atan2(Im(fhn$h), Re(fhn$h)), type="l", lwd=3, lty=2,
     xlab="frekvencija [Hz]", ylab="faza [rad]", 
     ylim = c(-3, 3), xlim=c(0, 500)) 
grid()
dev.off() 

# prikaz faze bez diskontinuiteta (korišćenje unwrap() funkcije)
jpeg(file = "faza-ma-hn-unwrap.jpg", width = 6.5, height = 3.5, 
     units = 'in', res = 300)
par(mar = c(4.2, 4.2, 0.3, 0.3)) 
plot(fma$f, unwrap(atan2(Im(fma$h), Re(fma$h))), type="l", lwd=3, lty=1,
     xlab="frekvencija [Hz]", ylab="faza [rad]", 
     ylim = c(-12, 12), xlim=c(0, 500)) 
par(new = T)
plot(fhn$f, unwrap(atan2(Im(fhn$h), Re(fhn$h))), type="l", lwd=3, lty=2,
     xlab="frekvencija [Hz]", ylab="faza [rad]", 
     ylim = c(-12, 12), xlim=c(0, 500)) 
grid()
legend(0, 10, legend=c('MA filtar','Hanov filtar'),  
       lty=1:2, lwd=3)
dev.off() 

# Poređenje tri MA filtra različite dužine prozora za filtriranje EMG signala
datEMGkablMA50 <- signal::filter(b50, a50, datEMGkabl)
datEMGkablMA100 <- signal::filter(b100, a100, datEMGkabl)
datEMGkablMA200 <- signal::filter(b200, a200, datEMGkabl)

jpeg(file = "ma-tri-reda-signal.jpg", width = 8, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 0.1, 0.1)) 
plot(timeAxisKabl, datEMGkabl, type = "l", lwd = 3,
     ylab = "amplituda [mV]", xlab = "vreme [s]",
     col = "cyan", ylim = c(0, 1.3), xlim = c(6.7, 7.5))
grid()
par(new = TRUE)
plot(timeAxisKabl, datEMGkablMA50, type = "l",  lty = 1, lwd = 2,
     ylab = "amplituda [mV]", xlab = "vreme [s]",
     ylim = c(0, 1.3), xlim = c(6.7, 7.5))
par(new = TRUE)
plot(timeAxisKabl, datEMGkablMA100, type = "l",  lty = 2, lwd = 2,
     ylab = "amplituda [mV]", xlab = "vreme [s]",
     ylim = c(0, 1.3), xlim = c(6.7, 7.5))
par(new = TRUE)
plot(timeAxisKabl, datEMGkablMA200, type = "l",  lty = 3, lwd = 2,
     ylab = "amplituda [mV]", xlab = "vreme [s]",
     ylim = c(0, 1.3), xlim = c(6.7, 7.5))
legend(6.8, 1.2, legend = c("ispravljen EMG", "50 ms", "100 ms", "200 ms"),  
       col = c("cyan", "black", "black", "black"),
       lty = c(1, 1, 2, 3))
dev.off()

# Frakcioni izvod u R-u
fs <- 50
t <- seq(0, 10, by = 1/fs)
f <- 1 + sin(0.5*t*2*pi) + 0.4*sin(2*pi*t) - 0.2*sin(4*t*2*pi)
alpha = 1.5 # red frakcionog izvoda

# preuzeti fod() funkciju prema uputstvu iz udžbenika

fod02 <- fod(t, f, fs, 0.2)
fod04 <- fod(t, f, fs, 0.4)
fod06 <- fod(t, f, fs, 0.6)
fod08 <- fod(t, f, fs, 0.8)
fod10 <- fod(t, f, fs, 1)
fod12 <- fod(t, f, fs, 1.2)

jpeg(file = "foc.jpg", width = 8, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 0.1, 0.1)) 
plot(t, f/max(f), type = "l", lwd = 3,
     ylab = "amplituda [a.u.]", xlab = "vreme [s]",
     col = "black", ylim = c(-0.5, 1.2))
grid()
par(new = TRUE)
plot(t, fod02/max(fod02, na.rm = T), type = "l", 
     ylab = "amplituda [a.u.]", xlab = "vreme [s]",
     col = "red", ylim = c(-0.5, 1.2))  
par(new = TRUE)
plot(t, fod04/max(fod04, na.rm = T), type = "l", 
     ylab = "amplituda [a.u.]", xlab = "vreme [s]",
     col = "cyan", ylim = c(-0.5, 1.2))  
par(new = TRUE)
plot(t, fod06/max(fod06, na.rm = T), type = "l", 
     ylab = "amplituda [a.u.]", xlab = "vreme [s]",
     col = "orange", ylim = c(-0.5, 1.2))  
par(new = TRUE)
plot(t, fod08/max(fod08, na.rm = T), type = "l", 
     ylab = "amplituda [a.u.]", xlab = "vreme [s]",
     col = "blue", ylim = c(-0.5, 1.2))  
par(new = TRUE)
plot(t, fod10/max(fod10, na.rm = T), type = "l", 
     ylab = "amplituda [a.u.]", xlab = "vreme [s]",
     col = "grey", ylim = c(-0.5, 1.2))  
par(new = TRUE)
plot(t, fod12/max(fod12, na.rm = T), type = "l", 
     ylab = "amplituda [a.u.]", xlab = "vreme [s]",
     col = "darkred", ylim = c(-0.5, 1.2))  
legend(7.5, 1.2, legend = c("originalni signal", "izvod 0.2. reda",
                            "izvod 0.4. reda", "izvod 0.6. reda",
                            "izvod 0.8. reda",
                            "izvod 1. reda", "izvod 1.2. reda"),  
       fill = c("black", "red", "cyan", "orange", "blue",
                "grey", "darkred"))
dev.off()

# izoelektrična linija
dat <- read.table("ecg_70.txt")
ecg7 <- dat$id7
ecg9 <- dat$id9

library(signal)
fs <- 2000
f1 <- butter(3, 1/(fs/2), "high")
ecg7 <- filtfilt(f1, ecg7)
ecg9 <- filtfilt(f1, ecg9)
t <- seq(0, length(ecg7)/fs - 1/fs, by = 1/fs)

# specijalizovana biblioteka za crtanje objekata na grafiku
library(plotrix)

jpeg(file = "isoelectric.jpg", width = 8, height = 4, 
     units = 'in', res = 300)
par(mar = c(4, 4, 0.1, 0.1), mfrow = c(1,2)) 
plot(t, ecg7, type = "l", lwd = 2,
     ylab = "amplituda [mV]", xlab = "vreme [s]",
     col = "black", xlim = c(10, 15))
draw.circle(x = 10.7, y = 0, radius=.3, border = "cyan", lwd = 2)
abline(h = 0, lty = 2, col = "black", lwd = 2)
grid()
plot(t, ecg9, type = "l", lwd = 2,
     ylab = "", xlab = "vreme [s]",
     col = "black", xlim = c(10, 15))
draw.circle(x = 10.5, y = 0, radius=.25, border = "cyan", lwd = 2)
abline(h = 0, lty = 2, col = "black", lwd = 2)
grid()
dev.off()

################################# POGLAVLJE 3. #################################
# 3. Simulirani biosignali
?rnorm
?set.seed

rnorm(5)
rnorm(5)
set.seed(38); rnorm(5)
set.seed(38); rnorm(5)
rnorm(5)
set.seed(45); rnorm(5)
set.seed(45); rnorm(5)

# primena summary() funkcije
y <- rnorm(10); y
summary(y)
options(digits = 2)
summary(y)

jpeg(file = "odbirci.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
plot(rnorm(100), ylim = c(-4, 4))
polygon(c(-4, 104, 104, -4), c(1, 1, -1, -1),
        col = rgb(0, 1, 1, alpha = 0.4))
grid()
text(50, 3.5, "oblast unutar standardne devijacije (68.2%)")
dev.off()

jpeg(file = "fgv.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
plot(seq(-4, 4, by = 0.1), 
     dnorm(seq(-4, 4, by = 0.1), mean = 0, sd = 1), type = "l",
     lwd = 2)
grid()
polygon(c(-1, 1, 1, -1), c(0.42, 0.42, -0.1, -0.1),
        col = rgb(0, 1, 1, alpha = 0.4))
text(0, 0.2, "68.2%")
dev.off()

# verovatnoća da je broj jednak 0.5 i koji broj ima verovatnoću od 0.69
pnorm(mean = 0, sd = 1, 0.5)
qnorm(0.69, mean = 0, sd = 1)
# verovatnoća da se pseudoslučajni broj nalazi u opsegu
pnorm(mean = 0, sd = 1, 1)
pnorm(mean = 0, sd = 1, 4)
pnorm(mean = 0, sd = 1, 400)

# Generisanje belog šuma
# signal je dostupan: 
# https://automatika.etf.bg.ac.rs/images/FAJLOVI_srpski/predmeti/izborni_kursevi_os/biomedicinsko_inzenjerstvo/TOBS/vezbe/puls.txt
# (pristupljeno 17.05.2024.)
dat <- read.table("puls.txt")
fs <- 1000

ecg <- dat$V1[1:(10*fs)]
set.seed(14)
bs <- 0.11 * rnorm(10*fs)

ecg_sum <- ecg + bs
vreme <- seq(0, length(ecg)/fs - 1/fs, by = 1/fs)

# prikaz EKG signala sa šumom u vremenskom domenu
jpeg(file = "ecg-sum.jpg", width = 8, height = 5, 
     units = 'in', res = 300)

plot(vreme, ecg_sum, type = "l", xlab = "vreme [s]", 
     ylab = "amplituda [mV]", xlim = c(0, 10), ylim = c(-0.7, 0.7),
     col = "cyan")
par(new = TRUE)
plot(vreme, ecg, type = "l", xlab = "vreme [s]", 
     ylab = "amplituda [mV]", xlim = c(0, 10), ylim = c(-0.7, 0.7),
     col = "black")
grid()
legend(6, 0.7, legend = c("EKG sa belim šumom", "EKG"),  
       fill = c("cyan", "black"))
dev.off()

# prikaz EKG signala sa šumom u frekvencijskom domenu
ecgFFT <- Mod(fft(ecg))
ecgFFT <- ecgFFT[1:length(ecgFFT)/2]
ecg_sumFFT <- Mod(fft(ecg_sum))
ecg_sumFFT <- ecg_sumFFT[1:length(ecg_sumFFT)/2]

fosa <- 1:length(ecgFFT) / (length(ecg)/fs)

jpeg(file = "ecg-sum-fft.jpg", width = 8, height = 5, 
     units = 'in', res = 300)

plot(fosa, ecg_sumFFT, type = "l", xlab = "frekvencija [Hz]", 
     ylab = "magnituda [a.u.]", xlim = c(0, 500), ylim = c(0, 150),
     col = "cyan")
par(new = TRUE)
plot(fosa, ecgFFT, type = "l", xlab = "frekvencija [Hz]", 
     ylab = "magnituda [a.u.]", xlim = c(0, 500), ylim = c(0, 150),
     col = "black")
grid()
legend(300, 150, legend = c("EKG sa belim šumom", "EKG"),  
       fill = c("cyan", "black"))

dev.off()

# generisanje pink / roze šuma
library(tuneR)
pinksum <- noise(kind = c("pink"))
pinksum
fs <- pinksum@samp.rate
vreme <- seq(0, length(pinksum@left)/fs - 1/fs, by = 1/fs)
pinksumFFT <- Mod(fft(pinksum@left))
pinksumFFT <- pinksumFFT[1:length(pinksumFFT)/2]
fosa <- 1:length(pinksumFFT) / (length(pinksum@left)/fs)

# prikaz pink šuma sa Furijeovom transformacijom na grafiku
jpeg(file = "pink-sum.jpg", width = 8, height = 5, 
     units = 'in', res = 300)
par(mar = c(4.2, 4.2, 1, 1)) 
par(mfrow = c(2,1))
plot(vreme, pinksum@left, type = "l", col = "pink",
     ylab = "amplituda [a.u.]", xlab = "vreme [s]")
grid()
plot(fosa, pinksumFFT, type = "l", xlab = "frekvencija [Hz]", 
     ylab = "magnituda [a.u.]", col = "pink")
grid()
dev.off()

# prikaz raspodele EMG signala
library(ggplot2)
library(gridExtra)
sig <- read.csv("EMG.csv")
fs <- 1000
sig$vreme <- seq(0, (length(sig$EMG)/fs - 1/fs), 1/fs)

b <- round(sqrt(length(sig$EMG)) + 1)
bw <- (max(sig$EMG) - min(sig$EMG)) / b
emgHist <- ggplot(sig, aes(x=EMG)) +
  geom_histogram(aes(y = ..density..), binwidth = bw, 
                 color = "#000000", fill = "white") + 
  geom_density(alpha = 0.2, color = "cyan", fill = "cyan") +
  xlab("napon [mV]") + ylab("br. merenja")
emgTime <- ggplot(sig, aes(x=vreme, y=EMG)) +
  geom_line(color = "black")+ xlab("vreme [s]") +
  ylab("napon [mV]")
g <- arrangeGrob(emgHist, emgTime)
ggsave("EMGhist.jpg", g, dpi = 300, units = "in",
       width = 6.5, height = 3.58)

# može se učitati i signal iz EMG.txt, kako ima puno šuma
# dobiće se Laplasova transformacija

# korišćenje ugrađene R funkcije za prikaz histograma
jpeg(file = "histogrami.jpg", width = 8, height = 6, 
     units = 'in', res = 300)
par(mfrow = c(2, 2))
hist(sig$EMG, 2, xlab = "napon [mV]", ylab = "frekvencija",
     main = paste(as.character(2), "intervala", sep = " "), col = "cyan")
hist(sig$EMG, xlab = "napon [mV]", ylab = "frekvencija",
     main = "Podrazumevan broj intervala", col = "cyan")
hist(sig$EMG, 40, xlab = "napon [mV]", ylab = "frekvencija",
     main = paste(as.character(40), "intervala", sep = " "), col = "cyan")
hist(sig$EMG, b, xlab = "napon [mV]", ylab = "frekvencija",
     main = paste(as.character(b), "intervala", sep = " "), col = "cyan")
dev.off()
dev.off()

# QQ grafik za procenu raspodele odbiraka EMG signala
jpeg(file = "QQceoEMG.jpg", width = 6, height = 4.5, 
     units = 'in', res = 300)
qqnorm(sig$EMG, col = "cyan", main = "QQ grafik",
       xlab = "teorijski kvantili", ylab = "kvantili odbiraka")
qqline(sig$EMG)
grid()
dev.off()

# QQ grafici za dva segmenta
jpeg(file = "QQsegmentiEMG.jpg", width = 8, height = 6, 
     units = 'in', res = 300)
par(mfrow = c(2, 1))
qqnorm(sig$EMG[(5*fs):(6*fs)], col = "cyan", main = "QQ grafik za kontrakciju",
       xlab = "teorijski kvantili", ylab = "kvantili odbiraka")
qqline(sig$EMG[(5*fs):(6*fs)], lwd = 2, lty = 2)
grid()
qqnorm(sig$EMG[1:fs], col = "cyan", main = "QQ grafik za relaksaciju",
       xlab = "teorijski kvantili", ylab = "kvantili odbiraka")
qqline(sig$EMG[1:fs], lwd = 2, lty = 2)
grid()
dev.off()
dev.off()

# poređenje histograma tokom kontrakcije i relaksacije primenom ggplot2 paketa
emg <- list()
emg$napon <- c(sig$EMG[1:fs], sig$EMG[(5*fs):(6*fs-1)])
emg$vreme <- c(sig$vreme[1:fs], sig$vreme[(5*fs):(6*fs-1)])
emg$stanje <- c(rep("relaksacija", fs), rep("kontrakcija", fs))
emg <- data.frame(emg)

ggplot(emg, aes(x = napon, color = stanje)) +
  geom_histogram(fill = "white") + 
  scale_color_manual(values = c("cyan", "black")) +
  xlab("br. merenja") + ylab("napon [mV]")

ggsave("poredjenje_histograma.jpg", dpi = 300)

# simulacija EMG signala
library(biosignalEMG)

jpeg(file = "simuliranEMG.jpg", width = 8, height = 6, 
     units = 'in', res = 300)
semg <- syntheticemg()
plot(semg$values, xlab = "vreme [ms]", ylab = "napon [mV]",
     type = "l", col = "cyan")
grid()
dev.off()

# prikaz trenutaka u kojima je mišić aktivan i u kojima nije aktivan
library(ggplot2)

emgDf <- list()
emgDf$values <- semg$values
emgDf$state <- semg$on.off 
emgDf$vreme <- (1:length(emgDf$values))/1000 # u sekundama
emgDf <- data.frame(emgDf)

ggplot(emgDf, aes(x = vreme, y = values)) + 
  geom_line(aes( group = 1, color= (state < 1) )) +
  xlab("vreme [s]") + ylab("napon [mV]") +
  scale_color_manual(values = c("cyan", "black"))
ggsave("sintetickiEMG.jpg", dpi = 300, units = "in",
       width = 6.5, height = 3.58)

# prikaz obvojnice primenom biosignalEMG paketa
emgr <- rectification(semg, rtype = "fullwave")
emgma <- envelope(semg, method = "MA", wsize = 60)

jpeg(file = "emg-obvojnica-ispravljanje.jpg", width = 6.5, height = 4, 
     units = 'in', res = 300)
plot(emgDf$vreme, emgr$values, main = "Sintetički EMG i obvojnica",
     xlab = "vreme [s]", ylab = "amplituda [mV]", 
     col = "cyan", type = "l")
lines(emgDf$vreme, emgma$values, main = "Sintetički EMG i obvojnica",
      xlab = "vreme [s]", ylab = "amplituda [mV]",
      col = "black", lwd = 2, lty = 2)
grid()
legend("topright", c("ispravljen EMG", "obvojnica"),
       col = c("cyan", "black"), lty = c(1, 2), lwd = c(1, 2), cex = 0.8)
dev.off()

# faza simuliranog EMG signala
library(signal)
library(biosignalEMG)
semg <- syntheticemg()
emgDf <- list()
emgDf$values <- semg$values
emgDf$state <- semg$on.off 
emgDf$vreme <- (1:length(emgDf$values))/1000 # u sekundama
emgDf <- data.frame(emgDf)

emgFFT <- fft(emgDf$values)
emgMag <- Mod(emgFFT)        
emgMag <- emgMag[1:(length(emgMag)/2)]
emgFaza <- unwrap(atan2(Im(emgFFT), Re(emgFFT)))
emgFaza <- emgFaza[1:(length(emgFaza)/2)]
fs <- 1000
fosa <- 1:length(emgMag) / (length(emgDf$values)/fs)

jpeg(file = "emg-fft-simulirani.jpg", width = 8, height = 6, 
     units = 'in', res = 300)
par(mfrow = c(3, 1))
par(mar = c(4, 4, 0.1, 0.1)) 
plot(emgDf$vreme, emgDf$values, type = "l",
     xlab = "vreme [s]", ylab = "napon [mV]")
grid()
plot(fosa, emgFaza, type = "l", col = "cyan",
     xlab = "frekvencija [Hz]", ylab = "faza [rad]")
grid()
plot(fosa, emgMag, type = "l",
     xlab = "frekvencija [Hz]", ylab = "magnituda [a.u.]")
grid()
dev.off()

jpeg(file = "histogram-faza-fft-simulirani.jpg", width = 6, height = 4.5, 
     units = 'in', res = 300)
hist(emgFaza, main = "", xlab = "faza EMG signala [rad]",
     ylab = "br. odbiraka", col = "cyan")
dev.off()

# faza realnog EMG signala
library(signal)
emgDf <- read.csv("EMG.csv")
emgDf$values <- emgDf$EMG
emgDf$vreme <- (1:length(emgDf$values))/1000 # u sekundama
emgDf <- data.frame(emgDf)

emgFFT <- fft(emgDf$values)
emgMag <- Mod(emgFFT)        
emgMag <- emgMag[1:(length(emgMag)/2)]
emgFaza <- unwrap(atan2(Im(emgFFT), Re(emgFFT)))
emgFaza <- emgFaza[1:(length(emgFaza)/2)]
fs <- 1000
fosa <- 1:length(emgMag) / (length(emgDf$values)/fs)

jpeg(file = "emg-fft-realni.jpg", width = 8, height = 6, 
     units = 'in', res = 300)
par(mfrow = c(3, 1))
# da se smanji prazan beli prostor koristi se mar argument od eng. margine
par(mar = c(4, 4, 0.1, 0.1))
plot(emgDf$vreme, emgDf$values, type = "l",
     xlab = "vreme [s]", ylab = "napon [mV]")
grid()
plot(fosa, emgFaza, type = "l", col = "cyan",
     xlab = "frekvencija [Hz]", ylab = "faza [rad]")
grid()
plot(fosa, emgMag, type = "l",
     xlab = "frekvencija [Hz]", ylab = "magnituda [a.u.]")
grid()
dev.off()
dev.off()

jpeg(file = "histogram-faza-fft-realni.jpg", width = 6, height = 4.5, 
     units = 'in', res = 300)
hist(emgFaza, main = "", xlab = "faza EMG signala [a.u.]",
     ylab = "br. odbiraka", col = "cyan")
dev.off()

# Hanov prozor
n <- 4;
hann <- 0.5 * (1 - cos(2*pi*(0:n)/n))
hann
hanning(5)

# inverzna Furijeova transformacija
emgMag <- c(emgMag, emgMag)
emgFaza <- c(emgFaza, emgFaza)

emgNovFFT = emgMag * exp(1i*emgFaza);
emgNov = Re(ifft(emgNovFFT));

emgNovFFT2 = emgMag * exp(1i*runif(length(emgFaza),
                                   min = min(emgFaza), max = max(emgFaza)));
emgNov2 = Re(ifft(emgNovFFT2));

emgNovFFT3 = runif(length(emgMag), min = min(emgMag), max = max(emgMag)) *
  exp(1i*emgFaza);
emgNov3 = Re(ifft(emgNovFFT3));

jpeg(file = "emg-ifft.jpg", width = 8, height = 6, 
     units = 'in', res = 300)
par(mfrow = c(3, 1))
par(mar = c(4.3, 4.3, 1.4, 1.4))
plot(emgDf$vreme, emgNov, type = "l",
     main = "originalan EMG signal",
     xlab = "vreme [s]", ylab = "napon [mV]")
grid()
plot(emgDf$vreme, emgNov2, type = "l",
     main = "EMG signal sa uniformnom fazom",
     xlab = "vreme [s]", ylab = "amplituda [a.u.]")
grid()
plot(emgDf$vreme, emgNov3, type = "l",
     main = "EMG signal sa uniformnom magnitudom",
     xlab = "vreme [s]", ylab = "amplituda [a.u.]")
grid()
dev.off()
dev.off()

# Transformacije kao predkorak u izdvajanju obeležja
# Sledeći deo koda je preuzet i modifikovan sa GitHub-a:
# https://github.com/NadicaSm/satRday-Belgrade-2018/
# koji je podeljen pod otvorenom GNU GPL licencom.
library(EMD)

# suma prostoperiodičnih signala x(t)
t <- seq(0, 10, length = 1000)
x <- sin(pi*t) + sin(2*pi*t) + sin(5*pi*t)  - 0.7*t

# primena EMD tehnike
rez <- emd(x, t, boundary = "wave")

jpeg('suma-sin-za-emd.jpg', units='in', width=6, height=4, res=400)
plot(t, x, type = "l", xlab = "t", 
     ylab = expression(sin(pi*t) + sin (2 * pi * t) +
                         sin(5 * pi * t) - 0.7*t))
grid()
dev.off()

jpeg('rezultatEMD.jpg', units='in', width=6, height=5, res=400)
par(mfrow = c(rez$nimf+1, 1), mar = c(2,1,2,1))
rangeimf <- range(rez$imf)
for(ind in 1:rez$nimf) {
  plot(t, rez$imf[, ind], type="l", 
       xlab = "", ylab= "", ylim = rangeimf,
       main = paste(ind, ". IMF", sep = ""))
  abline(h = 0, lty = 2)
  grid()
}
plot(t, rez$residue, xlab = "t",
     ylab = "", main = "rezidual", type = "l")
abline(h = 0, lty = 2)
grid()
dev.off()

# EMD primenjena na EMG signalima merenih sa Pectoralis major (lat.) mišića
dat <- read.table("EMGpectoralis.txt")
emg <- dat$V1[1:8001] # 8 s sa prvog kanala
fs <- 1000
vr <- seq(0, length(emg)/fs - 1/fs, by = 1/fs)

jpeg('EMG-Pectoralis.jpg', units='in', width=6, height=4, res=400)
par(mar = c(4.2,4.2,0.8,0.8))
plot(vr, emg, type = "l", 
     xlab = "vreme [s]", ylab = "amplituda [mV]")
grid()
dev.off()

rez2 <- emd(emg, vr, boundary = "wave") # apply EMD

jpeg('EMGEMD.jpg', units='in', width=8, height=6, res=400)
par(mfrow = c(rez2$nimf+1, 1), mar = c(0,0,0,0)) 
rangeimf <- range(rez2$imf)
for(ind in 1:rez2$nimf) {
  plot(vr, rez2$imf[, ind], type="l", 
       xlab = "", ylab= "", ylim = rangeimf,
       axes = FALSE); #abline(h = 0, lty = 2)
}
plot(vr, rez2$residue, xlab = "",
     ylab = "", type = "l",
     axes = FALSE)
dev.off()

# integracija EMG signala
dat <- read.table("EMGpectoralis.txt")
emg <- dat$V1[1:8000] # 8 s sa prvog kanala

# plot(vr, emg, type = "l")

library(pracma)
# plot(vr, cumtrapz(abs(emg)), type = "l")

max(cumtrapz(abs(emg)))
trapz(abs(emg))

length(emg)
length(cumtrapz(abs(emg)))

# određivanje početka i kraja mišićne kontrakcije
fs <- 1000
datEMGkablAll <- read.table("emg-sve.txt")

# iskorišćen je opseg signala na kome nema šuma
datEMGkabl <- datEMGkablAll$V1[(50*fs):(80*fs)]
timeAxisKabl <- seq(0, length(datEMGkabl)/fs - 1/fs, by = 1/fs)

# ispravljanje EMG signala
datEMGkabl <- abs(datEMGkabl)

# računanje obvojnice EMG signala
ma <- rep(1/200, 200)
datEMGkablMA <- stats::filter(datEMGkabl, ma, sides = 1)

# prikaz EMG signala
jpeg('EMG-onoff-1.jpg', units='in', width=8, height=6, res=400)
par(mar = c(4.2,4.2,0.8,0.8))
plot(timeAxisKabl, abs(datEMGkabl), type = "l", col = "cyan",
     xlab = "vreme [s]", ylab = "amplituda [mV]",
     ylim = c(0, 1.6))
par(new = T)
plot(timeAxisKabl, datEMGkablMA, type = "l", lwd = 2,
     xlab = "vreme [s]", ylab = "amplituda [mV]",
     ylim = c(0, 1.6))
grid()
dev.off()

# primena findpeaks() funkcije za detekciju početka i kraja kontrakcije
library(pracma)
emgObvojnica <- as.vector(datEMGkablMA)
emgObvojnica <- emgObvojnica[200:length(emgObvojnica)]
prag <- 2*sd(emgObvojnica)
emgObvojnica[emgObvojnica > prag] <- prag
onoff <- findpeaks(emgObvojnica, minpeakheight = prag/2,
                   minpeakdistance = 2*fs, nups = 0)
onoff

emgObvojnica <- as.vector(datEMGkablMA)
timeAxisKabl <- seq(0, length(emgObvojnica)/fs - 1/fs, by = 1/fs)

# prikaz detektovanih mišićnih kontrakcija
jpeg('EMG-onoff-2.jpg', units='in', width=8, height=6, res=400)
par(mar = c(4.2,4.2,0.8,0.8))
plot(timeAxisKabl, emgObvojnica, type = "l", col = "cyan",
     xlab = "vreme [s]", ylab = "amplituda [mV]",
     ylim = c(0, 0.45), lwd = 3)
abline(v = onoff[1,3]/fs, lty=2, lwd = 2)
abline(v = onoff[3,3]/fs, lty=2, lwd = 2)
abline(v = onoff[5,3]/fs, lty=2, lwd = 2)
abline(v = onoff[7,3]/fs, lty=2, lwd = 2)
abline(v = onoff[9,3]/fs, lty=2, lwd = 2)
abline(v = onoff[11,3]/fs, lty=2, lwd = 2)
abline(v = onoff[13,3]/fs, lty=2, lwd = 2)
grid()
dev.off()

# računanje medijane na EMG signalu
library(bspec)

fs <- 1000
datEMGkablAll <- read.table("emg-sve.txt")

datEMGkabl <- datEMGkablAll$V1[(50*fs):(80*fs)]
timeAxisKabl <- seq(0, length(datEMGkabl)/fs - 1/fs, by = 1/fs)

# računanje Velčovog spektra snage
datEMGkabl <- as.ts(datEMGkabl, frequency = fs)
psdw <- welchPSD(datEMGkabl, seglength = 100)

# određivanje medijane
psdw05 = sum(psdw$power)/2;
psdwmed = 0;
ind <- 1
while (psdwmed <= psdw05) {
  psdwmed = psdwmed + psdw$power[ind];
  ind <- ind + 1
}

fmed <- psdw$frequency[ind]*fs

jpeg('EMG-medijana.jpg', units='in', width=8, height=6, res=400)
par(mar = c(SS4.2,4.2,0.8,0.8))
plot(psdw$frequency*fs, psdw$power, type = "l",
     lwd = 2.5, col = "cyan", xlab = "frekvencija [Hz]",
     ylab = "magnituda [a.u.]")
abline(v = fmed, lwd = 2, lty = 2)
grid()
text(fmed + 40, 0.25, "medijana")
dev.off()

################################# POGLAVLJE 4. #################################
# 4. Eksplorativna analiza biosignala
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

dat <- read.table("primerPodataka.csv", sep = ";")
head(dat)
redniBr <- 1:length(dat$V1)

# prikaz jednostavnog grafika
jpeg(file = "ugao1.jpg", width = 6, height = 3, 
     units = 'in', res = 300)
plot(redniBr, dat$V1,
     xlab = "ID", ylab = "ugao u stepenima",
     main = "Uglovi (simulirani podaci)")
grid()
dev.off()

# prikaz svih signala iz data-frame-a
jpeg(file = "ugao2.jpg", width = 6, height = 3, 
     units = 'in', res = 300)
par(mar = c(4.2, 4.2, 1.7, 1.7)) 
plot(redniBr, dat$V1, 
     xlab = "ID", ylab = "ugao u stepenima", main = "ROM",
     xlim = c(1, 12), ylim = c(60, 140))
par(new=TRUE)
plot(redniBr, dat$V1, type = "l", lwd = 2,
     xlab = "ID", ylab = "ugao u stepenima", main = "ROM",
     xlim = c(1, 12), ylim = c(60, 140))
par(new = TRUE)
plot(redniBr, dat$V2, col = "cyan", 
     xlab = "ID", ylab = "ugao u stepenima", main = "ROM",
     xlim = c(1, 12), ylim = c(60, 140))
par(new=TRUE)
plot(redniBr, dat$V2, type = "l", col = "cyan", lwd = 2,
     xlab = "ID", ylab = "ugao u stepenima", main = "ROM",
     xlim = c(1, 12), ylim = c(60, 140))
grid()grid()colours()
legend("topright", c('skočni zglob', 'koleno'),
       col=c("black", "cyan"), lty = 1, cex = 0.8,  lwd = 2)
abline(h = 120, col = "black", lwd = 2, lty = 2)
dev.off()

# grafici u ggplot2 paketu
library(ggplot2)
example("geom_point")

# primer jednog bar grafika
ggplot(diamonds) + geom_bar(aes( x = cut, fill = clarity )) 
ggsave("barPrimer.jpg", dpi = 400, units = "in",
       width = 6, height = 4.5)

ggplot()
ggplot(diamonds)

ggplot(diamonds) + geom_bar(aes(x = clarity, fill = cut))
ggsave("barPrimer2.jpg", dpi = 400, units = "in",
       width = 6, height = 4.5)

ggplot(diamonds) + geom_bar(
  aes(x = clarity, fill = cut), width = .5) + coord_flip()
ggsave("barPrimer3.jpg", dpi = 400, units = "in",
       width = 6, height = 4.5)

# primer box plot grafika
library(ISwR)

# osnovne R funkcije
jpeg(file = "boxPlot.jpg", width = 2.7, height = 5, 
     units = 'in', res = 300)
boxplot(volume ~ method, data = lung,
        main = "Kapacitet pluća", xlab = "Metoda merenja",
        ylab = "Zapremina [L]")
grid()
dev.off()

# ggplot funkcije
g <- ggplot(lung, aes(x=method, y=volume)) +
  ggtitle("Kapacitet pluća") + xlab("Metoda merenja") +
  ylab("Zapremina [L]")
g + geom_boxplot()
ggsave("boxPlot2.jpg", dpi = 400, units = "in",
       width = 2.7, height = 5)

# dodavanje pojedinačnih merenja na grafiku
g + geom_boxplot() + geom_point() + theme_minimal()
ggsave("boxPlot3.jpg", dpi = 400, units = "in",
       width = 6, height = 4.5, bg = "white")

# violinski grafik
g + geom_violin() + theme_classic()
ggsave("violin.jpg", dpi = 400, units = "in",
       width = 6, height = 4.5, bg = "white")

# Bland-Altman grafik za poređenje dve metode
library(blandr)

ba <- blandr.statistics(lung$volume[lung$method == "A"],
                        lung$volume[lung$method == "B"])
blandr.plot.ggplot(ba, plotTitle = "Bland-Altman grafik",
                   ciDisplay = FALSE , ciShading = FALSE) +
  xlab("srednje vrednosti [L]") + ylab("razlike [L]") +
  theme_minimal()
ggsave("BA.jpg", dpi = 400, units = "in",
       width = 6, height = 4.5, bg = "white")

# Blood1 skup podataka: učitavanje i definisanje kategoričkih promenljivih
dat <- read.csv(url("https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv/Stat2Data/Blood1.csv"))
head(dat)

dat$Smoke <- ifelse(dat$Smoke == 1, "smoker", "non-smoker")
dat$Smoke = as.factor(dat$Smoke)

dat$Overwt <- ifelse(dat$Overwt == 0, "normal", dat$Overwt)
dat$Overwt <- ifelse(dat$Overwt == 1, "overweight", dat$Overwt)
dat$Overwt <- ifelse(dat$Overwt == 2, "obese", dat$Overwt)
dat$Overwt <- as.factor(dat$Overwt)

head(dat)

# prikaz violinskog grafika
p1 <- ggplot(dat, aes(x = Overwt, y = SystolicBP)) +
  ggtitle("Sistolni krvi pritisak") +
  ylab("pritisak [mmHg]") +
  xlab("BMI kategorije") + theme_minimal()
p1 + geom_violin() + geom_point(aes(col = Smoke)) +
  scale_color_manual(values = c("cyan", "black"))
ggsave("violin500.jpg", dpi = 400, units = "in",
       width = 6, height = 4.5, bg = "white")

# podešavanje boja na ggplot grafiku
library(ISwR)
library(ggplot2)

ggplot(tlc, aes(x = height, y = tlc)) +
  geom_point(aes(col = age)) +
  scale_color_continuous(name = "Starost [godine]",
                         breaks = c(11, 31, 52), 
                         labels = c(11, 31, 52),
                         low = "cyan", high = "black") +
  xlab("visina [cm]") + ylab("tlc [L]") +
  theme_minimal()
ggsave("scatterTLC.jpg", dpi = 400, units = "in",
       width = 6, height = 4.5, bg = "white")

# prikaz histograma - provera da li je raspodela Gausova?
p2 <- ggplot(dat, aes(x = SystolicBP)) + 
  xlab("pritisak [mmHg]") + ylab("normalizovan broj merenja")
p3 <- p2 + geom_histogram(aes(y = ..density..), binwidth = 10, 
                          color = "black", fill = "white") +
  geom_density(alpha = .2, fill = "cyan") + theme_minimal()
print(p3)
ggsave("histogram500.jpg", dpi = 400, units = "in",
       width = 6, height = 4.5, bg = "white")

# prikaz histograma samo za jednu kategoriju
library(dplyr)
datN <- dplyr::filter(dat, Overwt == "normal")
p4 <- ggplot(datN, aes(x = SystolicBP)) + 
  xlab("pressure [mmHg]") + ylab("normalizovan broj merenja")
p5 <- p4 + geom_histogram(aes(y = ..density..), binwidth = 10,
                          color = "black", fill = "white") +
  geom_density(alpha = .2, fill = "cyan") + theme_minimal()
print(p5)

# prikaz error bara (dijagram greške)
numNormal <- length(dat$SystolicBP[dat$Overwt == "normal"])
numObese <- length(dat$SystolicBP[dat$Overwt == "obese"])
numOverweight <- length(dat$SystolicBP[dat$Overwt == "overweight"])

meanN <- filter(dat, Overwt == "normal") %>% select(SystolicBP) %>% unlist %>% mean
meanOw <- filter(dat, Overwt == "overweight") %>% select(SystolicBP) %>% unlist %>% mean
meanOb <- filter(dat, Overwt == "obese") %>% select(SystolicBP) %>% unlist %>% mean

sdN <- filter(dat, Overwt == "normal") %>% select(SystolicBP) %>% unlist %>% sd
sdOw <- filter(dat, Overwt == "overweight") %>% select(SystolicBP) %>% unlist %>% sd
sdOb <- filter(dat, Overwt == "obese") %>% select(SystolicBP) %>% unlist %>% sd

pod <- list("numeric")
pod$BMI <- as.factor(c("normal", "obese", "overweight"))
pod$mn <- c(meanN, meanOw, meanOb)
pod$sd <- c(sdN, sdOw, sdOb)
pod <- as.data.frame(pod)
pod

barE <- ggplot(pod, aes(x = BMI, y = mn)) +
  xlab("BMI grupa") + ylab("pritisak [mmHg]") +
  geom_bar(stat = "identity", fill = "cyan")
barE + geom_errorbar(aes(ymin = mn - sd, ymax = mn + sd), width = 0.2, lwd = 1.2) +
  coord_flip(ylim = c(70, 200)) + theme_minimal()
ggsave("error500.jpg", dpi = 400, units = "in",
       width = 6, height = 3, bg = "white")

# prikaz vremenskih serija
library(ggplot2)
library(dplyr)
dat <- read.table("puls.txt")
fs <- 1000

pod <- list()
pod$vreme <- rep(seq(0, length(dat$V1)/fs - 1/fs, by = 1/fs), 2)
pod$odbirci <- c(dat$V1, dat$V2)
pod$signal <- c(rep("ekg", length(dat$V1)), 
                rep("ppg", length(dat$V2)))
pod <- as.data.frame(pod)
head(pod, 3)

ggplot(pod, aes(x = vreme, y = odbirci)) + geom_line(aes(col = signal)) +
  ylab("amplituda [a.u.]") + xlab('vreme [s]') +
  xlim(0, 10) + theme_minimal()
ggsave("ekg-ppg.jpg", dpi = 300, units = "in",
       width = 6, height = 3, bg = "white")

ggplot(pod, aes(x = vreme, y = odbirci)) + geom_line() +
  ylab("normalizovana amplituda") + xlab('vreme [s]') +
  xlim(0, 7) + facet_wrap(~signal)
ggsave("ekg-ppg-2.jpg", dpi = 300, units = "in",
       width = 6, height = 3, bg = "white")

# transformacija signala za bolje poređenje
pod$odbirci <- c(dat$V1/mean(dat$V1), ((dat$V2 - mean(dat$V2))*4000)+40)

ggplot(pod, aes(x = vreme, y = odbirci)) + geom_line() +
  ylab("normalizovana amplituda") + 
  xlab('vreme [s]') + ylim(-40, 40) +
  xlim(0, 4) + facet_wrap(~signal) + theme_minimal()
ggsave("ekg-ppg-3.jpg", dpi = 300, units = "in",
       width = 6, height = 3, bg = "white")

ggplot(pod, aes(x = vreme, y = odbirci)) + geom_line(aes(col = signal)) +
  ylab("normalizovana amplituda") +
  xlab('vreme [s]') + xlim(0, 4) + theme_minimal() + ylim(-40, 40) +
  scale_color_manual(values = c("cyan", "black"))
ggsave("ekg-ppg-4.jpg", dpi = 300, units = "in",
       width = 6, height = 3, bg = "white")

# zavisnost biomarkera
library(ISwR)
library(ggplot2)

gr <- ggplot(tlc, aes(x = height, y = tlc)) + geom_point(shape = 3) +
  xlab("visina [cm]") + ylab("tlc [L]") +
  theme_minimal()
print(gr)
ggsave("tlc-1.jpg", dpi = 300, units = "in",
       width = 6, height = 3, bg = "white")

# linearni model
gr + geom_smooth(method = lm)
ggsave("tlc-2.jpg", dpi = 300, units = "in",
       width = 6, height = 3, bg = "white")

# nelinearni model
gr + geom_smooth(method = loess)
ggsave("tlc-3.jpg", dpi = 300, units = "in",
       width = 6, height = 3, bg = "white")

# kros-korelacioni koeficijenti sa linearnim modelom
library(ggpubr)
ggscatter(tlc, x = "height", y = "tlc", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "visina [cm]", ylab = "tlc [L]")
ggsave("tlc-4.jpg", dpi = 300, units = "in",
       width = 6, height = 3, bg = "white")

# testovi za računanje kros-korelacionoih koeficijenata
cor.test(tlc$height, tlc$tlc, method = "pearson")
cor.test(tlc$height, tlc$tlc, method = "spearman")
cor(tlc$height, tlc$tlc, method = "pearson")
cor(tlc$height, tlc$tlc, method = "spearman")

# kros-korelaciona matrica
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dat <- read.table("ID1_PS_P1.txt")
corM <- cor(dat)

# osnovna funkcija u R-u za prikaz kros-korelacione matrice
image(corM)

# ggcorrplot funkcija za prikaz kros-korelacione matrice
library(ggcorrplot)
ggcorrplot(corM)
ggsave("corr1.jpg", dpi = 400, units = "in",
       width = 6, height = 3, bg = "white")

ggcorrplot(corM, type = "lower", outline.col = "white",
           colors = c("darkmagenta", "white", "gold3"))
ggsave("corr2.jpg", dpi = 400, units = "in",
       width = 6, height = 3, bg = "white")

# prikaz vrednosti kros-korelacionih koeficijenata na grafiku
ggcorrplot(corM, type = "lower", outline.col = "white",
           colors = c("darkmagenta", "white", "gold3"),
           lab = TRUE)

# interaktivna kros-korelaciona matrica
library(plotly)
ggcorrplot(corM, type = "lower", outline.col = "white",
           colors = c("darkmagenta", "white", "gold3"))
ggplotly(ggplot2::last_plot())

# prikaz akcionog potencijala muholovke - kompozitni grafik
library(jpeg)
library(ggplot2)
slika <- readJPEG("muholovka.jpg")

dat <- read.table("plant.csv")
fs <- 1000 # frekvencija odabiranja
dat$time <- seq( 0, (length(dat$V1) - 1)/fs, by = 1/fs )

# prikaz kompozitne slike (pojačanje je bilo 100 puta)
plot(dat$time, 10*dat$V1,  type = "l",
     xlab = "vreme [s]", ylab = "napon [mV]", 
     main = "Akcioni potencijal muholovke",
     xlim = c(0, 9), ylim = c(-400, 150), lwd = 2)
grid()
rasterImage(slika, 0, -400, 5, -100)

# ggplot2 paket za prikaz kompozitnog grafika
ggplot(dat, aes(x = time, y = 10*V1)) +
  geom_line() + ggtitle("Akcioni potencijal muholovke") +
  xlab("vreme [s]") + ylab("napon [mV]") +
  theme_minimal() +
  xlim(0, 9) + ylim(-400, 150) +
  annotation_raster(slika, 0, 5, -400, -100)
ggsave("muholovkaggplot.jpg", dpi = 400, units = "in",
       width = 6, height = 3, bg = "white")

# heat maps odnosno toplotne mape
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dat <- read.table("ID1_PS_P1.txt")
corM <- cor(dat)

jpeg(file = "heatmap1.jpg", width = 6, height = 7, 
     units = 'in', res = 300)
heatmap(corM)
dev.off()

mojeBoje <- colorRampPalette(c("cyan", "black"))(20)
jpeg(file = "heatmap2.jpg", width = 6, height = 7, 
     units = 'in', res = 300)
heatmap(corM, xlab = "EMG kanali",
        ylab = "EMG kanali",
        col = mojeBoje)
dev.off()

# interaktivna mapa sa barom sa bojama
library(heatmaply)
heatmaply_cor(corM)

# spektrogram EMG signala
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
dat <- read.table("ID1_PS_P1.txt")
fs <- 1000

# primena funkcije iz signal paketa
library(signal)
jpeg(file = "signal-spektrogram.jpg", width = 6, height = 3, 
     units = 'in', res = 300)
  par(mfrow = c(1, 2))
    specgram(dat$V1, n = 50, Fs = fs)
    specgram(dat$V1, n = 500, Fs = fs)
dev.off()

# primena seewave funkcije koja je kompatibilna sa ggplot paketom
library(seewave)

ggspectro(dat$V1, f = fs, tlab = "vreme [s]", wl = 500,
                   wn = "hamming", ovlp = 50,
                   flab = "frekvencija (kHz)", 
                   alab = "magnituda [a.u.]") +
  geom_tile(aes(fill = amplitude)) + theme_minimal()
ggsave("seewave-spektrogram.jpg", dpi = 400, units = "in",
       width = 6, height = 3, bg = "white")

# grafik vremenske linije
library(timevis)

dat <- data.frame(id = 1:11, content = c("internet", "MIT najava",
              "UNESCO forum: OER definicija", "MIT OpenCourseWare",
              "MIT izveštaj", "Vikiverziti", "Coursera", "edX",
              "Pravilnici", "Prvi OER na ETFu", "Platforma za otvorenu nauku"),
  start   = c("1995-01-01", "2001-01-01", "2002-01-01", "2002-01-01", "2006-01-01", 
              "2006-01-01", "2012-01-01", "2012-01-01", "2007-01-01", "2010-01-01",
              "2018-01-01"),
  end     = c(rep(NA, 11)),
  group = c(rep(1, 8), rep(2, 3)),
  style = c(rep(NA, 8), rep("background-color: white;", 3)))
            
timevis(dat, groups = data.frame(id = 1:2, content = c("svet", "ETF"),
                                 style = c("color: black; background-color: #d5ddf6;", 
                                           "color: black; background-color: white;")))

################################# POGLAVLJE 5. #################################
# 5. Statistička analiza biosignala
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# podaci sa sajta: 
# https://vincentarelbundock.github.io/Rdatasets/datasets.html
# pristupljeno 17.05.2024.
dat <- read.csv("sleepstudy.csv") 
head(dat)
dim(dat)
class(dat$rownames)
class(dat$Reaction)
class(dat$Days)
class(dat$Subject)
sum(is.na(dat$Reaction))

# alternativno, moguće je koristiti funkciju:
str(dat)

# prikaz deskriptivne statistike
summary(dat)

# prikaz Reaction parametra po danima korišćenjem dplyr paketa
library(dplyr)
filter(dat, Days == 0) %>% select(Reaction) %>% unlist %>% mean()
filter(dat, Days == 1) %>% select(Reaction) %>% unlist %>% mean()
filter(dat, Days == 9) %>% select(Reaction) %>% unlist %>% mean()

pod <- dat %>% group_by(Days) %>% 
  summarize(Average = mean(Reaction)) %>% as.data.frame()
pod

pod$SD <- dat %>% group_by(Days) %>% 
  summarize(SD = sd(Reaction)) %>% select(SD) %>% unlist()
pod

# primena Saphiro-Wilk testa za proveru da li su podaci normalno raspodeljeni
shapiro.test(dat$Reaction[dat$Days == 9])$p.value
shapiro.test(dat$Reaction[dat$Days == 8])$p.value
shapiro.test(dat$Reaction[dat$Days == 7])$p.value
shapiro.test(dat$Reaction[dat$Days == 6])$p.value
shapiro.test(dat$Reaction[dat$Days == 5])$p.value
shapiro.test(dat$Reaction[dat$Days == 4])$p.value
shapiro.test(dat$Reaction[dat$Days == 3])$p.value
shapiro.test(dat$Reaction[dat$Days == 2])$p.value
shapiro.test(dat$Reaction[dat$Days == 1])$p.value

# prikaz box plot grafika za vreme reakcije po danima
library(ggplot2)
ggplot(dat, aes(x = as.factor(Days), y = Reaction)) + 
  geom_boxplot() + xlab("Days") + ylab("Reaction [ms]")
ggsave("box_plot_sleep_study.jpg", dpi=300)

# primena t-testa za poređenje vremena reakcije 0. i 9. dana
t.test(dat$Reaction[dat$Days == 0], dat$Reaction[dat$Days == 9])
t.test(dat$Reaction[dat$Days == 9], dat$Reaction[dat$Days == 0])
t.test(dat$Reaction[dat$Days == 9], dat$Reaction[dat$Days == 0],
       paired = TRUE)

# neparametarski test
wilcox.test(dat$Reaction[dat$Days == 9], dat$Reaction[dat$Days == 0],
            paired = TRUE)

# objašnjenje t-testa
library(ISwR)
dat <- PlantGrowth
head(dat)

library(dplyr)
ctDat <- filter(dat, group == "ctrl")
t1Dat <- filter(dat, group == "trt1")
t2Dat <- filter(dat, group == "trt2")

pod <- list("numeric")
pod$mean <- c(mean(ctDat$weight), mean(t1Dat$weight),
              mean(t2Dat$weight))
pod$group <- c("control", "treatment1", "treatment2")
pod$sd <- c(sd(ctDat$weight), sd(t1Dat$weight),
            sd(t2Dat$weight))
pod <- as.data.frame(pod)

library(ggplot2)
barE <- ggplot(pod, aes(x = group, y = mean)) +
  ggtitle("Error bar za različite grupe biljaka") +
  xlab("grupa") + ylab("tezina [g]") +
  geom_bar(stat = "identity",
           fill = "cyan", col = "black")
barE + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                     width = 0.3) + coord_flip()

ggplot(dat, aes(x = group, y = weight)) +
  ggtitle("Box plot za različite grupe biljaka") +
  xlab("grupa") + ylab("tezina [g]") +
  geom_boxplot()

# kako napraviti dijagram greške sa p vrednostima za ctrl i trt2
pod2 <- list("numeric")
pod2$mean <- c(mean(ctDat$weight), mean(t2Dat$weight))
pod2$group <- c("control", "treatment2")
pod2$sd <- c(sd(ctDat$weight), sd(t2Dat$weight))
pod2 <- as.data.frame(pod2)

result <- t.test(ctDat$weight, t2Dat$weight, paired = F)$p.value

df_p_val <- data.frame(
  group1 = "control",
  group2 = "treatment2",
  label = round(result, 3),
  y.position = 7.5,
  p.adj.signif = "*"
)

barE <- ggplot(pod2, aes(x = group, y = mean)) +
  ggtitle("Error bar za dve grupe biljaka") +
  xlab("grupa") + ylab("tezina [g]") +
  geom_bar(stat = "identity",
           fill = "white", col = "black")
barE + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                     width = 0.3) + # coord_flip() + 
  add_pvalue(df_p_val, label = "p = {label}{p.adj.signif}") #,
#coord.flip = TRUE)

################################# POGLAVLJE 6. #################################
## 6.2 Kako izbeći greške (Debagovanje)
# profilisanje R koda
Rprof()
set.seed(38)
p <- rnorm(45)
Rprof(NULL) 

# učitava se sačuvana datoteka i sumira rezultat
summaryRprof()

# funkcija štoperice (ne sme se kombinovati sa Rprof funkcijama za profilisanje)
system.time({
  set.seed(38)
  p <- rnorm(4500000)
})