library(MARSS)
library(ggplot2)
library(forecast)
library(dplyr)
library(lubridate)
library(mgcv)
#### Reading in Data ####
data <-read.csv('Oceanography_DailyMean.csv')
head(data)
unique(data$site)
marine_reserves<-unique(data$marine_reserve)
unique(data$location)

####Herring####

yy=data%>%
  #filter(marine_reserve=='Otter Rock'| marine_reserve=='Redfish Rocks')%>% #filtering for reserves that we have SMURF data
  mutate(depth_bin=ifelse(depth_m==1,'surface',ifelse(depth_m==5|depth_m==7.5,"mid","bottom")))%>% # binning 5 and 7.5 depth to a single "mid" depth
 # group_by(depth_bin,date,julian,year,marine_reserve)%>%
  #summarise(temp=mean(temp_c))%>%
  ungroup()
yy$date<-as.Date(yy$date)
ggplot(yy, aes(x = date, y = temp_c, group=as.factor(depth_bin), col=as.factor(depth_bin)))+
  geom_line()+ 
  facet_wrap(~marine_reserve)

for(i in 1:length(marine_reserve)){
p<-ggplot(yy%>%filter(marine_reserve==marine_reserves[i]), aes(x = julian, y = temp_c, group=as.factor(depth_bin), col=as.factor(depth_bin)))+
    geom_line()+ 
    ggtitle(marine_reserve[i])+
    facet_wrap(~year)
plot(p)

}
yy<-yy%>%mutate(locationdepth=paste(location, depth_bin))

##### creating a continuouse covariate #####
dates<-seq(as.Date(min(yy$date)), as.Date(max(yy$date)), "days")
years<-lubridate::year(dates)
months<-lubridate::month(dates)
depth_bins<-c('bottom', 'mid','surface')
covariates_df<-data.frame(date=as.Date(rep(dates,3)),
           depth_bin=rep(depth_bins,each=length(dates)),
           year=years, months=months, julian=yday(dates)
)

covariates<-covariates_df%>%filter(julian>100&julian<275)
yyy<-covariates%>%left_join(yy%>%mutate(date=as.Date(date)))
unique(yy$depth_bin)
dat <- reshape2::acast(yyy%>%filter(marine_reserve=="Otter Rock"&julian>100&julian<275), depth_bin ~ date, value.var = "temp")
the.mean <- apply(dat, 1, mean, na.rm = TRUE)
the.sigma <- sqrt(apply(dat, 1, var, na.rm = TRUE))
dat <- (dat - the.mean) * (1/the.sigma)
covariates <- reshape2::acast(yyy%>%filter(marine_reserve=="Otter Rock"&julian>100&julian<275), depth_bin ~ date, value.var = "year")
#the.mean <- apply(covariates, 1, mean, na.rm = TRUE)
#the.sigma <- sqrt(apply(covariates, 1, var, na.rm = TRUE))
#covariates <- (covariates - the.mean) * (1/the.sigma)

d <-covariates
B <- "identity"
Q <- 'equalvarcov'
#"equalvarcov"
R <- "diagonal and equal"
U <- "unequal"
A <- "zero"
x0 <- "unequal"
D <- "unconstrained"
y <- dat  # to show relationship between dat & the equation
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, D = D, 
                   d = d, x0 = x0)
mod.list<-list(B = matrix(1), U = matrix("u"), Q = matrix("q"), 
     Z = matrix(1, 3, 1), A = "scaling", R = "diagonal and unequal", 
     x0 = matrix("mu"), tinitx = 0)
kem <- MARSS(dat, model = mod.list)

acf(dat, main = "flat level v(t)", na.action = na.pass)
library(broom)
library(ggplot2)
autoplot(kem, plot.type="fitted.ytT")
autoplot(fit)



par(mfrow = c(2, 2), mar = c(2, 2, 4, 2))
resids <- residuals(kem.0)
acf(resids$.resids, main = "flat level v(t)", na.action = na.pass)
resids <- residuals(kem.1)
acf(resids$.resids, main = "linear trend v(t)", na.action = na.pass)
resids <- residuals(kem.2)
acf(resids$.resids, main = "stoc level v(t)", na.action = na.pass)





d <- augment(fit, interval = "confidence")
#d$Year <- d$t + 1980
d$Station <- d$.rownames
p <- ggplot(data = d) + geom_line(aes(t, exp(.fitted))) + geom_ribbon(aes(x = t, 
                                                                          ymin = exp(.conf.low), ymax = exp(.conf.up)), linetype = 2, alpha = 0.5)
#p <- p + geom_point(data = yy, mapping = aes(x = Year, y = SWE))
p + facet_wrap(~Station) + xlab("") + ylab("SWE (demeaned)")



Fit.1 <- fit
Fit.1$states.se
tot.herring.bio <- colSums(exp(Fit.1$states))
value<- tot.herring.bio
herring.tot <- cbind(yr=seq(1973,2012,1),value)
write.csv(herring.tot, "Data/Compiled/Prey/herring.tot.csv")
pdf(file="Results/Figures/HerringMarss.pdf", width=12, height=11)
p + facet_wrap(~Station) + xlab("") + ylab("SWE (demeaned)")
dev.off()

pdf(file="Results/Figures/HerringRaw.pdf", width=12, height=11)
p <- ggplot(yy, aes(x = YEAR, y = Biomass)) + geom_line()
p + facet_wrap(~Group)
dev.off()

####Harbor Seal####

rm(list = ls()) 
data <- read.csv("Data/Processed/HarborSeal.csv")
yy=data

pdf(file="Results/Figures/HarborSealRaw.pdf", width=12, height=11)
p <- ggplot(yy, aes(x = Year, y = Count)) + geom_line()
p + facet_wrap(~Group)
dev.off()

ns <- length(unique(yy$Group))
B <- "identity"
Q <- 'equalvarcov'
#"equalvarcov"
R <- "diagonal and equal"
U <- "unequal"
A <- "zero"
x0 <- "unequal"

dat <- reshape2::acast(yy, Group ~ Year, value.var = "Count")


mod.list = list(B = B, Q = Q, R = R, U = U, x0 = x0, A = A)
dat.1 <- dat
dat.1[dat.1==0]<- NA

m <- apply(log(dat.1), 1, mean, na.rm = TRUE)
fit <- MARSS(log(dat.1), model = mod.list, control = list(maxit = 5000), 
             inits = list(A = matrix(m, ns, 1)))




library(broom)
library(ggplot2)
d <- augment(fit, interval = "confidence")
#d$Year <- d$t + 1980
d$Station <- d$.rownames

pdf(file="Results/Figures/HarborSealMARSS.pdf", width=12, height=11)
p <- ggplot(data = d) + geom_line(aes(t, exp(.fitted))) + geom_ribbon(aes(x = t, 
                                                                          ymin = exp(.conf.low), ymax = exp(.conf.up)), linetype = 2, alpha = 0.5)
#p <- p + geom_point(data = yy, mapping = aes(x = Year, y = SWE))
p + facet_wrap(~Station) + xlab("") + ylab("SWE (demeaned)")
dev.off()


Fit.1 <- fit
Fit.1$states.se
tot.harborseal.pop <- colSums(exp(Fit.1$states))
value<- tot.harborseal.pop
value <- c(value,rep(value[25], 16))
seal.tot <- cbind(yr=seq(1975,2015,1),value)
write.csv(seal.tot, "Data/Compiled/Prey/seal.tot.csv")
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(3,5,1,1))
plot(seal.tot[,1], seal.tot[,2], ylab="Population", xlab="Year", type='l')


####Chinook####

data <- read.csv("Data/Compiled/Prey/Salmon/SPS_Download_MAR102020_Chinook2.csv")
yy <- cbind(Year=data[,1], Count =rowSums(data[,2:28]))
write.csv(yy, "Data/Compiled/Prey/chinook.total.csv")

data <- read.csv("Data/Compiled/Prey/Salmon/SPS_Download_MAR102020_Chinook3.csv")
yy=data

dat <- reshape2::acast(yy, Group ~ Year, value.var = "Spawners")


mod.list = list(B = B, Q = Q, R = R, U = U, x0 = x0, A = A)
dat.1 <- dat
dat.1[dat.1==0]<- NA

m <- apply(log(dat.1), 1, mean, na.rm = TRUE)
fit <- MARSS(log(dat.1), model = mod.list, control = list(maxit = 5000), 
             inits = list(A = matrix(m, ns, 1)))




####Chum####

data <- read.csv("Data/Compiled/Prey/Salmon/SPS_Download_MAR102020_Chum.csv")
yy <- cbind(Year=data[4:46,1], Count =rowSums(data[4:46,2:3]))
write.csv(yy, "Data/Compiled/Prey/chum.total.csv")

data <- read.csv("Data/Compiled/Prey/Salmon/SPS_Download_MAR102020_Chum2.csv")

pdf(file="Results/Figures/ChumRaw.pdf", width=12, height=11)
p <- ggplot(data, aes(x = Year, y = Spawners)) + geom_line()
p + facet_wrap(~Group)
dev.off()

####Coho####

data <- read.csv("Data/Compiled/Prey/Salmon/SPS_Download_MAR102020_Coho.csv")


library(broom)
library(ggplot2)
d <- augment(fit, interval = "confidence")
#d$Year <- d$t + 1980
d$Station <- d$.rownames

pdf(file="Results/Figures/CohoMARSS.pdf", width=12, height=11)
p <- ggplot(data = d) + geom_line(aes(t, exp(.fitted))) + geom_ribbon(aes(x = t, 
                                                                          ymin = exp(.conf.low), ymax = exp(.conf.up)), linetype = 2, alpha = 0.5)
#p <- p + geom_point(data = yy, mapping = aes(x = Year, y = SWE))
p + facet_wrap(~Station) + xlab("") + ylab("SWE (demeaned)")
dev.off()

Fit.1 <- fit
Fit.1$states.se
tot.coho.pop <- colSums(exp(Fit.1$states))
value<- tot.coho.pop
coho.tot <- cbind(yr=seq(1957,2013,1),value)
write.csv(coho.tot, "Data/Compiled/Prey/coho.tot.csv")
par(mfrow=c(1,1), oma=c(0,0,0,0), mar=c(3,5,1,1))
plot(coho.tot[,1], coho.tot[,2], ylab="Population", xlab="Year", type='l')

pdf(file="Results/Figures/CohoRaw.pdf", width=12, height=11)
p <- ggplot(data, aes(x = Year, y = Spawners)) + geom_line()
p + facet_wrap(~Group)
dev.off()

pdf(file="Results/Figures/CohoRaw.pdf", width=12, height=11)
p <- ggplot(data, aes(x = Year, y = Spawners)) + geom_line()
p + facet_wrap(~Group)
dev.off()


####Wild Production####

data <- read.csv("Data/Processed/WildProductionChascoWASHINGTONcondensed.csv")
data<- subset(data, Year>=1973 & Year<=2013)
yy=data

p <- ggplot(yy, aes(x = Year, y = escapement)) + geom_line()
p + facet_wrap(~Group)

ns <- length(unique(yy$Group))
B <- "identity"
Q <- 'equalvarcov'
#"equalvarcov"
R <- "diagonal and equal"
U <- "unequal"
A <- "zero"
x0 <- "unequal"

dat <- reshape2::acast(yy, Group ~ Year, value.var = "escapement")


mod.list = list(B = B, Q = Q, R = R, U = U, x0 = x0, A = A)
dat.1 <- dat
dat.1[dat.1==0]<- NA

m <- apply(log(dat.1), 1, mean, na.rm = TRUE)
fit <- MARSS(log(dat.1), model = mod.list, control = list(maxit = 5000), 
             inits = list(A = matrix(m, ns, 1)))




library(broom)
library(ggplot2)
d <- augment(fit, interval = "confidence")
#d$Year <- d$t + 1980
d$Station <- d$.rownames
p <- ggplot(data = d) + geom_line(aes(t, exp(.fitted))) + geom_ribbon(aes(x = t, 
                                                                          ymin = exp(.conf.low), ymax = exp(.conf.up)), linetype = 2, alpha = 0.5)
#p <- p + geom_point(data = yy, mapping = aes(x = Year, y = SWE))
p + facet_wrap(~Station) + xlab("") + ylab("SWE (demeaned)")


Fit.1 <- fit
Fit.1$states.se
tot.wildproduction <- colSums(exp(Fit.1$states))
value<- tot.wildproduction
wildproduction.tot <- cbind(yr=seq(1973,2013,1),value)

write.csv(wildproduction.tot, "Data/Compiled/Prey/wildproduction.tot.csv")
pdf(file="Results/Figures/WildProductionMarss.pdf", width=12, height=11)
p + facet_wrap(~Station) + xlab("") + ylab("SWE (demeaned)")
dev.off()

pdf(file="Results/Figures/WildProduction.pdf", width=12, height=11)
p <- ggplot(yy, aes(x = Year, y = escapement)) + geom_line()
p + facet_wrap(~Group)
dev.off()