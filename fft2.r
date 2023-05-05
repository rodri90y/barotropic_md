rm(list=ls())

# load libraries
library(ggplot2)
library(grid)
# library(ggmap)
# library(reshape)
library(ncdf)
library(maptools)
library(plot3D)
library(animation)

path <- '~/workspace/barotropic_md/'
setwd(path)
files <- c('out_barotrp_jkb123_','out_barotrp_jkb1_',
           'out_barotrp_njkb.nc')


domain <- "160x40" #80x40_s1.5


grid1 <- open.ncdf(paste(path,files[1],domain,'.nc',sep=''))
lon <- get.var.ncdf(grid1,varid="lon")
lat <- get.var.ncdf(grid1,varid="lat")

psi <- get.var.ncdf(grid1,varid="streamfunction")
zeta <- get.var.ncdf(grid1,varid="vorticity")
u <- get.var.ncdf(grid1,varid="u")
v <- get.var.ncdf(grid1,varid="v")

nt    <- grid1$dim$time$len
nlat  <- grid1$dim$lat$len
nlon  <- grid1$dim$lon$len
dt <- 6

#cmd <- as.numeric(readline("cmd (1-field/2-conservation): "))

t  <- as.numeric(readline("t: "))

uw <- u[,,t]
vw <- v[,,t]
z  <- zeta[,,t]
p  <- psi[,,t]

dx <- dim(p)[1]
dy <- dim(p)[2]

ft2 <- mvfft(p, inverse=F)
Ck2 <- abs(ft2[2:dx/2,1:dy/2])
# image(Ck2)

ft <- fft(p, inverse=F)
Ck <- abs(ft[1:dx/2])**2

if(dx%%2==0){
  k <- seq(1,(dx-1))
}else{
  k <- seq(1,(dx-1))
}
freq <- k/dx
period <- 1./freq

#plot(Ck,typ='l')


p3 <- ggplot() +
  geom_line(aes(x=k,y=Ck),size=.20) +
  labs(y="Amplitudes espectrais", x="L (1,25 10⁵ m)") +  
  scale_x_continuous(limits=c(0, 160),breaks=seq(0,160,10)) +
  ggtitle(paste('time: ',t)) +
  theme_bw(base_size=12, base_family="Helvetica")  +
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"),legend.position="bottom")
print(p3)
dev.copy(png,paste(path,"img/baro_J123_AmplSpectral",domain,"_t",t,".png",sep=''),width=10, height=5,units="in",res=150)
dev.off()


