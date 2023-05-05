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
files <- c('out_barotrp_jkb123','out_barotrp_jkb1'
           )


domain <- '' #'80x40' #80x40_s1.5


grid1 <- open.ncdf(paste(path,files[1],domain,'.nc',sep=''))
grid2 <- open.ncdf(paste(path,files[2],domain,'.nc',sep=''))
lon <- get.var.ncdf(grid1,varid="lon")
lat <- get.var.ncdf(grid1,varid="lat")

psi <- get.var.ncdf(grid1,varid="streamfunction")
zeta <- get.var.ncdf(grid1,varid="vorticity")
u <- get.var.ncdf(grid1,varid="u")
v <- get.var.ncdf(grid1,varid="v")

zeta2 <- get.var.ncdf(grid2,varid="vorticity")

nt    <- grid1$dim$time$len
nlat  <- grid1$dim$lat$len
nlon  <- grid1$dim$lon$len
dt <- 48
io_dt <- 6


cmd <- as.numeric(readline("cmd (1-field/2-conservation): "))

if(cmd==1){
  cmd <- readline("animation? (n/y): ")
  
#  saveGIF({   

 for(t in seq(1,nt+1,dt)){ 
#  t<-as.numeric(readline("t: "))
    
  uw <- u[,,t]
  vw <- v[,,t]
  z <- zeta[,,t]
  p <- psi[,,t]
  
  df <- data.frame(lat=as.numeric(),
                   lon=as.numeric(), 
                   dlat=as.numeric(), 
                   dlon=as.numeric(),
                   zeta=as.numeric(),
                   psi=as.numeric()
                   )
  
  n <- nlat*nlon
  
  
  for(i in 1:nlat){
  
    for(j in 1:nlon){
      df <- rbind(df, data.frame(lat=lat[i], lon=lon[j],dlat=uw[j,i], dlon=vw[j,i], zeta=z[j,i], psi=p[j,i]))
  
    }
  }
  
  
  ft_s <- 1.5e1

  
  dlim <- max(abs(min(df$psi)),abs(max(df$psi)))

  jetcolors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  p1 <- ggplot(df, aes(x=lon, y=lat)) +
        geom_tile(aes(fill=psi),trans = 'log') +
         scale_fill_gradientn(name="PSI", colours=jetcolors(12),na.value="transparent",limits=c(-dlim,dlim)
                             ) + #limits=c(-3e4,3e4)
        scale_x_continuous(expand = c(0,0)) +
        scale_y_continuous(expand = c(0,0)) +
        scale_size_manual(values=c(1,2)) +
        geom_segment(aes(xend=lon+dlon*ft_s, yend=lat+dlat*ft_s), arrow=arrow(length = unit(0.1,"cm"))) +
        ggtitle(paste('time step:',t)) +
        theme(plot.margin=unit(c(0.4,0.5,0.1,0.1), "cm"),legend.position="right")
  
if(cmd=="y"){

  print(p1)

}else
{
  print(p1)
  dev.copy(png,paste(path,"img/baro_J123",domain,"_t",t,".png",sep=''),width=10, height=5,units="in",res=150)
  dev.off()
}
  
 
  # image2D(data=df,psi, clab=c("psi"),cex.main=0.95,
  #         rasterImage=F,
  #         contour=T, # col=gray(0:8 / 8),
  #         colkey = list(length = .75, shift = -0.05, width = 0.5, cex.axis = 0.6),# zlim = c(-1, 1)
  #)
  
 }
 
#  }, movie.name = paste('baro_vort_uv_jkb123_',domain,'.gif',sep=''), interval = 0.4, ani.width = 960, ani.height = 480)




}else if(cmd==2){

f <- 2.*7.292e-05*sin(60*pi/180.)
  
  df <- data.frame(t=as.numeric(),
                   entp1=as.numeric(),
                   entp2=as.numeric(),
                   vrt1=as.numeric(),
                   vrt2=as.numeric(),
                   k1=as.numeric(),
                   k2=as.numeric()
  )
  

for(t in seq(1,nt,1)){ 
  
  v1<-sum(zeta[,,t])
  v2<-sum(zeta2[,,t])
  
  e1<-sum((zeta[,,t])^2)
  e2<-sum((zeta2[,,t])^2)
  
  k1<-sum(0.5*abs(zeta[,,t])^2)
  k2<-sum(0.5*abs(zeta2[,,t])^2)
  
  df <- rbind(df, data.frame(t=t, vrt1=v1, vrt2=v2, entp1=e1, entp2=e2, k1=k1, k2=k2) )

}


cols <- c("(J1+J2+J3)/3"="black","J1"="red")

p2 <- ggplot() +
  geom_line(data=df, aes(x=t*io_dt,y=entp1,group=1,color='(J1+J2+J3)/3'),size=0.7) + 
  geom_line(data=df, aes(x=t*io_dt,y=entp2,group=1,color='J1'),size=.7) +
  labs(y = "Enstrofia", x = "Time steps") +
  ggtitle("Quadrado da vorticidade") +
#   scale_x_continuous(limits=c(0, 150),breaks=seq(0,150,10)) +
#   scale_y_continuous(limits=c(7.319e-11, 7.320e-11)) +
  theme_bw(base_size=12, base_family="Helvetica")  +
  scale_colour_manual(name="Esquemas",values=cols) +
  scale_size_manual(values=c(1,2)) +
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),legend.position="bottom")
print(p2)
dev.copy(png,paste(path,"img/baro_J123_enstro_",domain,".png",sep=''),width=12, height=8,units="in",res=150)
dev.off()

p3 <- ggplot() +
  geom_line(data=df, aes(x=t*io_dt,y=k1,group=1,color='(J1+J2+J3)/3'),size=0.7) + 
  geom_line(data=df, aes(x=t*io_dt,y=k2,group=1,color='J1'),size=.7) +
  labs(y = "Energia cinética", x = "Time steps") +
  ggtitle("Energia cinética") +
  theme_bw(base_size=12, base_family="Helvetica")  +
  scale_colour_manual(name="Esquemas",values=cols) +
  scale_size_manual(values=c(1,2)) +
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),legend.position="bottom")
print(p3)

 dev.copy(png,paste(path,"img/baro_J123_energ_",domain,".png",sep=''),width=12, height=8,units="in",res=150)
 dev.off()

p4 <- ggplot() +
  geom_line(data=df, aes(x=t*io_dt,y=vrt1,group=1,color='(J1+J2+J3)/3'),size=0.7) + 
  geom_line(data=df, aes(x=t*io_dt,y=vrt2,group=1,color='J1'),size=.7) +
  labs(y = "Vorticidade", x = "Time steps") + 
  #scale_y_continuous(limits=c(0, -1e-2)) +
  ggtitle("vorticidade") +
  theme_bw(base_size=12, base_family="Helvetica")  +
  scale_colour_manual(name="Esquemas",values=cols) +
  scale_size_manual(values=c(1,2)) +
  theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),legend.position="bottom")
print(p4)

 dev.copy(png,paste(path,"img/baro_J123_vort_",domain,".png",sep=''),width=12, height=8,units="in",res=150)
 dev.off()


}


