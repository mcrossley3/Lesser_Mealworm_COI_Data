library(vegan) #for rarefaction curves
library(viridis)
setwd("path/to/data")

# Get vectors of species and sites
species = c('H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12')
sites = c('DE Farm 1 House 1','DE Farm 1 House 2','DE Farm 1 House 3','DE Farm 1 House 4','DE Farm 2 House 1','DE Farm 2 House 2','AL Farm 1 House 1','AL Farm 1 House 2','AL Farm 2','AL Farm 3','Czech Republic','Georgia','Italy','MD Farm 1','MD Farm 2','Turkey')

mat1 = matrix(data=c(16,19,17,18,13,18,19,18,20,17,11,17,14,13,18,6,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
              ,nrow=length(sites),ncol=length(species)) #empty matrix to be filled with species (column) counts at each site (row)
rownames(mat1) = sites #add rownames
colnames(mat1) = Haplotypes #add column names

sort(rowSums(mat1)) #sort sites by beetle abundance
sort(colSums(mat1)) #sort least to most abundant species


png('Haplotype accumulation curves.png',res=300,height=480*4,width=480*4) #open graphics device for drawing
pal1 = viridis(length(sites)) #color palette
par(oma=c(0,0,0,0),mar=c(5,5,5,1)) #set some graphical parameters
rarecurve(mat1,step=1,sample=NULL,col=pal1,cex=1,lwd=4,label=F,lty=rep(c(1,3),length(sites)/2),main='Accumulation of Alphitobius Haplotypes',ylab='Haplotypes') #draw rarefaction curves
legend('topleft',legend=sites,col=pal1,lwd=3,lty=rep(c(1,3),length(sites)/2),ncol=2,cex=0.5) #add legend
# Trim x-axis to better visualize species accumulation at sites with lower beetle abundance
dev.off() #close graphics device for viewing

