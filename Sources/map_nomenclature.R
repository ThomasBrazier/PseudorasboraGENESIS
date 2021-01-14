scalebar <- function(loc,length,unit="km",division.cex=0.8,bg="white",border="black",...) {
  lablength <-length
  length<-length/111
  x <- c(0,length/c(4,2,4/3,1),length*1.1)+loc[1]
  y <- c(0,length/(10*3:1))+loc[2]
  rect(x[1]-diff(x)[1]/4,y[1]-(y[2]-y[1]),x[5]+strwidth(paste(round(x[5]*111-loc[1]*111,0),unit))/2+diff(x)[1]/4,y[4]+(y[4]-y[1])/2, col=bg,border=border)
  cols <- rep(c("black","white"),2)
  for (i in 1:4) rect(x[i],y[1],x[i+1],y[2],col=cols[i])
  for (i in 1:5) segments(x[i],y[2],x[i],y[3])
  labels <- round(x[c(1,3)]*111-loc[1]*111,0)
  labels <- append(labels, paste(round(x[5]*111-loc[1]*111,0),unit))
  text(x[c(1,3,5)],y[4],labels,adj=0.5,cex=division.cex)
}

north.arrow = function(x, y, h,lab="North",lab.pos="below") {
  polygon(c(x, x, x + h/2), c(y - (1.5*h), y, y - (1 + sqrt(3)/2) * h), col = "black", border = NA)
  polygon(c(x, x + h/2, x, x - h/2), c(y - (1.5*h), y - (1 + sqrt(3)/2) * h, y, y - (1 + sqrt(3)/2) * h))
  if(lab.pos=="below") text(x, y-(2.5*h), lab, adj = c(0.5, 0), cex = 1)
  else text(x, y+(0.25*h), lab, adj = c(0.5, 0), cex = 1.5)
}