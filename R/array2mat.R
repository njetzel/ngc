#################################################################################
#            SOME MATRIX MANIPULATION FUNCTIONS
#################################################################################
adj_right2left = function(adjmat){
ntime = dim(adjmat)[1]
p = dim(adjmat)[2]
adjout = array(0, c(p, p, ntime))
for (tt in 1:ntime)
  for (i in 1:p)
    for (j in 1:p)
      adjout[i,j,ntime-tt+1] = adjmat[tt, i,j]  #p-j+1

return(adjout)
}


adj_left2right = function(adjmat){
ntime = dim(adjmat)[3]
p = dim(adjmat)[2]
adjout = array(0, c(ntime, p, p))
for (tt in 1:ntime)
  for (i in 1:p)
    for (j in 1:p)
      adjout[tt, i,j] = adjmat[i,j,ntime-tt+1]   #p-j+1

return(adjout)
}

MatImage <- function(AA, maintitle=NULL, margin=c(2,1,2,1)) {
	d1 = dim(AA)[1]
	d2 = dim(AA)[2]

	par(mar=margin)
	image( t(AA[d1:1,]) , col=gray(100:0/100), axes=FALSE )
	box()
	title(maintitle, cex.main = 2,   font.main= 1)

	axis(2, at=c(d1,1)/d1, label=c(d1,1))

	if (d1 < d2){
		axis(1, at=seq(d2,1,by=-d1)/d2, label=seq(1,d2,by=d1))
	}else{
		axis(1, at=c(d2,1)/d2, label=c(1,d2))
	}
}

array2mat <-
function(
	myArray,		#input array
	xtnd=2		#dimension to extend the array over
){
####### START OF FN #######
	d1 = dim(myArray)[1]
	d2 = dim(myArray)[2]
	d3 = dim(myArray)[3]

	if (is.na(d3)){
		myMat = myArray
		return(myMat)
	}

	if (xtnd==1){
		myMat = matrix(0,d1*d3,d2)
	}else{
		myMat = matrix(0,d1,d2*d3)
	}

	for(i in 1:d3){
		low1 = (i-1)*d1+1
		high1 = i*d1
		low2 = (i-1)*d2+1
		high2 = i*d2

		if (xtnd==1){
			myMat[low1:high1,] = myArray[,,i]
		}else{
			myMat[,low2:high2] = myArray[,,i]
		}
	}

	return(myMat)
}

