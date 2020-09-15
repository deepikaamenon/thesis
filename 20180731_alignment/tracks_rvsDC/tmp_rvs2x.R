
w1_files = dir( pattern = 'w1.txt$' )
w2_files = dir( pattern = 'w2.txt$' )

pdf( 'tmp.pdf' )

for ( i in seq( length( w1_files ) ) )  {

	w1 = read.table( w1_files[ i ] , comment.char = "%" )
	w2 = read.table( w2_files[ i ] , comment.char = "%" )
	
	w2_t=w2[,1]-w2[,1][1]
	w2_fi=w2[,4]-w2[,4][4]
	
	w1_t=w1[,1]-w2[,1][1]
	w1_fi=w1[,4]-w1[,4][4]
	
	plot( w2_t , w2_fi , col = 'green' , t = 'l' , xlim=c(-20,150),ylim =c(-1,3), xlab = "frame", ylab = "fl.int" )
	lines( w1_t , w1_fi , col = 'red' )
	title(w1_files[i])
}

dev.off()
#xlim=c(-5,25), ylim=c(0,1.2))