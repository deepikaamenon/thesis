# 17.09.2012
# a rewritten and optimized version of the code used to align SC tracks usgin DW tracks, or DW tracks alone if no reference track is given
library(MASS)

#tracks can have missing time points, introduceNA find them and replaces them with NA
introduceNA<-function(data){
	#create the frame numbers, including those missing 
	framesn<-seq(data[1,1],data[length(data[,1]),1],by=1)
	#create a matrix of framen (first column) and NA numbers for all the other coulumns of data
	M<-cbind(framesn,matrix(NA,length(framesn),length(data[1,])-1))
	#pick the rows of M which are in data...
	subst<-pmatch(as.integer(data[,1]),M[,1])
	#... and assign data to them
	M[subst,]<-as.matrix(data)
	#return M with NA for the points that were not in data
	return(M)
}
#tracks can have missing time points, introduceNA find them and replaces them with NA
#introduceNA<-function(data){
#	output<-c();
#	Ddata<-data[-1,1]-data[-length(data[,1]),1]
#	gaps<-cbind(which(Ddata>1),Ddata[which(Ddata>1)])
#	if (length(gaps[,1])>0){
#		old_gap_end<-1;
#		for (i in 1:length(gaps[,1])){
#			output<-rbind(output,data[old_gap_end:gaps[i,1],]);
#			m<-cbind(seq(1,gaps[i,2]-1)+data[gaps[i,1],1],matrix(NA,gaps[i,2]-1,length(data[1,])-1))
#			output<-rbind(output,m);
#			old_gap_end<-gaps[i,1]+gaps[i,2]-1;
#		}
#		return(rbind(output,data[old_gap_end:length(data[,1]),]));
#	} else return(data)
#}
#import tracks as ts
importts<-function(path,dt,comment.char="%",na.put=TRUE,five.cols=TRUE){
	data<-read.table(path,comment.char=comment.char);
	if (na.put) output<-introduceNA(data)
	else output<-data
	t0=(output[1,1])*dt
	names<-dimnames(output)[[2]]
	names[c(2,3,4)]<-c("x","y","FI")
	dimnames(output)<-list(dimnames(output)[[1]],names)
	if (five.cols) return(ts(output[,2:6],start=t0, deltat=dt))
	else return(ts(output[,2:length(output[1,])],start=t0, deltat=dt))
}

#center of mass: center a track on it's center of mass
center_mass<-function(ts){
	ts[,1]<-ts[,1]-mean(ts[,1],na.rm=TRUE)
	ts[,2]<-ts[,2]-mean(ts[,2],na.rm=TRUE)
	return(ts)
}
	
#load.data load SC data
#load.data<-function(path,center.mass=FALSE){
load.data<-function(path){

	d<-get(load(path))
	t<-window(d$ts,start=d$start["mean"],end=d$end["mean"])
	
	return(t)
}

#cross correlation of time series having different time steps. ref is the second track, candidate to align the first (the oppsite of the first version)
cc<-function(ts1,ts2,whichFI.1=3,whichFI.2="mean_sFI",filter.span=5){

	#filtering of the data, NAs, if any, are discarded as they would propagate 2*filter.span - 1 times	
	ts1[!is.na(ts1[,whichFI.1]),whichFI.1]<-filter(ts1[!is.na(ts1[,whichFI.1]),whichFI.1],rep(1/filter.span,filter.span))
	ts2[!is.na(ts2[,whichFI.2]),whichFI.2]<-filter(ts2[!is.na(ts2[,whichFI.2]),whichFI.2],rep(1/filter.span,filter.span))

	#filtering introduce NAs at the beginning and end of the track in the FI, which need to be removed before the spline
	ts1.t<-time(ts1)[!is.na(ts1[,whichFI.1])]
	ts2.t<-time(ts2)[!is.na(ts2[,whichFI.2])]
	
	#ts1 end will match ts2 start first, this is also good because we wanna discretize the two coourves in the same intervals
	lag0<-ts2.t[1]-ts1.t[length(ts1.t)]	
	ts1.t<-ts1.t+lag0
	
	#choose the smallest dt avaliable and generate the splines	
	if (tsp(ts1)[3]>tsp(ts2)[3]) dt<-1/tsp(ts1)[3]
	else 	dt<-1/tsp(ts2)[3]
	f1<-spline(ts1.t,ts1[!is.na(ts1[,whichFI.1]),whichFI.1],xout=rev(seq(ts1.t[length(ts1.t)],ts1.t[1],-dt)))
	f2<-spline(ts2.t,ts2[!is.na(ts2[,whichFI.2]),whichFI.2],xout=seq(ts2.t[1],ts2.t[length(ts2.t)],dt))
	#compute the cross correlation
	cc.output<-c()
	i<-0
	fs<-sort(c(min(f1$x),max(f1$x),min(f2$x),max(f2$x)))
	while (!min(f1$x)>max(f2$x)){
		tmatch1<-which(f1$x>=fs[2] & f1$x<=fs[3])
		tmatch2<-which(f2$x>=fs[2] & f2$x<=fs[3])
		if (isTRUE(all.equal(f1$x[tmatch1[1]-1],fs[2]))) tmatch1<-c(tmatch1[1]-1,tmatch1)
		if (isTRUE(all.equal(f1$x[tmatch1[length(tmatch1)]+1],fs[3]))) tmatch1<-c(tmatch1[length(tmatch1)]+1,tmatch1)
		if (isTRUE(all.equal(f2$x[tmatch2[1]-1],fs[2]))) tmatch2<-c(tmatch2[1]-1,tmatch2)
		if (isTRUE(all.equal(f2$x[tmatch2[length(tmatch2)]+1],fs[3]))) tmatch2<-c(tmatch2[length(tmatch2)]+1,tmatch2)
		if (length(tmatch1) != length(tmatch2)){
			save(f1,f2,tmatch1,tmatch2,fs,file="try.Rdata")
			stop("error in defining the overlap of the curves")
		}
		cc.output<-rbind(cc.output,c(dt*i,sum(f1$y[tmatch1]*f2$y[tmatch2])))
		i<-i+1
		f1$x<-f1$x+dt
		fs<-sort(c(min(f1$x),max(f1$x),min(f2$x),max(f2$x)))
	}
	#llag=lag0+cc.output[which(cc.output[,2]==max(cc.output[,2])),1]
	#save(ts1,ts2,f1,f2,lag0,cc.output,ts1.t,ts2.t,llag,dt,file="tmp3.Rdata");stop()
	return(list(lag=lag0+cc.output[which(cc.output[,2]==max(cc.output[,2])),1],cc=cc.output))
}
#compute the weighted mean

wmean<-function(x,w,na.rm=TRUE){
	#na.rm because there might be less points than weight (some points might be NA while the weight not)
	return(sum(x*w,na.rm=na.rm)/sum(w[!is.na(x)],na.rm=na.rm))
}

Rmatrix<-function(theta){	
	R<-matrix(c(
	cos(theta),-sin(theta),
	sin(theta),cos(theta)
	),nrow=2,ncol=2,byrow=TRUE)
	R
}

rototrasl<-function(x,R,trasl=c(0,0),transf.err=c()){

	#if the track has std_x and std_y propagate the error of the transformation if there
	if (!is.null(transf.err) & length(grep("std_x",dimnames(x)[[2]]))){
		sR<-R^2
		sx<-x[,"std_x"]/sqrt(x[,"n"])
		sy<-x[,"std_y"]/sqrt(x[,"n"])
		theta.err<-cbind(transf.err[1]*x[,"mean_y"],transf.err[1]*x[,"mean_x"])
		trasl.err<-cbind(rep(transf.err[3],length(x[,"mean_x"])),rep(transf.err[4],length(x[,"mean_y"])))
		x[,c("std_x","std_y")]<-sqrt(t(sR%*%t(cbind(sx^2,sy^2)))+t(sR%*%t(theta.err^2))+trasl.err^2)
	}	
	x[,c(1,2)]<-t(R%*%t(x[,c(1,2)])+as.vector(trasl))
	return(x)
}
	

#trackalign aligns a track to its template
trackalign<-function(track,template,doweight=TRUE,input.filter=4,filter.method="convolution"){
	#tracks are assumed to be axis simmetric. Filter the track to reduce the noise and interpolate it with a spline. This preserves the main features of the track, spline dt will then match the template dt.
	deltat<-1/tsp(template)[3]

	tmp.names<-dimnames(track)[[2]]
	track<-cbind(time(track)[!is.na(track[,"x"])],track[!is.na(track[,"x"]),])
	dimnames(track)<-list(c(),c("t",tmp.names))
	tmp.names<-dimnames(template)[[2]]
	template<-cbind(time(template)[!is.na(template[,"mean_x"])],template[!is.na(template[,"mean_x"]),])
	dimnames(template)<-list(c(),c("t",tmp.names))

	if (is.null(input.filter)){
		#select only the matching points in the two tracks. Drawbacks: very few points for the alignment and most worse the points in the template do not match the exact time of those in the track
		#match time points
		sel.template<-c()
		for (t in track[,"t"]) sel.template<-c(sel.template,which(abs(template[,"t"]-t)==min(abs(template[,"t"]-t)))[1])
	} else {
		filter.vector<-rep(1/input.filter,input.filter)
		filtered.track<-cbind(track[,"t"],filter(track[,"x"],filter.vector,method=filter.method),filter(track[,"y"],filter.vector,method=filter.method),filter(track[,"FI"],filter.vector,method=filter.method))
		filtered.track<-filtered.track[!is.na(filtered.track[,2]),]

		sel.template<-seq(which(abs(template[,"t"]-filtered.track[1,1])==min(abs(template[,"t"]-filtered.track[1,1]),na.rm=TRUE)),which(abs(template[,"t"]-filtered.track[length(filtered.track[,1]),1])==min(abs(template[,"t"]-filtered.track[length(filtered.track[,1]),1]),na.rm=TRUE)),by=1)
		splined.track.t<-seq(template[,"t"][sel.template[1]],template[,"t"][sel.template[length(sel.template)]],deltat)
		splined.track.x<-spline(filtered.track[,1],filtered.track[,2],xout=splined.track.t)
		splined.track.y<-spline(filtered.track[,1],filtered.track[,3],xout=splined.track.t)
		splined.track.FI<-spline(filtered.track[,1],filtered.track[,4],xout=splined.track.t)
	
		#substitute the track with the filtered and splined one
		track<-cbind(splined.track.x$x,splined.track.x$y,splined.track.y$y,splined.track.FI$y)
		dimnames(track)<-list(c(),c("t","x","y","FI"))
	}
	#define the weight if doweight=TRUE. the weight uses the FIs from the template and the track (cross correlating functions maximize FI product) and the standart deviation from the template
	if (doweight) weight<-as.vector(template[sel.template,"mean_sFI"])*as.vector(track[,"FI"])/(as.vector(template[sel.template,"std_x"])*as.vector(template[sel.template,"std_y"]))
	else weight<-rep(1,length(sel.template))


	#nomenclature: template is right coordinate while track is left coordinate, Horn's (1987) nomencalture
	#right coordinates center of mass
	rc<-c(wmean(template[sel.template,"mean_x"],weight),wmean(template[sel.template,"mean_y"],weight))	
	#left coordinates center of mass
	lc<-c(wmean(track[,"x"],weight),wmean(track[,"y"],weight))	
	r<-t(t(template[sel.template,c("mean_x","mean_y")])-rc)
	l<-t(t(track[,c("x","y")])-lc)
	
	Sxx<-wmean(l[,1]*r[,1],weight) 
	Sxy<-wmean(l[,1]*r[,2],weight) 
	Syx<-wmean(l[,2]*r[,1],weight) 
	Syy<-wmean(l[,2]*r[,2],weight) 
	
	A<-(Syx-Sxy)
	B<-(Sxx+Syy)
	M<--A/B
	theta1<-atan(M)
	theta2<-theta1+pi
	if (B*cos(theta1) >= A*sin(theta1)) {
		theta<-theta1
	} else {
		theta<-theta2
	}
	
	R<-Rmatrix(theta)
	trasl<-rc-R%*%lc
	return(list(R=R,trasl=trasl,cm=c(lc,rc)))
}
	

theta<-function(R){
	a<-c(1,0)
	b<-R%*%a
	atan2(a[1]*b[2]-a[2]*b[1],a[1]*b[1]+a[2]*b[2])
}

ctrl.plot.transf<-function(t,sel,outputlabel){
		#interpolation
		theta.fit<-lm(W2.theta[sel]~W1.theta[sel],data=t)
		lag.fit<-lm(W2.lag[sel]~W1.lag[sel],data=t)
		trasl.x.fit<-lm(W2.rtrasl.x[sel]~W1.trasl.x[sel],data=t)
		trasl.y.fit<-lm(W2.rtrasl.y[sel]~W1.trasl.y[sel],data=t)
		
		mythetas<-t[sel,"W1.theta"]-t[sel,"W2.theta"]
		mylags<-t[sel,"W1.lag"]-t[sel,"W2.lag"]
		mytrasls.x<-t[sel,"W1.trasl.x"]-t[sel,"W2.rtrasl.x"]
		mytrasls.y<-t[sel,"W1.trasl.y"]-t[sel,"W2.rtrasl.y"]

		pdf(paste(outputlabel,"fits.pdf",sep='_'))
		par(mfcol=c(2,2))
		plot(t[,"W1.theta"],t[,"W2.theta"],xlab=expression(theta[W1]),ylab=expression(theta[W2]))
		points(t[sel,"W1.theta"],t[sel,"W2.theta"],pch=16)
		lines(t[sel,"W1.theta"],predict(theta.fit),lty=2,col="red")
		plot(t[,"W1.lag"],t[,"W2.lag"],xlab=expression(Delta*t[W1]),ylab=expression(Delta*t[W2]))
		points(t[sel,"W1.lag"],t[sel,"W2.lag"],pch=16)
		lines(t[sel,"W1.lag"],predict(lag.fit),lty=2,col="red")
		plot(t[,"W1.trasl.x"],t[,"W2.rtrasl.x"],xlab=expression(Delta*x[W1]),ylab=expression(Delta*x[W2]))
		points(t[sel,"W1.trasl.x"],t[sel,"W2.rtrasl.x"],pch=16)
		lines(t[sel,"W1.trasl.x"],predict(trasl.x.fit),lty=2,col="red")
		plot(t[,"W1.trasl.y"],t[,"W2.rtrasl.y"],xlab=expression(Delta*y[W1]),ylab=expression(Delta*y[W2]))
		points(t[sel,"W1.trasl.y"],t[sel,"W2.rtrasl.y"],pch=16)
		lines(t[sel,"W1.trasl.y"],predict(trasl.y.fit),lty=2,col="red")

		#hists
		hist(mythetas,xlab=expression(theta))
		hist(mylags,xlab=expression(Delta*t))
		hist(mytrasls.x,xlab=expression(Delta*x))
		hist(mytrasls.y,xlab=expression(Delta*y))
		dev.off()	
		return(list(theta=theta.fit,lag=lag.fit,trasl.x=trasl.x.fit,trasl.y=trasl.y.fit))
}

ctrl.plot.align<-function(W1,W2,r,t,W1.transf,W2.transf){
	W1.rot<-rototrasl(W1,W1.transf$R,W1.transf$trasl)
#	plot(as.vector(W1[,"x"]),as.vector(W1[,"y"]),type="l",col="green",xlim=c(min(c(W1[,"x"],W1.rot[,"x"],r[,"mean_x"]),na.rm=TRUE),max(c(W1[,"x"],W1.rot[,"x"],r[,"mean_x"]),na.rm=TRUE)),ylim=c(min(c(W1[,"y"],W1.rot[,"y"],r[,"mean_y"]),na.rm=TRUE),max(c(W1[,"y"],W1.rot[,"y"],r[,"mean_y"]),na.rm=TRUE)),xlab="[pixels]",ylab="[pixels]")
#	lines(as.vector(r[,"mean_x"]),as.vector(r[,"mean_y"]),col="black")
#	lines(as.vector(W1.rot[,"x"]),as.vector(W1.rot[,"y"]),col="red")
#	plot.new()
	plot(as.vector(time(W1))[!is.na(W1[,"x"])],as.vector(W1[!is.na(W1[,"x"]),"FI"]),type="l",col="red",xlim=c(min(c(time(W1),time(r)),na.rm=TRUE),max(c(time(W1),time(r)),na.rm=TRUE)),ylim=c(min(c(W1[,3],r[,"mean_sFI"]),na.rm=TRUE),max(c(W1[,3],r[,"mean_sFI"]),na.rm=TRUE)),xlab="[s]",ylab="[pixels]",lwd=2)
	lines(as.vector(time(r)),as.vector(r[,"mean_sFI"]),col="black")
	plot(as.vector(time(W1.rot))[!is.na(W1[,"x"])],as.vector(W1.rot[!is.na(W1[,"x"]),"x"]),type="l",col="red",xlim=c(min(c(time(W1.rot),time(r)),na.rm=TRUE),max(c(time(W1.rot),time(r)),na.rm=TRUE)),ylim=c(min(c(W1.rot[,"x"],r[,"mean_x"]),na.rm=TRUE),max(c(W1.rot[,"x"],r[,"mean_x"]),na.rm=TRUE)),xlab="[s]",ylab="[pixels]",lwd=2)
	lines(as.vector(time(r)),as.vector(r[,"mean_x"]),col="black")
	W2.rot<-rototrasl(W2,W2.transf$R,W2.transf$trasl)
#	plot(as.vector(W2[,"x"]),as.vector(W2[,"y"]),type="l",col="green",xlim=c(min(c(W2[,"x"],W2.rot[,"x"],t[,"mean_x"]),na.rm=TRUE),max(c(W2[,"x"],W2.rot[,"x"],t[,"mean_x"]),na.rm=TRUE)),ylim=c(min(c(W2[,"y"],W2.rot[,"y"],t[,"mean_y"]),na.rm=TRUE),max(c(W2[,"y"],W2.rot[,"y"],t[,"mean_y"]),na.rm=TRUE)),xlab="[pixels]",ylab="[pixels]")
#	lines(as.vector(t[,"mean_x"]),as.vector(t[,"mean_y"]),col="black")
#	lines(as.vector(W2.rot[,"x"]),as.vector(W2.rot[,"y"]),col="red")
#	plot.new()
	if (is.null(t)) t<-r
	plot(as.vector(time(W2))[!is.na(W2[,"x"])],as.vector(W2[!is.na(W2[,"x"]),"FI"]),type="l",col="red",xlim=c(min(c(time(W2),time(t)),na.rm=TRUE),max(c(time(W2),time(t)),na.rm=TRUE)),ylim=c(min(c(W2[,3],t[,"mean_sFI"]),na.rm=TRUE),max(c(W2[,3],t[,"mean_sFI"]),na.rm=TRUE)),xlab="[pixels]",ylab="[pixels]",lwd=2)
	lines(as.vector(time(t)),as.vector(t[,"mean_sFI"]),col="black")
	plot(as.vector(time(W2.rot))[!is.na(W2[,"x"])],as.vector(W2.rot[!is.na(W2[,"x"]),"x"]),type="l",col="red",xlim=c(min(c(time(W2.rot),time(r)),na.rm=TRUE),max(c(time(W2.rot),time(r)),na.rm=TRUE)),ylim=c(min(c(W2.rot[,"x"],t[,"mean_x"]),na.rm=TRUE),max(c(W2.rot[,"x"],t[,"mean_x"]),na.rm=TRUE)),xlab="[pixels]",ylab="[pixels]",lwd=2)
	lines(as.vector(time(t)),as.vector(t[,"mean_x"]),col="black")
}

scale_FI<-function(ts1,ts2,which_FI=3,which_err=c(),which_FI2=c(),label="scaled_FI",label.err="scaled_err"){
	if (is.null(which_FI2)) which_FI2<-which_FI
	#load dimnames of the track which will be scaled. dimnames will be restored at the end of the process
	input_names<-dimnames(ts2)
	
	w.1<-ts1[,which_FI]/sum(ts1[,which_FI],na.rm=TRUE)
	w.2<-ts2[,which_FI2]/sum(ts2[,which_FI2],na.rm=TRUE)
	
	s<-wmean(x=ts1[,which_FI],w=w.1)/wmean(x=ts2[,which_FI2],w=w.2)
	
	ts2.6<-s*ts2[,which_FI2]
	if (!is.null(which_err)) {
		ts2<-cbind(ts2,ts2.6,s*ts2[,which_err])
		dimnames(ts2)<-list(input_names[[1]],c(input_names[[2]],label,label.err))
	} else {
		ts2<-cbind(ts2,ts2.6)
		if (is.null(input_names[[2]])){
			input_names<-dimnames(ts2)
			input_names[[2]]<-input_names[[2]][-length(input_names[[2]])]
		}
		dimnames(ts2)<-list(input_names[[1]],c(input_names[[2]],label))#substitute "ts2.6" dimname with "scaled_FI"
	}
	ts2
}
	
#main function
perform_alignment<-function(path.ref,path.target,inputfolder,dt,outputlabel,file_patternW1="W1data$",file_patternW2="W2data$",W1_fimax=FALSE,which_W1FI="mean_sFI",W2_fimax=FALSE,which_W2FI=3,filter.cc=5,cm=TRUE,max.angle=pi/3){
	#load SC data
	ref<-load.data(path.ref)
	ref<-center_mass(ref)
	if (!is.null(path.target)) {
		target<-load.data(path.target)
		#the target protein has to be centered on the center of mass, this is to minimize the dicrepancy between the rotation of the W2 traslation, done with the angle of the transformation for "that" DW pair
		#and the final transformation whose rotation is derived from the average of all angles:
		#	a' = R a + < t1 - R't2 > = R a0 + < t1 - R' t2> + R a_cm
		#	a' = R a0 + < t1 - R't2_0 >. 
		#	t2_0 = t2 - a_cm
		#however
		#	< t1 - R't2_0 > != <t1 - R' t2> + R a_cm
		#unless a_cm = 0
		target<-center_mass(target)
		transfs<-c()
	} else rotatedW2<-c() 
	if (length(inputfolder)!=length(dt)) stop("you need a dt for each input folder")
	all_filesW1<-list()
	all_filesW2<-list()
	for (i in 1:length(inputfolder)){ 
		all_filesW1[[i]]<-dir(inputfolder[i],pattern=file_patternW1,full.names=TRUE);
		all_filesW2[[i]]<-dir(inputfolder[i],pattern=file_patternW2,full.names=TRUE);
	}
	#variables for DW aligmnent without SC track	
	is_the_first_w2track=TRUE
	startend<-c()
	
	#load DW data and align them to ref (and target if present) to get the all transformation estimates
	pdf(paste(outputlabel,"alignments.pdf",sep="_"))
	par(mfcol=c(2,2))
	for (f in 1:length(inputfolder)){	
		for ( i in 1:length(all_filesW1[[f]])){
			print(all_filesW1[[f]][i])
			
			W1.track<-importts(all_filesW1[[f]][i],dt[f])
			#cross correlation with the reference
			W1.shift<-cc(W1.track,ref,filter.span=filter.cc)
			W1.track<-ts(W1.track,start=tsp(W1.track)[1]+W1.shift$lag,deltat=deltat(W1.track))
			if(W1_fimax) W1.transf<-trackalign(W1.track,window(ref,start=tsp(ref)[1],end=time(ref)[which(ref[,which_W1FI]==max(ref[,which_W1FI],na.rm=TRUE))]),doweight=TRUE)
			else W1.transf<-trackalign(W1.track,ref,doweight=TRUE)
			
			W2.track<-importts(all_filesW2[[f]][i],dt[f])
			if (!is.null(path.target)){
				W2.shift<-cc(W2.track,target,filter.span=filter.cc)
				W2.track<-ts(W2.track,start=tsp(W2.track)[1]+W2.shift$lag,deltat=deltat(W2.track))
				if(W2_fimax) W2.transf<-trackalign(W2.track,window(target,start=tsp(target)[1],end=time(target)[which(target[,which_W2FI]==max(target[,which_W2FI],na.rm=TRUE))]),doweight=TRUE)
				else W2.transf<-trackalign(W2.track,target,doweight=TRUE)
				#compute the angles and reduce them if they are multiple of 2*pi	
				W1.theta<-theta(W1.transf$R)
				W2.theta<-theta(W2.transf$R)
				while(abs(W1.theta-W2.theta)>pi) W1.theta<-W1.theta-sign(W1.theta-W2.theta)*2*pi
				#comment this.A
				#a and b are the ref/target, x and y the dw tracks W2 and W1, we want to firn a'' which match b
				#a = R1*x + t1
				#b = R2*y + t2
				#a' = R1^-1*a - R1^-1*t1
				#a'' = R2*R1^-1*a - R2*R1^-1*t1 + t2
				inverseR<-W1.transf$R%*%solve(W2.transf$R)	
				W2.rtrasl<-inverseR%*%W2.transf$trasl#see comment about center_mass of the target
				#collect the transfs
				transfs<-rbind(transfs,c(W1.shift$lag,W1.theta,W1.transf$trasl[1],W1.transf$trasl[2],W2.shift$lag,W2.theta,W2.transf$trasl[1],W2.transf$trasl[2],W2.rtrasl[1],W2.rtrasl[2]))
				ctrl.plot.align(W1.track,W2.track,ref,target,W1.transf,W2.transf)
			} else {
				if(is_the_first_w2track){
					first.W2.track<-W2.track
					is_the_first_w2track=FALSE
				}
				W2.track<-ts(W2.track,start=tsp(W2.track)[1]+W1.shift$lag,deltat=deltat(W2.track))
				startend<-rbind(startend,tsp(W2.track)[c(1,2)])
				W2.track<-scale_FI(first.W2.track,W2.track,which_FI="FI",label="sFI")
				#generate transformation to rotate W2
				rotatedW2<-rbind(rotatedW2,cbind(time(W2.track),rototrasl(W2.track,R=W1.transf$R,trasl=W1.transf$trasl)))
				ctrl.plot.align(W1.track,W2.track,ref,t=c(),W1.transf,W1.transf)
			}
		}
	}
	if (is.null(path.target)) dimnames(rotatedW2)<-list(c(),c("t",dimnames(W2.track)[[2]]))
	dev.off()
	#combine all the transformations together
	if (!is.null(path.target)) {
		dimnames(transfs)<-list(c(),c("W1.lag","W1.theta","W1.trasl.x","W1.trasl.y","W2.lag","W2.theta","W2.trasl.x","W2.trasl.y","W2.rtrasl.x","W2.rtrasl.y"))
		transfs<-data.frame(transfs)

		#estimate the overall theta and reject outlier (abs(theta)>pi/3)
		a<-sin(transfs[,"W1.theta"])*cos(transfs[,"W2.theta"])-cos(transfs[,"W1.theta"])*sin(transfs[,"W2.theta"])
		b<-cos(transfs[,"W1.theta"])*cos(transfs[,"W2.theta"])+sin(transfs[,"W1.theta"])*sin(transfs[,"W2.theta"])
		thetas<-atan2(a,b)
		sel<-which(abs(thetas)<max.angle)
		
		mythetas<-transfs[sel,"W1.theta"]-transfs[sel,"W2.theta"]
		mylags<-transfs[sel,"W1.lag"]-transfs[sel,"W2.lag"]
		mytrasls.x<-transfs[sel,"W1.trasl.x"]-transfs[sel,"W2.rtrasl.x"]
		mytrasls.y<-transfs[sel,"W1.trasl.y"]-transfs[sel,"W2.rtrasl.y"]
		
		#output plots named outputlabel_fits.pdf, with the quality control for the transformation and the fitting	
		fit<-ctrl.plot.transf(transfs,sel,outputlabel)
	
		angle<-as.numeric(mean(mythetas))
		mylag<-as.numeric(mean(mylags))
		trasl<-c(as.numeric(mean(mytrasls.x)),as.numeric(mean(mytrasls.y)))

		angle.err<-sd(mythetas)/sqrt(length(mythetas))
		mylags.err<-sd(mylags)/sqrt(length(mylags))
		trasl.err<-c(sd(mytrasls.x)/sqrt(length(mytrasls.x)),sd(mytrasls.y)/sqrt(length(mytrasls.y)))

		mytrasf<-cbind(c(angle,mylag,trasl),c(angle.err,mylags.err,trasl.err))
		dimnames(mytrasf)<-list(c("angle","lag","trasl1","trasl2"),c("trasf","trasf.err"))
						
		print(sprintf("angle: %0.2fpi lag: %0.2f trasl: (%0.2f,%0.2f)",angle/pi,mylag,trasl[1],trasl[2]))
		print("fimax:")
		print(W1_fimax)
		print(W2_fimax)
		print("sel:--------------------------------------------------------")
		print("")
		print(sel)
		print("lag:--------------------------------------------------------")
		print(summary(fit[["lag"]]))
		print("theta:------------------------------------------------------")
		print(summary(fit[["theta"]]))
		print("trasl1:-----------------------------------------------------")
		print(summary(fit[["trasl.x"]]))
		print("trasl2:-----------------------------------------------------")
		print(summary(fit[["trasl.y"]]))
		print("transf:-----------------------------------------------------")
		print(mytrasf)
	
		save(mytrasf,file=paste(outputlabel,"_mytrasf.Rdata",sep=""))
	} else {
		
		rotatedW2<-rotatedW2[order(rotatedW2[,1]),]
		time.span<-seq(min(rotatedW2[,1]),max(rotatedW2[,1]),by=min(dt)) #taking the min dt available, maximise the estimate of the standard error (i.e. the biggest possible), as less tracks might be bracketed in the time step
		averageW2<-c()
		sel.time.old<-1
		for (i in 2:length(time.span)){
			sel.time.span<-which(rotatedW2[,1]>=time.span[i-1] &rotatedW2[,1]<time.span[i])

			for (j in seq(sel.time.span[1]-1,sel.time.old[length(sel.time.old)])) {
				if (isTRUE(all.equal(rotatedW2[j,1],time.span[i-1]))) sel.time.span<-c(j,sel.time.span)
			}

			#while (isTRUE(all.equal(rotatedW2[[sel.time.span[1]-1,1]],time.span[1]))) sel.time.span<-c(sel.time.span[1]-1,sel.time.span)
			#while (i == length(time.span) & isTRUE(all.equal(rotatedW2[[sel.time.span[length(sel.time.span)]+1,1]],time.span[i]))) sel.time.span<-c(sel.time.span,sel.time.span[length(sel.time.span)]+1)
			#there might be some rounding errors in determining whetere a timepoint is >= than time.span[i-1] so let's double check them:	
			if (length(sel.time.span)>1) averageW2<-rbind(averageW2,c(apply(rotatedW2[sel.time.span,],2,mean,na.rm=TRUE),apply(rotatedW2[sel.time.span,],2,sd,na.rm=TRUE)/sqrt(length(sel.time.span)),length(sel.time.span)))
			else averageW2<-rbind(averageW2,c(rotatedW2[sel.time.span,],rep(NA,length(rotatedW2[sel.time.span,])),length(sel.time.span)))
			sel.time.old<-sel.time.span
		}
		averageW2<-averageW2[which(averageW2[,"t"]>mean(startend[,1]) & averageW2[,"t"]<mean(startend[,2])),]
		dimnames(averageW2)<-list(c(),c(paste("mean",dimnames(rotatedW2)[[2]],sep="_"),paste("std",dimnames(rotatedW2)[[2]],sep="_"),"n"))
		save(averageW2,file=paste(outputlabel,"_W2.Rdata",sep=""))

	}
}

gen.data<-function(data,data_trasf,n,sn,rescale.n=1,zero=0,t0=0,propagate.err=TRUE,pixel_size=100){
	#rescale the number of molecules according to Salmon correction
	n<-n*rescale.n
	sn<-sn*rescale.n
	
	if (!is.null(data_trasf)){
		d<-window(data$ts,start=data$start["mean"],end=data$end["mean"])
		d<-center_mass(d)
		if (propagate.err) d<-rototrasl(d,Rmatrix(data_trasf["angle","trasf"]),c(data_trasf["trasl1","trasf"],data_trasf["trasl2","trasf"]),transf.err=data_trasf[,"trasf.err"])
		else d<-rototrasl(d,Rmatrix(data_trasf["angle","trasf"]),c(data_trasf["trasl1","trasf"],data_trasf["trasl2","trasf"]))
	} else {
		d<-data
	}

	#VECTORS
	x<-as.vector(d[,"mean_x"])*pixel_size-zero
	y<-as.vector(d[,"mean_y"])*pixel_size
	sx<-as.vector(d[,"std_x"])*pixel_size
	sy<-as.vector(d[,"std_y"])*pixel_size
	if (!is.null(data_trasf)){
		t<-as.vector(time(d))-t0+data_trasf["lag","trasf"]
		st<-rep(data_trasf["lag","trasf.err"],length(t))
	} else {
		t<-as.vector(d[,"mean_t"])-t0
		st<-as.vector(d[,"std_t"])
	}
	fi<-as.vector(d[,"mean_sFI"])
	sfi<-as.vector(d[,"std_sFI"])
	N<-as.vector(d[,"n"])
	
	#number of molecules
	m<-mean(fi-min(fi,na.rm=TRUE),na.rm=TRUE)
	sm<-sd(fi-min(fi,na.rm=TRUE),na.rm=TRUE)/sqrt(length(fi[!is.na(fi)]))#se of the mean
	sf<-sfi/sqrt(as.vector(d[,"n"]))#se of fi
	dn<-n*(fi-min(fi,na.rm=TRUE))/m#number of molecules
	dsen<-sqrt((n*sf/m)^2+(n*(fi-min(fi,na.rm=TRUE))*sm/(m^2))^2+((fi-min(fi,na.rm=TRUE))*sn/m)^2)

	output<-cbind(t,st,x,y,sx,sy,fi,sf,dn,dsen,N)
	dimnames(output)<-list(c(),c("t","t.err","x","y","x.err","y.err","fi","fi.err","n","n.err","N"))
	return(output)
}	

#fancy ploting functio with shades for the errors		
linearize_errors<-function(e1,e2,e3,minmax="min",deltat=0.15){
	a<-rbind(e1,e2,e3)
	a.order<-order(a[,1])
	a.ordered<-a[a.order,]
	
	t<-seq(min(a.ordered[,1],na.rm=TRUE),max(a.ordered[,1],na.rm=TRUE),deltat)
	output<-c()	
	for (i in 1:(length(t)-1)){
		sel<-which(a.ordered[,1]>=t[i] & a.ordered[,1]<=t[i+1])
		if (minmax=="min" & length(sel)) output<-rbind(output,c((t[i]+t[i+1])/2,min(a.ordered[sel,2],na.rm=TRUE)))
		if (minmax=="max" & length(sel)) output<-rbind(output,c((t[i]+t[i+1])/2,max(a.ordered[sel,2],na.rm=TRUE)))
	}
	dimnames(output)<-list(c(),c("t","x"))
	return(output)	
}
#myplot<-function(t,x,e.t,e.x,line.col="#000000",err.col=c(),err.col.percent="50",line.lwd=3,hold_on=FALSE,deltat=0.15,...){
#	if (is.null(deltat)) deltat=min(t[-1]-t[-length(t)],na.rm=TRUE)
#	if(is.null(err.col)) err.col=paste(line.col,err.col.percent,sep="")
#	if (!hold_on) plot(t,x,col=line.col,lwd=line.lwd,type="l",...)
#	e.x.low<-x-e.x
#	e.x.high<-x+e.x
#	xlow<-cbind(t,e.x.low)
#	xhigh<-cbind(t,e.x.high)
#	
#	e.t.left<-t-e.t
#	tleft<-cbind(e.t.left,x)
#	e.t.right<-t+e.t
#	tright<-cbind(e.t.right,x)
#	low<-linearize_errors(xlow,tleft,tright,minmax="min",deltat=deltat)
#	high<-linearize_errors(xhigh,tleft,tright,minmax="max",deltat=deltat)
#	
#	polygon(c(low[,"t"],rev(high[,"t"])),c(low[,"x"],rev(high[,"x"])),col=err.col,border=NA)
#
#	
#	if (hold_on) lines(t,x,col=line.col,lwd=line.lwd,type="l")
#}
myplot<-function(t,x,e.t,e.x,line.col="#000000",input.filter=2,err.col=c(),err.col.percent="50",line.lwd=3,hold_on=FALSE,max.counts=15,deltat=c(),rotate90=FALSE,normalize=FALSE,...){

	if(is.null(deltat)) deltat<-min(t[-1]-t[-length(t)],na.rm=T)
	low<-c()
	high<-c()
	#flip x and y coords to compute the polygon, alg do not work then polygon is vertical, therefore rotation->polygon->rotation back
	if (rotate90){
		t.orig<-t
		x.orig<-x
		e.t.orig<-e.t
		e.x.orig<-e.x
		t<-x.orig
		x<-t.orig
		e.t<-e.x.orig
		e.x<-e.t.orig
	}
	if (normalize){
		x<-(x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
		e.x<-e.x/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
	}
	if(is.null(err.col)) err.col=paste(line.col,err.col.percent,sep="")
	terr<-cbind(x,t-e.t,t+e.t)
	xerr<-cbind(t,x-e.x,x+e.x)

	e1<-cbind(t,x-e.x)
	e2<-cbind(t,x+e.x)
	e3<-cbind(t-e.t,x)
	e4<-cbind(t+e.t,x)
	e.all<-rbind(e1,e2,e3,e4)

	e.all.order<-order(e.all[,1])
	e.all[e.all.order,]->e.all
	
	e.times<-seq(e.all[1,1],e.all[length(e.all[,1]),1],by=deltat)
	e.times<-c(e.times[1],e.times+deltat/2)
	for (j in 1:(length(e.times))){
		sel.e<-which(e.all[,1]>=e.times[j] & e.all[,1]<e.times[j+1])
		if(length(sel.e)){
			low<-rbind(low,c(e.all[sel.e,1][which(e.all[sel.e,2]==min(e.all[sel.e,2],na.rm=T))],min(e.all[sel.e,2],na.rm=T)))
			high<-rbind(high,c(e.all[sel.e,1][which(e.all[sel.e,2]==max(e.all[sel.e,2],na.rm=T))],max(e.all[sel.e,2],na.rm=T)))
		}
	}
	sel.e<-which(e.all[,1]>=e.times[j] & e.all[,1]<=e.all[length(e.all[,1]),1])
	if(length(sel.e)){
		low<-rbind(low,c(e.all[length(e.all[,1]),1],min(e.all[sel.e,2],na.rm=T)))
		high<-rbind(high,c(e.all[length(e.all[,1]),1],max(e.all[sel.e,2],na.rm=T)))
	}


	dimnames(low)<-list(c(),c("t","x"))
	dimnames(high)<-list(c(),c("t","x"))


	if(rotate90) {
		t<-t.orig
		x<-x.orig
	}
	if (!hold_on) {
		plot(t,x,col=line.col,lwd=line.lwd,type="l",...)
		if (rotate90) polygon(c(low[,"x"],rev(high[,"x"])),c(low[,"t"],rev(high[,"t"])),col=err.col,border=NA)
		else polygon(c(low[,"t"],rev(high[,"t"])),c(low[,"x"],rev(high[,"x"])),col=err.col,border=NA)
	} else {
		lines(t,x,col=line.col,lwd=line.lwd)
		if (rotate90) polygon(c(low[,"x"],rev(high[,"x"])),c(low[,"t"],rev(high[,"t"])),col=err.col,border=NA)
		else polygon(c(low[,"t"],rev(high[,"t"])),c(low[,"x"],rev(high[,"x"])),col=err.col,border=NA)
	}

	return(list(l=low,h=high))
}
#
