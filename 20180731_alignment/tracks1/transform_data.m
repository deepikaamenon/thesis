
load('/Users/deepikaa/Desktop/data_desktop/Tracking/DC/wt/Rvs167/2017_04_06_3131/tracks/transform.mat')
W1files=dir('*_W1.txt');
for i=1:length(W1files)
	sprintf(W1files(i).name)
	rawW1=load(W1files(i).name,'-ascii'); 
	
	newW1=tforminv(warp,rawW1(:,2:3));
	outputW1=[rawW1(:,1),newW1,rawW1(:,4:length(rawW1(1,:)))];
	save(W1files(i).name,'-ascii','outputW1')
	save(strcat(W1files(i).name(1:(length(W1files(i).name))),'_original'),'-ascii','rawW1')
end


