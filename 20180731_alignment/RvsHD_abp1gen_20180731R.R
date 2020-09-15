## Align Sla1/Abp1 (DV) to Sla1 (SC), to get Abp1-mC average

source("/Users/kaksonenlab/Desktop/20180731_alignment/align_tracks.R")

##
#path.abp1="/Users/deepikaa/Desktop/data_desktop/Tracking/SC/WT/Abp1/2104_03_29_abp1egfp/deepikaa_abp1egfp_1.Rdata"
path.abp1="/Users/kaksonenlab/Desktop/20180731_alignment/tracksAbp.Rdata"

#path.target="/Users/deepikaa/Desktop/data_desktop/Tracking/SC/WT/Inp52gfp/2017_03_13/tracks.Rdata"
path.target="/Users/kaksonenlab/Desktop/20180731_alignment/data2014_03_26_rvs.Rdata"

output.label="RvsHD_abpgen"

##folder where the dc tracks are that need to be aligned
#input.folder="/Users/deepikaa/Desktop/data_desktop/Tracking/DC/wt/inp52GFP_abpmCH/images/stacks/tracks/"
input.folder="/Users/kaksonenlab/Desktop/20180731_alignment/tracks1/"
dt=0.2715

perform_alignment(path.ref=path.abp1,path.target=path.target,inputfolder=input.folder,dt=dt,outputlabel=output.label,W1_fimax=FALSE,W2_fimax=FALSE,file_patternW1="W1.txt$",file_patternW2="W2.txt$")
