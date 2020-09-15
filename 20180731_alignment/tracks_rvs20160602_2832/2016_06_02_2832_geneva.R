#run R
source("/Volumes/MarkoKaksonenLab/Deepika/data_desktop/Doc/align_SC_tracks.R")

dir<-"tracks_rvs20160602_2832"
dt<-0.1035
movie_length=500
min_w=80
output_name<-"20180731_redoalign"

remove<-c(14,53)

all_align_all(dir=dir,min_w=min_w,Delta=dt,output_name=output_name,perform_alignment=FALSE,remove=remove,movie_length=movie_length,fimax=FALSE)
