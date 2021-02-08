input_file = "FastQC_stats.txt"

input.table = read.table(input_file, head=T, sep="\t")

flag2num = function(flag_arr){
	num_arr=c()
	
	for(i in 1:length(flag_arr)){
		flag_value = flag_arr[i]
		
		if(is.na(flag_value)){
			num_arr[i]=NA
		}else if(flag_value == "PASS"){
			num_arr[i]=1
		}else{
			num_arr[i]=0
		}
	}#end for(i in 1:length(flag_arr))
	
	return(num_arr)
}#end def bool2num

status.mat = input.table[,3:ncol(input.table)]
status.mat = apply(status.mat, 2, flag2num)

R1_mat = status.mat[grep("_1$",input.table$Sample),]
R2_mat = status.mat[grep("_2$",input.table$Sample),]

png("FastQC_flag_summary.png", height=400, width=800)
par(mfcol=c(1,2))

R1_status = matrix(ncol=ncol(R1_mat),nrow=2)
rownames(R1_status)=c("PASS","FAIL")
colnames(R1_status)=colnames(R1_mat)
for(i in 1:ncol(R1_status)){
	temp_cat = R1_mat[,i]
	pass_percent = length(temp_cat[temp_cat==1]) / nrow(R1_mat)
	fail_percent = length(temp_cat[temp_cat==0]) / nrow(R1_mat)
	R1_status[1,i]=pass_percent
	R1_status[2,i]=fail_percent
}#end for(i in 1:ncol(R1_status))
par(mar=c(10,5,3,1))
barplot(R1_status, las=2, col=c("black","red"),
			ylab="Pass Frequency (Forward Read)", cex.names=0.9)
legend("top",legend=c("PASS","FAIL"),col=c("black","red"), pch=15,
		ncol=2, xpd=T, inset=-0.2)

R2_status = matrix(ncol=ncol(R2_mat),nrow=2)
rownames(R2_status)=c("PASS","FAIL")
colnames(R2_status)=colnames(R2_mat)
for(i in 1:ncol(R2_status)){
	temp_cat = R2_mat[,i]
	pass_percent = length(temp_cat[temp_cat==1]) / nrow(R2_mat)
	fail_percent = length(temp_cat[temp_cat==0]) / nrow(R2_mat)
	R2_status[1,i]=pass_percent
	R2_status[2,i]=fail_percent
}#end for(i in 1:ncol(R1_status))
barplot(R2_status, las=2, col=c("black","red"),
			ylab="Pass Frequency (Reverse Read)", cex.names=0.9)
legend("top",legend=c("PASS","FAIL"),col=c("black","red"), pch=15,
		ncol=2, xpd=T, inset=-0.2)
dev.off()