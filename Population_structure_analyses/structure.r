#!/usr/bin/env Rscript

get_genes<-function(gff){
    genes=read.table(gff,sep="\t",quote="")
    sub=genes$V3=="CDS"
    genes=genes[sub,]
    genes2=data.frame(Contig=genes$V1,Start=genes$V4,End=genes$V5)
    s1=genes$V9
    s2=unlist(regmatches(s1,gregexpr("ID=.*?;",s1))) #extract gene names
    genes2$Name=substr(s2,4,nchar(s2)-1)
    genes2
}

get_samples<-function(samplefile){

	data=read.table(samplefile)
	samples=as.character(data$V2)
	names(samples)=data$V1
	samples
}

get_snps<-function(samples){

	print("Reading data")
	alldata=list()
	for (s in names(samples)){
		print(s)
		alldata[[s]]=read.table(samples[s],header=T)
	}

	data=alldata[[1]]
	for (i in seq(2,length(samples))){
		data=cbind(data,alldata[[i]][,c("A","C","G","T")])
	}

	names=c("Contig","Pos","Ref")
	for (s in names(samples)){
		names=c(names,paste(s,c("A","C","G","T"),sep="."))
	}
	names(data)=names
	data
}

add_pi<-function(snps,samples){
	
	ns=dim(snps)[1]
	print ("Pi calculation")
	for (s in samples){
		print(s)
		snps[,paste("pi_",s,sep="")]=rep(0,ns)
		names=paste(s,c("A","C","G","T"),sep=".")
		for (i in seq(1,ns)){
			subd=snps[i,names]
			sd=rowSums(subd)
			dem=1/sd/(sd-1)
			dem[!is.finite(dem)]=0 #do not count the one with 1 or 0 entries
			pi=2*(subd[,1]*subd[,2]*dem+subd[,1]*subd[,3]*dem+subd[,1]*subd[,4]*dem+subd[,2]*subd[,3]*dem+subd[,2]*subd[,4]*dem+subd[,3]*subd[,4]*dem)
			snps[i,paste("pi_",s,sep="")]=pi
		}
	}
	snps
}

add_pi2<-function(snps,samples){

	print("Between-sample pi calculation")
	ns=dim(snps)[1]
	print(ns)
	for (i1 in seq(1,length(samples)-1)){
		s1=samples[i1]
		for (i2 in seq(i1+1,length(samples))){
			s2=samples[i2]
			print(c(s1,s2))
			name=paste("pi",s1,s2,sep="_")
			snps[,name]=rep(0,ns)
			names1=paste(s1,c("A","C","G","T"),sep=".")
			names2=paste(s2,c("A","C","G","T"),sep=".")
			for (i in seq(1,ns)){
				sub1=snps[i,names1]
				sub2=snps[i,names2]
				sum1=rowSums(sub1)
				sum2=rowSums(sub2)
				dem=1/sum1/sum2
				dem[!is.finite(dem)]=0
				pi=sub1[,1]*sub2[,2]*dem+sub1[,1]*sub2[,3]*dem+sub1[,1]*sub2[,4]*dem+sub1[,2]*sub2[,3]*dem+sub1[,2]*sub2[,4]*dem+sub1[,3]*sub2[,4]*dem+sub2[,1]*sub1[,2]*dem+sub2[,1]*sub1[,3]*dem+sub2[,1]*sub1[,4]*dem+sub2[,2]*sub1[,3]*dem+sub2[,2]*sub1[,4]*dem+sub2[,3]*sub1[,4]*dem
				snps[i,name]=pi
			}
		}
	}
	snps
}

add_fst<-function(snps,samples){
	
	print("Position Fst calculation")
	
	nsamp=length(samples)
	names1=paste("pi",samples,sep="_") #sample pi
	n=names(snps)
	n=n[substr(n,1,3)=="pi_"]
	names2=setdiff(n,names1)  #between sample pi
	
	pi1=snps[,names1]
	pi2=snps[,names2]
	pi1a=rowSums(pi1)/nsamp
	pi2a=rowSums(pi2)/(nsamp*(nsamp-1)/2)
	snps$Fst=1-(pi1a/pi2a)
	
	snps
}

print_fst_pos<-function(snps,samplefile){

	name=paste(samplefile,"Fst_pos.txt",sep="")
	print(paste("Writing",name))
	write.table(snps,name,sep="\t",quote=F,row.names=F)
}

add_pi_genes<-function(snps,genes,samples){

	print("Pi gene calculation")

	ng=dim(genes)[1]
	genes$SNPs=rep(0,ng)
	for (g in seq(1,ng)){
		start=genes[g,"Start"]
		end=genes[g,"End"]
		sub=snps$Contig==as.character(genes[g,"Contig"])&snps$Pos>=start&snps$Pos<=end
		genes[g,"SNPs"]=sum(sub)
	}

	for (s in samples){
		name=paste("pi_",s,sep="")
		genes[,name]=rep(0,ng)
	}

	for (g in seq(1,ng)){
		if(!genes[g,"SNPs"]){next;}
		start=genes[g,"Start"]
		end=genes[g,"End"]
		sub=snps$Contig==as.character(genes[g,"Contig"])&snps$Pos>=start&snps$Pos<=end
		
		for (s in samples){
			name=paste("pi_",s,sep="")
			pi=sum(snps[sub,name])
			pi=pi/(end-start+1)
			genes[g,name]=pi
		}
	}
	genes
}

add_pi_genes2<-function(snps,genes,samples){
	
	print("Between-sample pi gene calculation")

	ng=dim(genes)[1]
	for (i1 in seq(1,length(samples)-1)){
		s1=samples[i1]
		for (i2 in seq(i1+1,length(samples))){
			s2=samples[i2]
			print(c(s1,s2))
			name=paste("pi",s1,s2,sep="_")
			genes[,name]=rep(0,ng)
		}
	}

	for (g in seq(1,ng)){
		if(!genes[g,"SNPs"]){next;}
		start=genes[g,"Start"]
		end=genes[g,"End"]
		sub=snps$Contig==as.character(genes[g,"Contig"])&snps$Pos>=start&snps$Pos<=end
		
		for (i1 in seq(1,length(samples)-1)){
			s1=samples[i1]
			for (i2 in seq(i1+1,length(samples))){
				s2=samples[i2]
				name=paste("pi",s1,s2,sep="_")
				pi=sum(snps[sub,name])
				pi=pi/(end-start+1)
				genes[g,name]=pi
			}
		}
	}
	genes
}

add_fst_genes<-function(genes,samples){

	print("Gene Fst calculation")
	
	nsamp=length(samples)
	names1=paste("pi",samples,sep="_")
	n=names(genes)
	n=n[substr(n,1,3)=="pi_"]
	names2=setdiff(n,names1)  #between sample pi
	pi1=genes[,names1]
	pi2=genes[,names2]
	pi1a=rowSums(pi1)/nsamp
	pi2a=rowSums(pi2)/(nsamp*(nsamp-1)/2)
	genes$Fst=1-(pi1a/pi2a)

	for (i1 in seq(1,length(samples)-1)){
		s1=samples[i1]
		for (i2 in seq(i1+1,length(samples))){
			s2=samples[i2]
			name=paste("Fst",s1,s2,sep="_")
			genes[,name]=1-(genes[,paste("pi",s1,sep="_")]+genes[,paste("pi",s2,sep="_")])/2/genes[,paste("pi",s1,s2,sep="_")]
		}
	}
	genes
}

print_fst_genes<-function(genes,samplefile){

	name=paste(samplefile,"Fst.txt",sep="")
	print(paste("Writing",name))
	write.table(genes,name,sep="\t",quote=F,row.names=F)
}

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Usage: structure.r <samples> <gff>", call.=FALSE)
}

samplefile=args[1]
gff=args[2]

samples=get_samples(samplefile)
snps=get_snps(samples)
snps=add_pi(snps,names(samples))
snps=add_pi2(snps,names(samples))
snps=add_fst(snps,names(samples))
print_fst_pos(snps,samplefile)

genes=get_genes(gff)
genes=add_pi_genes(snps,genes,names(samples))
genes=add_pi_genes2(snps,genes,names(samples))
genes=add_fst_genes(genes,names(samples))
print_fst_genes(genes,samplefile)








