library(ape)

max_snv_distance=5

get_outbreaks_for<-function(tree,table,col_id_name) {
	colname<-"Outbreak.number"
	values<-as.vector(table[match(tree$tip.label,table[[col_id_name]]),colname])
	values[is.na(values)]<-4
	values[values == "reference"]<-4
	return (as.numeric(values))
}

plot_tree<-function(tree,label,table,outbreaks,snv_matrix,coresize,snvs_used) {
	edgecolors<-rep("black",nrow(tree$edge))

	#colors<-c("blue","blue","blue","black")
	#tipcolors<-colors[get_outbreaks_for(tree,table,"Strain")]

	boxcolor<-"black"
	failure_label<-""
	failure_snv_distance<-0
	failed_monophyletic<-FALSE
	failed_snv_distance<-FALSE
	for(i in 1:length(outbreaks)) {

		if (!is_valid_cluster(tree,outbreaks[[i]],snv_matrix,max_snv_distance)) {
			snv_distance<-max(snv_matrix[outbreaks[[i]],outbreaks[[i]]])
			if (!is.monophyletic(tree,outbreaks[[i]],reroot=FALSE)) {
				edgecolors[which.edge(tree,outbreaks[[i]])]<-"red"
				boxcolor<-"red"
				failed_monophyletic<-TRUE
			}

			if (snv_distance > max_snv_distance) {
				edgecolors[which.edge(tree,outbreaks[[i]])]<-"red"
				boxcolor<-"red"
				failed_snv_distance<-TRUE
				failure_snv_distance<-max(failure_snv_distance,snv_distance)
			}
		} else {
			edgecolors[which.edge(tree,outbreaks[[i]])]<-"blue"
		}
	}

	if (failed_monophyletic) {
		failure_label<-"Failed monophyletic\n"
	}
	if (failed_snv_distance) {
		failure_label<-paste(failure_label,"Failed Distance ",failure_snv_distance," > ",max_snv_distance)
	}

	#plot(tree,cex=0.5,edge.color=edgecolors,edge.width=3,tip.color=tipcolors,main=label,type="unrooted",show.tip.label=FALSE)
	plot(tree,cex=0.5,edge.color=edgecolors,edge.width=2,type="phylogram",show.tip.label=FALSE)
	title(main=label,sub=paste(coresize['all','Percentage.in.core'],"% core\nSNVs Used ",snvs_used,"\n",failure_label,sep=''))

	nodelabels("1",getMRCA(tree,outbreaks[[1]]),frame="circle",bg="white")
	nodelabels("2",getMRCA(tree,outbreaks[[2]]),frame="circle",bg="white")
	nodelabels("3",getMRCA(tree,outbreaks[[3]]),frame="circle",bg="white")

	box(which="plot", lty="solid", lwd="2", col=boxcolor)
}

plot_all_trees<-function(trees,labels,table,snv_matrices,coresizes,snvs_used_list) {
	outbreak1<-as.vector(subset(table,Outbreak.number=="1")$Strain)
	outbreak2<-as.vector(subset(table,Outbreak.number=="2")$Strain)
	outbreak3<-as.vector(subset(table,Outbreak.number=="3")$Strain)
	outbreaks<-list(outbreak1,outbreak2,outbreak3)

	numtrees<-length(trees)
	size<-numtrees + (numtrees %% 4) # make multiple of 4
	plot.new()
	frame()
	pdf("figure3_trees.pdf",width=11,height=8.5)
	layout(t(matrix(1:size,4,size/4)))

	for (i in 1:length(trees)) {
		plot_tree(trees[[i]],labels[[i]],table,outbreaks,snv_matrices[[i]],coresizes[[i]],snvs_used_list[[i]])
	}

	dev.off()
}

root_on_tip<-function(tree,tip_label) {
	m<-match(tip_label,tree$tip.label)
	return(root(tree,outgroup=m,resolve.root=TRUE))
}

read_snv_matrix<-function(file) {
	m<-read.delim(file,header=FALSE,row.names=c(1),skip=1)
	colnames(m)<-rownames(m)
	return(m)
}

is_valid_cluster<-function(tree,isolates,snv_matrix,max_distance) {
	return(is.monophyletic(tree,isolates,reroot=FALSE) && max(snv_matrix[isolates,isolates]) <= max_distance)
}

all_valid_clusters<-function(tree,outbreaks,snv_matrix,max_distance) {
	is_valid<-TRUE
	for(i in 1:length(outbreaks)) {
		is_valid<-is_valid && is_valid_cluster(tree,outbreaks[[i]],snv_matrix,max_distance)
	}

	return(is_valid)
}

## MAIN ##

strain_table<-read.delim("strain_table.txt")
outbreak1<-as.vector(subset(strain_table,Outbreak.number=="1")$Strain)
outbreak2<-as.vector(subset(strain_table,Outbreak.number=="2")$Strain)
outbreak3<-as.vector(subset(strain_table,Outbreak.number=="3")$Strain)

experiments_cov<-sort(dir(c("experiments/cov"), full.names=TRUE))
experiments_alt<-sort(dir(c("experiments/alt"), full.names=TRUE))
experiments_contamination<-sort(dir(c("experiments/contamination"), full.names=TRUE))
experiments=c(experiments_cov,experiments_alt,experiments_contamination)
#experiments<-dir("experiments/scov2", full.names=TRUE)
trees<-list()
mdists<-list()
coresizes<-list()
cases<-list()
snvs_used_list<-list()

for(i in 1:length(experiments)) {
	tree_name<-list.files(experiments[i], pattern="pseudoalign.phy_phyml_tree.txt")
	mdist_name<-list.files(experiments[i], pattern="snp_matrix.tsv")
	vcf2core_name<-list.files(experiments[i], pattern="vcf2core.csv")
	filter_name<-list.files(experiments[i], pattern="filter-stats.tsv")
	title_name<-list.files(experiments[i], pattern="title")

	tree<-read.tree(paste(experiments[i],tree_name,sep="/"))
	mdist<-read_snv_matrix(paste(experiments[i],mdist_name,sep="/"))
	vcf2core<-read.delim(paste(experiments[i],vcf2core_name,sep="/"),row.names=1)
	case<-paste(readLines(paste(experiments[i],title_name,sep="/")))

	filter_stats<-paste(readLines(paste(experiments[i],filter_name,sep="/")))
	snvs_used<-sub("Number of sites used to generate phylogeny: ([0-9]+)$","\\1",grep("^Number of sites used to generate phylogeny: ",value=TRUE,filter_stats))

	trees[[length(trees)+1]]<-tree
	mdists[[length(mdists)+1]]<-mdist
	coresizes[[length(coresizes)+1]]<-vcf2core
	cases[[length(cases)+1]]<-case
	snvs_used_list[[length(snvs_used_list)+1]]<-snvs_used
}

plot_all_trees(trees,cases,strain_table,mdists,coresizes,snvs_used_list)
