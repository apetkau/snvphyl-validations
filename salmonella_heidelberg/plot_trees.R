library(ape)

max_snv_distance=5

normal_color<-"black"
outbreak_color<-"blue"
failure_color<-"red"


get_outbreaks_for<-function(tree,table,col_id_name) {
	colname<-"Outbreak.number"
	values<-as.vector(table[match(tree$tip.label,table[[col_id_name]]),colname])
	values[is.na(values)]<-4
	values[values == "reference"]<-4
	return (as.numeric(values))
}

plot_tree<-function(tree,label,table,outbreaks,snv_matrix,coresize,snvs_used) {
	edgecolors<-rep("black",nrow(tree$edge))

	boxcolor<-normal_color
	boxtype<-"solid"
	failure_label<-""
	failure_snv_distance<-0
	failed_monophyletic<-FALSE
	failed_snv_distance<-FALSE
	for(i in 1:length(outbreaks)) {

		if (!is_valid_cluster(tree,outbreaks[[i]],snv_matrix,max_snv_distance)) {
			snv_distance<-max(snv_matrix[outbreaks[[i]],outbreaks[[i]]])
			if (!is.monophyletic(tree,outbreaks[[i]],reroot=FALSE)) {
				edgecolors[which.edge(tree,outbreaks[[i]])]<-failure_color
				boxcolor<-failure_color
				boxtype<-"dashed"
				failed_monophyletic<-TRUE
			}

			if (snv_distance > max_snv_distance) {
				edgecolors[which.edge(tree,outbreaks[[i]])]<-failure_color
				boxcolor<-failure_color
				boxtype<-"dashed"
				failed_snv_distance<-TRUE
				failure_snv_distance<-max(failure_snv_distance,snv_distance)
			}
		} else {
			edgecolors[which.edge(tree,outbreaks[[i]])]<-outbreak_color
		}
	}

	if (failed_monophyletic || failed_snv_distance) {
		failure_label<-"Failed: "
		if (failed_monophyletic) {
			failure_label<-paste(failure_label,"M",sep='')
		}
		if (failed_snv_distance) {
			failure_label<-paste(failure_label," D",failure_snv_distance,">",max_snv_distance, sep='')
		}
	}

	coresize_rounded=round(as.numeric(coresize['all','Percentage.of.all.positions.that.are.valid..included..and.part.of.the.core.genome']))

	plot(tree,cex=0.5,edge.color=edgecolors,edge.width=2,type="phylogram",show.tip.label=FALSE)
	title(main=label,adj=0.5)
	title(sub=paste(snvs_used, " SNVs",sep=''),adj=0,line=0.25)
	title(sub=paste(coresize_rounded,"% core",sep=''),adj=1,line=0.25)
	title(sub=failure_label,adj=0,line=1.25,col.sub=failure_color)

	nodelabels("1",getMRCA(tree,outbreaks[[1]]),frame="circle",bg="white",cex=1.2)
	nodelabels("2",getMRCA(tree,outbreaks[[2]]),frame="circle",bg="white",cex=1.2)
	nodelabels("3",getMRCA(tree,outbreaks[[3]]),frame="circle",bg="white",cex=1.2)

	box(which="plot", lty=boxtype, lwd="2", col=boxcolor)
}

reset <- function() {
	par(mfrow=c(1, 1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=TRUE)
	plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
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
	par(mar=c(2.5,0.5,2,0.5))
	par(oma=c(5,0,3,0))

	for (i in 1:length(trees)) {
		plot_tree(trees[[i]],labels[[i]],table,outbreaks,snv_matrices[[i]],coresizes[[i]],snvs_used_list[[i]])
	}
	mtext("Phylogenetic trees",line=1,outer=TRUE)

	reset()
	legend("bottom",horiz=TRUE,cex=0.75,legend=c("Normal","Outbreak (1,2,3)","Failure"),fill=c(normal_color,outbreak_color,failure_color),xpd=NA)

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
experiments_scov<-sort(dir(c("experiments/scov"), full.names=TRUE))
experiments_alt<-sort(dir(c("experiments/alt"), full.names=TRUE))
experiments_contamination<-sort(dir(c("experiments/contamination"), full.names=TRUE))
experiments=c(experiments_cov,experiments_scov,experiments_alt,experiments_contamination)
#experiments<-dir("experiments/scov2", full.names=TRUE)
trees<-list()
mdists<-list()
coresizes<-list()
cases<-list()
snvs_used_list<-list()

for(i in 1:length(experiments)) {
	tree_name<-list.files(experiments[i], pattern="phylogeneticTree.newick")
	mdist_name<-list.files(experiments[i], pattern="snvMatrix.tsv")
	vcf2core_name<-list.files(experiments[i], pattern="vcf2core.tsv")
	filter_name<-list.files(experiments[i], pattern="filterStats.txt")
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
