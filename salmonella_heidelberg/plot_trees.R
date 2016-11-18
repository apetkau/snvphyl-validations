library(ape)

max_snv_distance=5

normal_color<-"black"
outbreak_color<-"black"
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
	outbreak_colors=c(normal_color,normal_color,normal_color)
	for(i in 1:length(outbreaks)) {

		if (!is_valid_cluster(tree,outbreaks[[i]],snv_matrix,max_snv_distance)) {
			snv_distance<-max(snv_matrix[outbreaks[[i]],outbreaks[[i]]])
			if (!is.monophyletic(tree,outbreaks[[i]],reroot=FALSE)) {
				edgecolors[which.edge(tree,outbreaks[[i]])]<-failure_color
				boxcolor<-failure_color
				outbreak_colors[i]<-failure_color
				boxtype<-"dashed"
				failed_monophyletic<-TRUE
			}

			if (snv_distance >= max_snv_distance) {
				edgecolors[which.edge(tree,outbreaks[[i]])]<-failure_color
				boxcolor<-failure_color
				outbreak_colors[i]<-failure_color
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
			failure_label<-paste(failure_label,"Not monophyletic",sep='')
		}
		if (failed_snv_distance) {
			failure_label<-paste(failure_label,"Distance not within ",max_snv_distance, " SNVs", sep='')
		}
	}

	coresize_rounded=round(as.numeric(coresize['all','Percentage.of.all.positions.that.are.valid..included..and.part.of.the.core.genome']))

	plot(tree,cex=0.5,edge.color=edgecolors,edge.width=2,type="phylogram",show.tip.label=FALSE)
	title(main=label,adj=0.5)
	title(sub=paste(snvs_used, " SNVs",sep=''),adj=0,line=0.25)
	title(sub=paste(coresize_rounded,"% core",sep=''),adj=1,line=0.25)
	title(sub=failure_label,adj=0,line=1.25,col.sub=failure_color)

	nodelabels("1",getMRCA(tree,outbreaks[[1]]),frame="circle",bg="white",cex=1.2,font=2,col=outbreak_colors[1])
	nodelabels("2",getMRCA(tree,outbreaks[[2]]),frame="circle",bg="white",cex=1.2,font=2,col=outbreak_colors[2])
	nodelabels("3",getMRCA(tree,outbreaks[[3]]),frame="circle",bg="white",cex=1.2,font=2,col=outbreak_colors[3])

	box(which="plot", lty=boxtype, lwd="2", col=boxcolor)
}

plot_all_trees<-function(experiment,figure_num,figure_label,trees,labels,table,snv_matrices,coresizes,snvs_used_list) {
	outbreak1<-as.vector(subset(table,Outbreak.number=="1")$Strain)
	outbreak2<-as.vector(subset(table,Outbreak.number=="2")$Strain)
	outbreak3<-as.vector(subset(table,Outbreak.number=="3")$Strain)
	outbreaks<-list(outbreak1,outbreak2,outbreak3)

	for (i in 1:length(trees)) {
		plot_tree(trees[[i]],labels[[i]],table,outbreaks,snv_matrices[[i]],coresizes[[i]],snvs_used_list[[i]])
	}
	mtext(paste("Figure S2\n",figure_num,") ",figure_label,sep=''),line=0,outer=TRUE)
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
	return(is.monophyletic(tree,isolates,reroot=FALSE) && max(snv_matrix[isolates,isolates]) < max_distance)
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

experiment_names<-c("cov","scov","alt","contamination")
experiment_labels<-c("Minimum Coverage","Subsample coverage level","SNV Abundance Ratio","Contamination")
experiment_letters<-c("a","b","c","d")

pdf("Supplementary_Figure_S2.pdf",width=11,height=8.5)
layout(t(matrix(1:4,2,2)))
par(mar=c(2.5,0.5,2,0.5))
par(oma=c(5,0,3,0))

for (i in 1:length(experiment_names)) {
	experiment<-experiment_names[i]
	experiment_dirs<-sort(dir(paste("experiments",experiment,sep="/"), full.names=TRUE))
	
	trees<-list()
	mdists<-list()
	coresizes<-list()
	cases<-list()
	snvs_used_list<-list()
	
	for(j in 1:length(experiment_dirs)) {
		tree_name<-list.files(experiment_dirs[j], pattern="phylogeneticTree.newick")
		mdist_name<-list.files(experiment_dirs[j], pattern="snvMatrix.tsv")
		vcf2core_name<-list.files(experiment_dirs[j], pattern="vcf2core.tsv")
		filter_name<-list.files(experiment_dirs[j], pattern="filterStats.txt")
		title_name<-list.files(experiment_dirs[j], pattern="title")
	
		tree<-read.tree(paste(experiment_dirs[j],tree_name,sep="/"))
		mdist<-read_snv_matrix(paste(experiment_dirs[j],mdist_name,sep="/"))
		vcf2core<-read.delim(paste(experiment_dirs[j],vcf2core_name,sep="/"),row.names=1)
		case<-paste(readLines(paste(experiment_dirs[j],title_name,sep="/")))
	
		filter_stats<-paste(readLines(paste(experiment_dirs[j],filter_name,sep="/")))
		snvs_used<-sub("Number of sites used to generate phylogeny: ([0-9]+)$","\\1",grep("^Number of sites used to generate phylogeny: ",value=TRUE,filter_stats))
	
		trees[[length(trees)+1]]<-tree
		mdists[[length(mdists)+1]]<-mdist
		coresizes[[length(coresizes)+1]]<-vcf2core
		cases[[length(cases)+1]]<-case
		snvs_used_list[[length(snvs_used_list)+1]]<-snvs_used
	}
	
	plot_all_trees(experiment,experiment_letters[i],experiment_labels[i],trees,cases,strain_table,mdists,coresizes,snvs_used_list)
}
