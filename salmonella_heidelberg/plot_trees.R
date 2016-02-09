library(ape)

max_snv_distance=5

get_outbreaks_for<-function(tree,table,col_id_name) {
	colname<-"Outbreak.number"
	values<-as.vector(table[match(tree$tip.label,table[[col_id_name]]),colname])
	values[is.na(values)]<-4
	values[values == "reference"]<-4
	return (as.numeric(values))
}

plot_tree<-function(tree,label,table,outbreaks,snv_matrix) {
	edgecolors<-rep("black",nrow(tree$edge))

	#colors<-c("blue","blue","blue","black")
	#tipcolors<-colors[get_outbreaks_for(tree,table,"Strain")]

	boxcolor<-"green"
	for(i in 1:length(outbreaks)) {
		if (!is_valid_cluster(tree,outbreaks[[i]],snv_matrix,max_snv_distance)) {
			edgecolors[which.edge(tree,outbreaks[[i]])]<-"red"
			boxcolor<-"red"
		} else {
			edgecolors[which.edge(tree,outbreaks[[i]])]<-"blue"
		}
	}

	#plot(tree,cex=0.5,edge.color=edgecolors,edge.width=3,tip.color=tipcolors,main=label,type="unrooted",show.tip.label=FALSE)
	plot(tree,cex=0.5,edge.color=edgecolors,edge.width=3,main=label,type="radial",show.tip.label=FALSE)

	nodelabels("O1",getMRCA(tree,outbreaks[[1]]),frame="circle",bg="white")
	nodelabels("O2",getMRCA(tree,outbreaks[[2]]),frame="circle",bg="white")
	nodelabels("O3",getMRCA(tree,outbreaks[[3]]),frame="circle",bg="white")

	box(which="plot", lty="solid", lwd="2", col=boxcolor)
}

plot_all_trees<-function(trees,labels,table,snv_matrices) {
	outbreak1<-as.vector(subset(table,Outbreak.number=="1")$Strain)
	outbreak2<-as.vector(subset(table,Outbreak.number=="2")$Strain)
	outbreak3<-as.vector(subset(table,Outbreak.number=="3")$Strain)
	outbreaks<-list(outbreak1,outbreak2,outbreak3)

	numtrees<-length(trees)
	size<-numtrees + (numtrees %% 2) # make size even
	plot.new()
	frame()
	layout(matrix(1:size,2,size/2))

	for (i in 1:length(trees)) {
		plot_tree(trees[[i]],labels[[i]],table,outbreaks,snv_matrices[[i]])
	}
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

get_invalid_isolates<-function(tree,isolates,snv_matrix,max_distance) {
	snv_matrix
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

experiments<-dir(c("experiments/alt","experiments/cov"), full.names=TRUE)
#experiments<-dir("experiments/alt", full.names=TRUE)
trees<-list()
mdists<-list()
cases<-list()

for(i in 1:length(experiments)) {
	tree_name<-list.files(experiments[i], pattern="pseudoalign.phy_phyml_tree.txt")
	mdist_name<-list.files(experiments[i], pattern="snp_matrix.tsv")
	case<-mdist_name

	tree<-read.tree(paste(experiments[i],tree_name,sep="/"))
	mdist<-read_snv_matrix(paste(experiments[i],mdist_name,sep="/"))

	trees[[length(trees)+1]]<-tree
	mdists[[length(mdists)+1]]<-mdist
	cases[[length(cases)+1]]<-case
}

plot_all_trees(trees,cases,strain_table,mdists)

#is.monophyletic
#nodepath
#which.edge
