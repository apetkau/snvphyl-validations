get_outbreaks_for<-function(tree,table,col_id_name) {
	colname<-"Outbreak"
	values<-as.vector(table[match(tree$tip.label,table[[col_id_name]]),colname])
	values[is.na(values)]<-4
	values[values == "Sporadic"]<-4
	values[values == "reference"]<-4
	return (as.numeric(values))
}

rename_tips<-function(tree,table) {
	sh_id<-table[match(tree$tip.label,table[["NLEP.ID"]]),"SH.ID"]
	tree$tip.label<-as.vector(sh_id)
	return(tree)
}

plot_tree<-function(tree,label,table) {
	plot(tree,cex=0.5,tip.color=colors[get_outbreaks_for(tree,table,"SH.ID")],main=label,type="unrooted")
	box(which="plot", lty="solid", lwd="2", col="green")
}

plot_all_trees<-function(trees,labels,table) {
	numtrees<-length(trees)
	size<-numtrees + (numtrees %% 2) # make size even
	plot.new()
	layout(matrix(1:size,2,size/2))

	for (i in 1:length(trees)) {
		plot_tree(trees[[i]],labels[[i]],table)
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

is_valid_cluster<-function(isolates,snv_matrix,max_distance) {
	# First, check if isolates are within max distance
	return(max(snv_matrix[isolates,isolates]) <= max_distance)
}

all_valid_clusters<-function(outbreaks,snv_matrix,max_distance) {
	is_valid<-TRUE
	for(i in 1:length(outbreaks)) {
		is_valid<-is_valid && is_valid_cluster(outbreaks[[i]],snv_matrix,max_distance)
	}

	return(is_valid)
}
