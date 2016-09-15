library(ape)

plot_tree<-function(original_tree,tree,label) {
	association<-cbind(original_tree$tip.label,original_tree$tip.label)
	cophyloplot(original_tree,tree,assoc=association,gap=1,space=0,use.edge.length=TRUE,show.tip.label=TRUE)
	#plot(original_tree,cex=1.0,edge.color="black",edge.width=2,type="phylogram",show.tip.label=TRUE)
	#plot(tree,cex=1.0,edge.color="black",edge.width=2,type="phylogram",show.tip.label=TRUE,d="l")
	title(main=label,adj=0.5)
	box(which="plot", lty="solid", lwd="2", col="black")
}

reset <- function() {
	par(mfrow=c(1, 1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=TRUE)
	plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

plot_all_trees<-function(original_tree,trees,labels) {

	frame()
	plot.new()
	file_name<-"figure-S2.pdf"
	pdf(file_name,width=11,height=8.5)
	layout(t(matrix(1:4,2,2)))
	par(mar=c(2.5,0.5,2,0.5))
	par(oma=c(5,0,3,0))

	for (i in 1:length(trees)) {
		plot_tree(original_tree,trees[[i]],labels[[i]])
	}

	reset()

	dev.off()
}

## MAIN ##

tree_dir<-"tree-distances-with-gaps-n"
original_tree<-read.tree(paste(tree_dir,"original_gubbins.phy_phyml_tree.txt",sep="/"))
original_tree<-root(original_tree,"reference",resolve.root=TRUE)
files<-c("snvphyl-no-filter.snvAlignment.phy_phyml_tree.txt", "snvphyl-2-500.snvAlignment.phy_phyml_tree.txt", "snvphyl-gubbins.phy_phyml_tree.txt","original_gubbins.phy_phyml_tree.txt")
labels<-list("No filter", "2 SNVs in 500 bp","SNVPhyl then Gubbins","Original Gubbins")
#files<-list.files(tree_dir, pattern="tree.txt$")

trees<-list()

for(i in 1:length(files)) {
	tree<-read.tree(paste(tree_dir,files[i],sep="/"))
	tree<-root(tree,"reference",resolve.root=TRUE)
	tree<-rotateConstr(tree,original_tree$tip.label)
	trees[[length(trees)+1]]<-tree
}

plot_all_trees(original_tree,trees,labels)
