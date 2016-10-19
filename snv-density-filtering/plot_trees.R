library(ape)

# using phytools from https://github.com/liamrevell/phytools as of September 15, 2016
library(phytools)

plot_tree<-function(original_tree,tree,label) {
	plot(cophylo(original_tree,tree,rotate=TRUE),mar=c(2,1,2,1))
	title(main=label,adj=0.5)
	mtext("Original alignment",side=1, adj=0, line=0.3)
	mtext("SNVPhyl alignment",side=1, adj=1, line=0.3)
	box(which="plot", lty="solid", lwd="2", col="black")
}

reset <- function() {
	par(mfrow=c(1, 1), oma=c(0,0,0,0), mar=c(0,0,0,0), new=TRUE)
	plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

plot_all_trees<-function(original_tree,trees,labels) {
	n<-1
	for (i in 1:length(trees)) {
		if (i %% 4 == 1) {
			reset()
			dev.off()
			frame()
			plot.new()
			file_name<-paste("figure-S2-",n,".pdf",sep="")
			n<-n+1
			pdf(file_name,width=11,height=8.5)
			layout(t(matrix(1:4,2,2)))
			par(mar=c(2.5,0.5,2,0.5))
			par(oma=c(5,0,3,0))
		}

		plot_tree(original_tree,trees[[i]],labels[[i]])

		if (i %% 4 == 1) {
			mtext("Figure S2",line=1,outer=TRUE)
		}
	}
}

## MAIN ##

tree_dir<-"tree-distances"
original_tree<-read.tree(paste(tree_dir,"original_gubbins.phy_phyml_tree.txt",sep="/"))
original_tree<-root(original_tree,"reference",resolve.root=TRUE)
files<-c("snvphyl-no-filter.snvAlignment.phy_phyml_tree.txt", "snvphyl-2-20.snvAlignment.phy_phyml_tree.txt","snvphyl-2-100.snvAlignment.phy_phyml_tree.txt", "snvphyl-2-500.snvAlignment.phy_phyml_tree.txt", "snvphyl-2-1000.snvAlignment.phy_phyml_tree.txt","snvphyl-2-2000.snvAlignment.phy_phyml_tree.txt","snvphyl-gubbins.phy_phyml_tree.txt")
labels<-list("a) No filter", "b) 2 SNVs in 20 bp", "c) 2 SNVs in 100 bp", "d) 2 SNVs in 500 bp", "e) 2 SNVs in 1000 bp", "f) 2 SNVs in 2000 bp", "g) SNVPhyl then Gubbins")

trees<-list()

for(i in 1:length(files)) {
	tree<-read.tree(paste(tree_dir,files[i],sep="/"))
	tree<-root(tree,"reference",resolve.root=TRUE)
	trees[[length(trees)+1]]<-tree
}

plot_all_trees(original_tree,trees,labels)
