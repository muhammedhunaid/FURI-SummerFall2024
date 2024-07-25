
# Load required library
library(ape)  ## generate neighbor joining tree

# Read the input data
ffc='/Users/MuhammedH/github/FURI-SummerFall2024/ClusteringBenchmarks/GIANA_To_GIANA/tmp_query.txt'  ## Adjust the path as needed
ffm='TCRantigen--RotationEncodingBL62.txt_EncodingMatrix.txt'  ## Placeholder for the matrix file path

d1=read.table(ffc, header=F, sep='\t', stringsAsFactors = F)

# Assuming ffm is available, if not, this part will be skipped or a placeholder matrix needs to be created.
# Mat=read.table(ffm, header=F, sep='\t', stringsAsFactors = F)

# Combining antigen annotation
tmp=strsplit(d1[,6],'_')
tmp=unlist(sapply(tmp,function(x)x[length(x)]))
d1=cbind(d1, Ag=tmp)
d1[grep('Antigen:Human gammaherpesvirus 4', d1[,10]), 10]='GLCTLVAML'
d1[grep('Antigen:EBV', d1[,10]), 10]='GLCTLVAML'
d1[grep('Antigen:CMV', d1[,10]), 10]='NLVPMVATV'
d1[grep('Antigen:HIV-1', d1[,10]), 10]='KAFSPEVIPMF'
d1[grep('Antigen:Influenza A virus', d1[,10]), 10]='GILGFVFTL'
d1[grep('Antigen:Yellow fever virus', d1[,10]), 10]='LLWNGPMAV'

# Match Embedding Matrix
# Assuming Mat is properly loaded and formatted, otherwise, skip this step or create a placeholder matrix.
# Mat=Mat[,c(1,6:101)]
# Mat=unique(Mat)
# rownames(Mat)=Mat[,1]
# Mat=Mat[,2:97]

# d1=cbind(d1, Mat[d1[,1],])

# Selection process
tt0=table(d1[,2], d1[,10])
LL=apply(tt0,1,function(x)length(which(x>0)))
names(which(LL==1))-> vv.tcrG
nn.pure=rownames(tt0[vv.tcrG,])
vv.select=grep('GIL', d1[,6])

nna=unique(d1[vv.select,2])
nna=intersect(nna,nn.pure)
nna1=sample(nna, 20)
ddc.s=d1[which(d1[,2] %in% nna1),]

# Distance matrix and neighbor-joining tree
Mat.cor=cor(t(ddc.s[, 8:103]))
Mat.dist=sqrt(1-Mat.cor)

tr = nj(Mat.dist)
tr$tip.label = gsub('\\.[1-9]{1}$','', tr$tip.label)

# High-contrast colors
COLs= c("dodgerblue2","gold", "green4", "#6A3D9A", "#FF7F00", "cyan", "skyblue2","#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon","orchid1","deeppink1","blue1","steelblue4", "darkturquoise","green1","yellow4","yellow3", "darkorange4","brown")

col.gr=COLs[1:20]
names(col.gr)=nna1
par(mar=c(0,0,0,0))
plot.phylo(tr, font=2, cex=0.9, tip.color=col.gr[as.character(ddc.s[,2])], align.tip.label=TRUE, x.lim=c(0,1.1))
