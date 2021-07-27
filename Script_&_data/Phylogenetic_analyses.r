library(diagram)
library(lubridate)
library(seraphim)
library(treeio)

"Exploiting genomic surveillance to map the spatio-temporal dispersal of SARS-CoV-2 spike mutations"

# 1. Loading and preparing the different GIS files used in the scripts
# 2. Preparing the input files for the discrete phylogeographic analyses 
# 3. Analysing the outputs of the preliminary discrete phylogeographic analysis 
# 4. Identifying the different clusters (clades following introduction events)
# 5. Preparing the continuous phylogeographic analysis (RRW, Cauchy model)
# 6. Running BEAST and building the maximum clade consensus (MCC) tree
# 7. Extracting spatio-temporal information embedded in MCC and posterior trees
# 8. Generating the dispersal history graphs (mapped MCC trees, 80% HPD polygons)
# 9. Exploring the spatio-temporal distribution of specific mutations

analysis = "TreeTime_011220"
	
data1 = read.csv("Sequences_metadata/SARS-CoV-2_KULeuven_011220.csv", sep=",")
data2 = read.csv("Sequences_metadata/SARS-CoV-2_ULiegeSeq_011220.csv", sep=",")
writingFiles = FALSE; showingPlots = FALSE

# 1. Loading and preparing the different GIS files used in the scripts

communes = shapefile("Shapefile_communes/Shapefile_post_codes.shp")
provinces = spTransform(raster::getData("GADM", country="BEL", level=2), crs(communes))
belgium = spTransform(raster::getData("GADM", country="BEL", level=0), crs(communes))
pop = projectRaster(raster("WorldPop_pop_raster.tif"), crs=crs(communes)); pop[] = log(pop[])
polIndex1 = which(provinces@data[,"NAME_2"]=="Liège"); maxArea = 0; polIndex2 = 0
for (i in 1:length(provinces@polygons[[polIndex1]]@Polygons))
	{
		if (maxArea < provinces@polygons[[polIndex1]]@Polygons[[i]]@area)
			{
				maxArea = provinces@polygons[[polIndex1]]@Polygons[[i]]@area; polIndex2 = i
			}
	}
pol = provinces@polygons[[polIndex1]]@Polygons[[polIndex2]]
p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
pol = sps; proj4string(pol) = crs(pop); pop_liege = raster::mask(crop(pop,pol),pol)
provinces_WGS = raster::getData("GADM", country="BEL", level=2)
polIndex1 = which(provinces_WGS@data[,"NAME_2"]=="Liège"); maxArea = 0; polIndex2 = 0
for (i in 1:length(provinces_WGS@polygons[[polIndex1]]@Polygons))
	{
		if (maxArea < provinces_WGS@polygons[[polIndex1]]@Polygons[[i]]@area)
			{
				maxArea = provinces_WGS@polygons[[polIndex1]]@Polygons[[i]]@area; polIndex2 = i
			}
	}
pol = provinces_WGS@polygons[[polIndex1]]@Polygons[[polIndex2]]
p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
pol = sps; proj4string(pol) = crs(raster("WorldPop_pop_raster.tif"))
pop_liege_WGS = raster::mask(crop(raster("WorldPop_pop_raster.tif"),pol),pol)

# 2. Preparing the input files for the discrete phylogeographic analyses 

tree = read.nexus(paste0(analysis,".tre"))
collection_dates = read.csv(paste0(analysis,".csv"))
if (grepl("TreeTime",analysis) == TRUE)
	{
		seqIDs = tree$tip.label; countries = rep(NA, length(seqIDs)); collectionDates = rep(NA, length(seqIDs))
		for (i in 1:length(seqIDs))
			{
				if (grepl("hCoV-19",seqIDs[i]))
					{
						countries[i] = unlist(strsplit(seqIDs[i],"\\/"))[2]
					}	else	{
						countries[i] = unlist(strsplit(seqIDs[i],"\\/"))[2]
					}
				collectionDates[i] = collection_dates[which(collection_dates[,"sequence_ID"]==seqIDs[i]),"collection_date"]
			}
		data = cbind(seqIDs,countries,collectionDates); colnames(data) = c("sequence_ID","country","collection_date")
	}	else	{
		data = read.csv(paste0(analysis,".csv"), sep=";")
	}
numberOfDifferentCountries = length(unique(data[,"country"]))
txt = c(); tab = c()
for (i in 1:length(tree$tip.label))
	{
		index = which(data[,1]==tree$tip.label[i])
		date = as.character(data[index,"collection_date"])
		if (date != "")
			{
				txt = c(txt, paste0(">",tree$tip.label[i]),"NNNN")
				if (grepl("hCoV-19",tree$tip.label[i]))
					{
						location = unlist(strsplit(tree$tip.label[i],"\\/"))[2]
					}	else	{
						location = unlist(strsplit(tree$tip.label[i],"\\/"))[2]
					}
				if (location != "Belgium") location = "other"
				tab = rbind(tab, cbind(tree$tip.label[i],location,date))
			}
	}
colnames(tab) = c("trait","location","collection_date")
if (writingFiles) write.table(tab, paste0(analysis,".txt"), row.names=F, quote=F, sep="\t")
if (writingFiles) write(txt, paste0(analysis,".fasta"))
txt = c(); tab = c()
for (i in 1:length(tree$tip.label))
	{
		index = which(data[,1]==tree$tip.label[i])
		date = as.character(data[index,"collection_date"])
		if (date != "")
			{
				txt = c(txt, paste0(">",tree$tip.label[i]),"NNNN")
				if (grepl("hCoV-19",tree$tip.label[i]))
					{
						location = unlist(strsplit(tree$tip.label[i],"\\/"))[2]
					}	else	{
						location = unlist(strsplit(tree$tip.label[i],"\\/"))[2]
					}
				if (grepl("ULG-",tree$tip.label[i])) location = "Liege"
				if (location != "Liege") location = "other"
				tab = rbind(tab, cbind(tree$tip.label[i],location,date))
			}
	}
colnames(tab) = c("trait","location","collection_date")
if (writingFiles) write.table(tab, paste0("DTA_on_Liege_province/",analysis,".txt"), row.names=F, quote=F, sep="\t")
if (writingFiles) write(txt, paste0("DTA_on_Liege_province/",analysis,".fasta"))

# 3. Analysing the outputs of the preliminary discrete phylogeographic analysis 

burnIn = 31; computingHPDInterval = FALSE # N.B.: long analysis
if (computingHPDInterval)
	{
		trees = scan(paste0(analysis,".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		indices1 = which(!grepl("tree STATE_",trees)); indices2 = which(grepl("tree STATE_",trees))
		mostRecentSamplingDate = max(decimal_date(ymd(collection_dates[,"collection_date"])))
		belgianBranches_list = rep(NA,length(trees)); belgianIntroductions_list = rep(NA,length(trees))
		belgianTipBranches_list = rep(NA,length(trees)); belgian_tMRCAs_list = list()
		for (i in (burnIn+1):length(indices2))
			{
				tree1 = trees[c(indices1[1:(length(indices1)-1)],indices2[i],indices1[length(indices1)])]
				write(tree1, paste0(analysis,"_sampled_tree_",i,".tree"))
				tree2 = readAnnotatedNexus(paste0(analysis,"_sampled_tree_",i,".tree"))
				belgianBranches = 0; belgianIntroductions = 0; belgianTipBranches = 0; belgian_tMRCAs = c()
				for (j in 1:dim(tree2$edge)[1])
					{
						if (tree2$annotations[[j]]$location == "Belgium")
							{
								belgianBranches = belgianBranches + 1
								index = which(tree2$edge[,2]==tree2$edge[j,1])
								if (tree2$annotations[[index]]$location != "Belgium")
									{
										belgianIntroductions = belgianIntroductions + 1
										tMRCA = mostRecentSamplingDate-nodeheight(tree2,tree2$edge[j,1])
										belgian_tMRCAs = c(belgian_tMRCAs, tMRCA)
									}
								if (!tree2$edge[j,2]%in%tree2$edge[,1])
									{
										belgianTipBranches = belgianTipBranches + 1
									}
							}
					}
				belgianBranches_list[i] = belgianBranches
				belgianIntroductions_list[i] = belgianIntroductions
				belgianTipBranches_list[i] = belgianTipBranches
				belgian_tMRCAs_list[[i]] = belgian_tMRCAs
				file.remove(paste0(analysis,"_sampled_tree_",i,".tree"))
			}
		saveRDS(belgian_tMRCAs_list, file="Belgian_tMRCAs.rds")
		quantiles = quantile(belgianIntroductions_list[!is.na(belgianIntroductions_list)],probs=c(0.025,0.975))
		cat("A minimum number of ",median(belgianIntroductions_list[!is.na(belgianIntroductions_list)])," lineage introductions (95% HPD interval = [",
			quantiles[1],"-",quantiles[2],"])"," identified from the global phylogenetic analysis of ",belgianTipBranches," SARS-CoV-2 sampled in Belgium",sep="")
		# Results for the analysis based on the IQTREE+TreeTime tree of the 01-12-20:
			# a minimum number of 338 lineage introductions (95% HPD interval = [312-353]) identified from the global phylogenetic analysis of 2114 SARS-CoV-2 sampled in Belgium
		trees = scan(paste0("DTA_on_Liege_province/",analysis,".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
		indices1 = which(!grepl("tree STATE_",trees)); indices2 = which(grepl("tree STATE_",trees))
		mostRecentSamplingDate = max(decimal_date(ymd(collection_dates[,"collection_date"])))
		liegeBranches_list = rep(NA,length(trees)); liegeIntroductions_list = rep(NA,length(trees))
		liegeTipBranches_list = rep(NA,length(trees)); liege_tMRCAs_list = list()
		for (i in (burnIn+1):length(indices2))
			{
				tree1 = trees[c(indices1[1:(length(indices1)-1)],indices2[i],indices1[length(indices1)])]
				write(tree1, paste0(analysis,"_sampled_tree_",i,".tree"))
				tree2 = readAnnotatedNexus(paste0(analysis,"_sampled_tree_",i,".tree"))
				liegeBranches = 0; liegeIntroductions = 0; liegeTipBranches = 0; liege_tMRCAs = c()
				for (j in 1:dim(tree2$edge)[1])
					{
						if (tree2$annotations[[j]]$location == "Liege")
							{
								liegeBranches = liegeBranches + 1
								index = which(tree2$edge[,2]==tree2$edge[j,1])
								if (tree2$annotations[[index]]$location != "Liege")
									{
										liegeIntroductions = liegeIntroductions + 1
										tMRCA = mostRecentSamplingDate-nodeheight(tree2,tree2$edge[j,1])
										liege_tMRCAs = c(liege_tMRCAs, tMRCA)
									}
								if (!tree2$edge[j,2]%in%tree2$edge[,1])
									{
										liegeTipBranches = liegeTipBranches + 1
									}
							}
					}
				liegeBranches_list[i] = liegeBranches
				liegeIntroductions_list[i] = liegeIntroductions
				liegeTipBranches_list[i] = liegeTipBranches
				liege_tMRCAs_list[[i]] = liege_tMRCAs
				file.remove(paste0(analysis,"_sampled_tree_",i,".tree"))
			}
		saveRDS(liege_tMRCAs_list, file="Belgian_tMRCAs.rds")
		quantiles = quantile(liegeIntroductions_list[!is.na(liegeIntroductions_list)],probs=c(0.025,0.975))
		cat("A minimum number of ",median(liegeIntroductions_list[!is.na(liegeIntroductions_list)])," lineage introductions (95% HPD interval = [",
			quantiles[1],"-",quantiles[2],"])"," identified from the global phylogenetic analysis of ",liegeTipBranches," SARS-CoV-2 sampled in Liège",sep="")
		# Results for the analysis based on the IQTREE+TreeTime tree of the 01-12-20:
			# a minimum number of 244 lineage introductions (95% HPD interval = [239-250]) identified from the global phylogenetic analysis of 869 SARS-CoV-2 sampled in Liège
	}
if (showingPlots)
	{
		tree = readAnnotatedNexus(paste0("DTA_on_Liege_province/",analysis,".tree"))
		samplingDates = decimal_date(ymd(gsub("\\/","-",tab[,"collection_date"]))); mostRecentSamplingYear = max(samplingDates)
		cols = rep("gray30",dim(tree$edge)[1]); lwds = rep(0.1,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="Liege") & (tree$annotations[[i]]$location=="Liege"))
							{
								cols[i] = rgb(222,67,39,255,maxColorValue=255); lwds[i] = 0.4
							}
					}
			}
		pdf("Figure_1a_NEW.pdf", width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		plot(tree, type="fan", show.tip.label=F, show.node.label=F, edge.width=lwds, cex=0.6, align.tip.label=3, col="gray30", edge.color=cols)
		for (i in 1:dim(tree$edge)[1])
			{
				if ((!tree$edge[i,2]%in%tree$edge[,1]) & (grepl("ULG-",tree$tip.label[tree$edge[i,2]])))
					{
						nodelabels(node=tree$edge[i,2], pch=16, cex=0.3, col=rgb(222,67,39,255,maxColorValue=255))
						nodelabels(node=tree$edge[i,2], pch=1, cex=0.3, col="gray30", lwd=0.5)
					}
				if (tree$annotations[[i]]$location == "Liege")
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if (tree$annotations[[index]]$location != "Liege")
							{
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col=rgb(222,67,39,255,maxColorValue=255))
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
							}
					}
			}
		add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
		cols = rep("gray50",dim(tree$edge)[1]); lwds = rep(0.05,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="Liege") & (tree$annotations[[i]]$location=="Liege"))
							{
								cols[i] = rgb(222,67,39,255,maxColorValue=255); lwds[i] = 0.4
							}
					}
			}
		dev.off()
		pdf("Figure_1b_NEW.pdf", width=4, height=4); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		ISO3_Europe = c("AUT","BEL","BGR","HRV","CYP","CZE","DNK","EST","FIN","FRA","DEU","GRC","HUN","IRL","ITA",
						"LVA","LTU","LUX","MLT","NLD","POL","PRT","ROU","SVK","SVN","ESP","SWE","GBR","NOR","CHE")
		shp = gSimplify(subset(shapefile("World_map_shapefile/World_map_shapefile.shp"),ISO3%in%ISO3_Europe),0.02)
		index = which(subset(shapefile("World_map_shapefile/World_map_shapefile.shp"),ISO3%in%ISO3_Europe)@data[,"ISO3"]=="BEL")
		liege = subset(raster::getData("GADM",country="BEL",level=2),NAME_2=="Liège")
		cols = rep("gray90", length(ISO3_Europe)); cols[index] = "gray70"
		plot(shp, border="gray30", lwd=0.2, col=cols); plot(liege, add=T, col=rgb(222,67,39,255,maxColorValue=255), border=NA)
		dev.off()
		pdf("Figure_1c_NEW.pdf", width=3.0, height=2.5); par(oma=c(0,0,0,0), mar=c(1,2,0,1), lwd=0.1)
		ats = c(decimal_date(ymd("2020-01-01")), decimal_date(ymd("2020-03-01")), decimal_date(ymd("2020-05-01")), decimal_date(ymd("2020-07-01")),
		decimal_date(ymd("2020-09-01")), decimal_date(ymd("2020-11-01")), decimal_date(ymd("2020-12-01")))
		labels = c("","01-03","01-05","01-07","01-09","01-11","")
		tab1 = read.table(paste0("DTA_on_Liege_province/",analysis,".txt"), head=T)
		tab1 = tab1[grepl("ULG-",tab1[,"trait"]),]; dates1 = decimal_date(ymd(tab1[,"collection_date"]))		
		tab2 = read.csv("https://epistat.sciensano.be/Data/COVID19BE_HOSP.csv")
		tab2 = tab2[which(tab2[,"PROVINCE"]=="Liège"),]; dates2 = decimal_date(ymd(tab2[,"DATE"]))
		values1 = tab2[,"NEW_IN"]; values2 = values1
		for (i in 1:length(values1))
			{
				indices = seq(i-3,i+3,1); indices = indices[which((indices[]>0)&(indices[]<=length(values1)))]; values2[i] = mean(values1[indices])
			}
		hist(dates1, breaks=50, border="gray30", col=rgb(222,67,39,150,maxColorValue=255), axes=F, ann=F, freq=T)
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30", at=ats, label=labels)
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		title(ylab="sampling frequency", cex.lab=0.85, mgp=c(1.2,0,0), col.lab=rgb(222,67,39,255,maxColorValue=255))
		par(new=T); lines(dates2, values2, col="gray30", lwd=1)
		axis(side=4, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.15,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		dev.off()
	}
if (showingPlots)
	{
		tree = readAnnotatedNexus(paste0(analysis,".tree"))
		samplingDates = decimal_date(ymd(gsub("\\/","-",tab[,"collection_date"]))); mostRecentSamplingYear = max(samplingDates)
		cols = rep("gray30",dim(tree$edge)[1]); lwds = rep(0.1,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="Belgium") & (tree$annotations[[i]]$location=="Belgium"))
							{
								cols[i] = "chartreuse3"; lwds[i] = 0.4
							}
					}
			}
		pdf("Figure_S1_NEW.pdf", width=7, height=7); par(oma=c(0,0,0,0), mar=c(0,0,0,0.0), lwd=0.1)
		plot(tree, type="fan", show.tip.label=F, show.node.label=F, edge.width=lwds, cex=0.6, align.tip.label=3, col="gray30", edge.color=cols)
		for (i in 1:dim(tree$edge)[1])
			{
				if ((!tree$edge[i,2]%in%tree$edge[,1]) & (grepl("Belgium",tree$tip.label[tree$edge[i,2]])))
					{
						nodelabels(node=tree$edge[i,2], pch=16, cex=0.3, col="chartreuse3")
						nodelabels(node=tree$edge[i,2], pch=1, cex=0.3, col="gray30", lwd=0.5)
					}
				if (tree$annotations[[i]]$location == "Belgium")
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if (tree$annotations[[index]]$location != "Belgium")
							{
								nodelabels(node=tree$edge[i,2], pch=16, cex=0.6, col="chartreuse3")
								nodelabels(node=tree$edge[i,2], pch=1, cex=0.6, col="gray30", lwd=0.5)
							}
					}
			}
		add.scale.bar(x=0.0, y=-0.01, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.7)
		cols = rep("gray50",dim(tree$edge)[1]); lwds = rep(0.05,dim(tree$edge)[1])
		for (i in 1:dim(tree$edge)[1])
			{
				if (tree$edge[i,1]%in%tree$edge[,2])
					{
						index = which(tree$edge[,2]==tree$edge[i,1])
						if ((tree$annotations[[index]]$location=="Belgium") & (tree$annotations[[i]]$location=="Belgium"))
							{
								cols[i] = "chartreuse3"; lwds[i] = 0.4
							}
					}
			}
		dev.off()
	}

# 4. Identifying the different clusters (clades following introduction events)

tree = readAnnotatedNexus(paste0(analysis,".tree"))
belgianBranches = c(); belgianIntroductions = c()
belgianTipBranches = c(); sampledSequences = c()
for (i in 1:dim(tree$edge)[1])
	{
		if (tree$annotations[[i]]$location == "Belgium")
			{
				belgianBranches = c(belgianBranches,i)
				index = which(tree$edge[,2]==tree$edge[i,1])
				if (tree$annotations[[index]]$location != "Belgium")
					{
						belgianIntroductions = c(belgianIntroductions, i)
					}
				if (!tree$edge[i,2]%in%tree$edge[,1])
					{
						belgianTipBranches = c(belgianTipBranches, i)
						sampledSequences = c(sampledSequences, tree$tip.label[tree$edge[i,2]])
					}
			}
	}
for (i in 1:length(belgianIntroductions))
	{
		if (i == 1) clusters1 = list()
		if (tree$edge[belgianIntroductions[i],2]%in%tree$edge[,1])
			{
				subtree = tree_subset(tree, tree$edge[belgianIntroductions[i],2], levels_back=0)
				clusters1[[i]] = gsub("'","",subtree$tip.label)
			}	else		{
				clusters1[[i]] = gsub("'","",tree$tip.label[tree$edge[belgianIntroductions[i],2]])
			}
	}
for (i in 2:length(clusters1))
	{
		for (j in 1:(i-1))
			{
				if (sum(clusters1[[i]]%in%clusters1[[j]]) == length(clusters1[[i]]))
					{
						clusters1[[j]] = clusters1[[j]][which(!clusters1[[j]]%in%clusters1[[i]])]
					}
				if (sum(clusters1[[j]]%in%clusters1[[i]]) == length(clusters1[[j]]))
					{
						clusters1[[i]] = clusters1[[i]][which(!clusters1[[i]]%in%clusters1[[j]])]
					}
			}
	}
sampledSequences = gsub("'","",sampledSequences)
if (!file.exists(paste0("Sampling_Belgium.csv")))
	{
		samplingData = matrix(nrow=length(sampledSequences), ncol=5)
		colnames(samplingData) = c("sequence_ID","collection_date","ZIP","longitude","latitude")
		samplingData[,"sequence_ID"] = sampledSequences
		for (i in 1074:dim(samplingData)[1])
			{
				index = which(data[,"sequence_ID"]==samplingData[i,"sequence_ID"])
				date = ymd(gsub("\\/","-",data[index,"collection_date"]))
				samplingData[i,"collection_date"] = decimal_date(date)
				if (grepl("rega",gsub("Rega","rega",samplingData[i,"sequence_ID"])))
					{
						temp = unlist(strsplit(samplingData[i,"sequence_ID"],"\\/"))
						ID = temp[which(grepl("rega",gsub("Rega","rega",temp)))]
					}
				if (!grepl("rega",gsub("Rega","rega",samplingData[i,"sequence_ID"]))) ID = unlist(strsplit(samplingData[i,"sequence_ID"],"\\/"))[3]
				index1 = which(grepl(gsub("Rega","rega",ID),data1[,"sequence_ID"]))
				if (length(index1) == 1)
					{
						samplingData[i,"ZIP"] = data1[index1,"ZIP"]
						indices = which(communes@data[,"nouveau_PO"]==data1[index1,"ZIP"])
						if (length(indices) > 0)
							{
								maxArea = 0; polIndex1 = 0; polIndex2 = 0
								for (j in 1:length(indices))
									{
										for (k in 1:length(communes@polygons[[indices[j]]]@Polygons))
											{
												if (maxArea < communes@polygons[[indices[j]]]@Polygons[[k]]@area)
													{
														maxArea = communes@polygons[[indices[j]]]@Polygons[[k]]@area; polIndex1 = indices[j]; polIndex2 = k
													}
											}
									}
								pol = communes@polygons[[polIndex1]]@Polygons[[polIndex2]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = communes@proj4string
								# samplingData[i,c("longitude","latitude")] = coordinates(pol) # to avoid a jitter:
								samplingData[i,c("longitude","latitude")] = spsample(pol, 1, type="random")@coords
							}
					}
				ID = unlist(strsplit(samplingData[i,"sequence_ID"],"\\/"))[3]
				index1 = which(grepl(ID,data2[,"sequence_ID"]))
				if (length(index1) == 1)
					{
						samplingData[i,"ZIP"] = as.character(data2[index1,"ZIP"])
						indices = which(communes@data[,"nouveau_PO"]==data2[index1,"ZIP"])
						if (length(indices) > 0)
							{
								maxArea = 0; polIndex1 = 0; polIndex2 = 0
								for (j in 1:length(indices))
									{
										for (k in 1:length(communes@polygons[[indices[j]]]@Polygons))
											{
												if (maxArea < communes@polygons[[indices[j]]]@Polygons[[k]]@area)
													{
														maxArea = communes@polygons[[indices[j]]]@Polygons[[k]]@area; polIndex1 = indices[j]; polIndex2 = k
													}
											}
									}
								pol = communes@polygons[[polIndex1]]@Polygons[[polIndex2]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; proj4string(pol) = communes@proj4string
								# samplingData[i,c("longitude","latitude")] = coordinates(pol) # to avoid a jitter:
								samplingData[i,c("longitude","latitude")] = spsample(pol, 1, type="random")@coords
							}
					}
			}
		print(samplingData[which(is.na(samplingData[,"ZIP"])),"sequence_ID"])
		write.csv(samplingData, "Sampling_Belgium.csv", quote=F, row.names=F)
	}	
samplingData = read.csv("Sampling_Belgium.csv", head=T)
for (i in 1:length(belgianIntroductions))
	{
		tab = c()
		if (i == 1)
			{
				clusters2 = list(); centroids = list()
			}
		for (j in 1:length(clusters1[[i]]))
			{
				index = which(samplingData[,"sequence_ID"]==clusters1[[i]][j])
				if (length(index) == 1)
					{
						line = cbind(as.numeric(samplingData[index,"collection_date"]),as.numeric(samplingData[index,"longitude"]),as.numeric(samplingData[index,"latitude"]))
						row.names(line) = clusters1[[i]][j]; tab = rbind(tab, line)
					}
			}
		colnames(tab) = c("collection_date","longitude","latitude"); clusters2[[i]] = tab
		centroids[[i]] = cbind(mean(tab[!is.na(tab[,"longitude"]),"longitude"]), mean(tab[!is.na(tab[,"latitude"]),"latitude"]))
	}
clusterSizes = rep(NA, length(clusters1))
collectionDates = c()
for (i in 1:length(clusters1))
	{
		clusterSizes[i] = length(clusters1[[i]])
		collectionDates = c(collectionDates, clusters2[[i]][,"collection_date"])
	}
if (showingPlots)
	{
		collectionDates_filetered = collectionDates
		dev.new(width=3.3, height=8); par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(2,2,1,1), lwd=0.2, col="gray30")
		hist(clusterSizes, breaks=50, axes=F, ann=F, title=NULL, col="#66CD0099", border="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		hist(collectionDates_filetered, breaks=65, axes=F, ann=F, title=NULL, col="#66CD0099", border="gray30")
		axis(side=2, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.65, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=decimal_date(ymd(c("2020-03-01","2020-05-01","2020-07-01","2020-09-01","2020-11-01"))),
			 labels=c("01-03","01-05","01-07","01-09","01-11"))
	}

# 5. Preparing the continuous phylogeographic analysis (RRW, Cauchy model)

template = scan("RRW_template_2.xml", what="", sep="\n", quiet=T, blank.lines.skip=F); xml = c()
sink(file=paste0("Phylogeographic_runs/All_clades_NEW.xml"))
for (i in 1:length(template))
	{
		cat(template[i],"\n")
		if (grepl("Insert taxa blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<taxa id=\"taxa_",j,"\">","\n"))
								for (k in 1:dim(clusters2[[j]])[1])
									{
										if (!is.na(clusters2[[j]][k,"longitude"]))
											{
												cat(paste0("\t\t<taxon id=\"",row.names(clusters2[[j]])[k],"\">","\n"))
												cat(paste0("\t\t\t<date value=\"",clusters2[[j]][k,"collection_date"],"\" direction=\"forwards\" units=\"years\"/>","\n"))
												cat("\t\t\t<attr name=\"latitude\">\n")
												cat(paste0("\t\t\t\t",clusters2[[j]][k,"latitude"],"\n"))
												cat("\t\t\t</attr>\n")
												cat("\t\t\t<attr name=\"longitude\">\n")
												cat(paste0("\t\t\t\t",clusters2[[j]][k,"longitude"],"\n"))
												cat("\t\t\t</attr>\n")
												cat("\t\t\t<attr name=\"coordinates\">\n")
												cat(paste0("\t\t\t\t",clusters2[[j]][k,"latitude"]," ",clusters2[[j]][k,"longitude"],"\n"))
												cat("\t\t\t</attr>\n")
												cat("\t\t</taxon>\n")
											}
									}
								cat("\t</taxa>","\n")
							}
					}
			}
		if (grepl("Insert alignment blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<alignment id=\"alignment_",j,"\" dataType=\"nucleotide\">","\n"))
								for (k in 1:dim(clusters2[[j]])[1])
									{
										if (!is.na(clusters2[[j]][k,"longitude"]))
											{
												cat("\t\t<sequence>\n")
												cat(paste0("\t\t\t<taxon idref=\"",row.names(clusters2[[j]])[k],"\"/>","\n"))
												cat("\t\t\tNNNN\n")
												cat("\t\t</sequence>\n")
											}
									}
								cat("\t</alignment>","\n")
							}
					}
			}
		if (grepl("Insert pattern blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<patterns id=\"patterns_",j,"\" from=\"1\" strip=\"false\">","\n"))
								cat(paste0("\t\t<alignment idref=\"alignment_",j,"\"/>","\n"))
								cat("\t</patterns>","\n")
							}
					}
			}
		if (grepl("Insert starting tree blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								tre = tree_subset(tree, tree$edge[belgianIntroductions[j],2], levels_back=0)
								tips = row.names(clusters2[[j]]); tips = tips[which(!is.na(clusters2[[j]][,"longitude"]))]
								tips_to_drop = tre$tip.label[which(!gsub("'","",tre$tip.label)%in%tips)]
								if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
								write.tree(tre, paste0("Phylogeographic_runs/Clade_",i,".tre"))
								tre = scan(paste0("Phylogeographic_runs/Clade_",i,".tre"), what="", sep="\n", quiet=T)
								txt = c("#NEXUS","begin trees;",paste0("\ttree tree_1 = [&R] ",tre),"end;")
								write(txt, paste0("Phylogeographic_runs/Clade_",j,".tre"))
								cat(paste0("\t<empiricalTreeDistributionModel id=\"treeModel_",j,"\" fileName=\"Clade_",j,".tre\">","\n"))
								cat(paste0("\t\t<taxa idref=\"taxa_",j,"\"/>","\n"))
								cat("\t</empiricalTreeDistributionModel>","\n")
							}
					}
			}
		if (grepl("Insert tree model blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<treeModel id=\"treeModel_",j,"\">","\n"))
								cat(paste0("\t\t<coalescentTree idref=\"startingTree_",j,"\"/>","\n"))
								cat("\t\t<rootHeight>","\n")
								cat(paste0("\t\t\t<parameter id=\"treeModel.rootHeight_",j,"\"/>","\n"))
								cat("\t\t</rootHeight>","\n")
								cat("\t\t<nodeHeights internalNodes=\"true\">","\n")
								cat(paste0("\t\t\t<parameter id=\"treeModel.internalNodeHeights_",j,"\"/>","\n"))
								cat("\t\t</nodeHeights>","\n")
								cat("\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">","\n")
								cat(paste0("\t\t\t<parameter id=\"treeModel.allInternalNodeHeights_",j,"\"/>","\n"))
								cat("\t\t</nodeHeights>","\n")
								cat("\t</treeModel>","\n")
							}
					}
			}
		if (grepl("Insert arbitraryBranchRates blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<arbitraryBranchRates id=\"coordinates.diffusion.branchRates",j,"\">","\n"))
								cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
								cat("\t\t<rates>","\n")
								cat(paste0("\t\t\t<parameter id=\"coordinates.diffusion.rates",j,"\" lower=\"0.0\"/>","\n"))
								cat("\t\t</rates>","\n")
								cat("\t</arbitraryBranchRates>","\n")
							}
					}
			}
		if (grepl("Insert distributionLikelihood blocks 1",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<distributionLikelihood id=\"coordinates.diffusion.prior",j,"\">","\n"))
								cat("\t\t<data>","\n")
								cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
								cat("\t\t</data>","\n")
								cat("\t\t<distribution>","\n")
								cat(paste0("\t\t\t<onePGammaDistributionModel>","\n"))
								cat("\t\t\t\t<shape>","\n")
								cat("\t\t\t\t\t<parameter value=\"0.5\"/>","\n")
								cat("\t\t\t\t</shape>","\n")
								cat("\t\t\t</onePGammaDistributionModel>","\n")
								cat("\t\t</distribution>","\n")
								cat("\t</distributionLikelihood>","\n")
							}
					}
			}
		if (grepl("Insert coordinates.traitLikelihood blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<multivariateTraitLikelihood id=\"coordinates.traitLikelihood",j,"\" traitName=\"coordinates\" useTreeLength=\"true\" scaleByTime=\"true\" reportAsMultivariate=\"true\" reciprocalRates=\"true\" integrateInternalTraits=\"true\">","\n"))
								cat("\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
								cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>"))
								cat("\t\t<traitParameter>","\n")
								cat(paste0("\t\t\t<parameter id=\"leaf.coordinates",j,"\"/>","\n"))
								cat("\t\t</traitParameter>","\n")
								cat("\t\t<conjugateRootPrior>","\n")
								cat("\t\t\t<meanParameter>","\n")
								cat("\t\t\t\t<parameter value=\"0.0 0.0\"/>","\n")
								cat("\t\t\t</meanParameter>","\n")
								cat("\t\t\t<priorSampleSize>","\n")
								cat("\t\t\t\t<parameter value=\"0.000001\"/>","\n")
								cat("\t\t\t</priorSampleSize>","\n")
								cat("\t\t</conjugateRootPrior>","\n")
								cat(paste0("\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
								cat("\t</multivariateTraitLikelihood>","\n")
							}
					}
			}
		if (grepl("Insert continuousDiffusionStatistic blocks 1",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<continuousDiffusionStatistic id=\"coordinates.diffusionRate",j,"\" greatCircleDistance=\"true\">","\n"))
								cat(paste0("\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
								cat("\t</continuousDiffusionStatistic>","\n")
							}
					}
			}
		if (grepl("Insert scaleOperator blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"30\">","\n"))
								cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
								cat("\t\t</scaleOperator>","\n")
							}
					}
			}
		if (grepl("Insert precisionGibbsOperator blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t<precisionGibbsOperator weight=\"2\">","\n"))
								cat(paste0("\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
								cat("\t\t\t<multivariateWishartPrior idref=\"coordinates.precisionPrior\"/>","\n")
								cat("\t\t</precisionGibbsOperator>","\n")
							}
					}
			}
		if (grepl("Insert distributionLikelihood blocks 2",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<distributionLikelihood idref=\"coordinates.diffusion.prior",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("Insert multivariateTraitLikelihood blocks 1",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("Insert continuousDiffusionStatistic blocks 2",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<continuousDiffusionStatistic idref=\"coordinates.diffusionRate",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("Insert multivariateTraitLikelihood blocks 2",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("<!-- Insert logTree blocks -->",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t<logTree id=\"treeFileLog",j,"\" logEvery=\"1000\" nexusFormat=\"true\" fileName=\"Clade_",j,".trees\" sortTranslationTable=\"true\">","\n"))
								cat(paste0("\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
								cat("\t\t\t<joint idref=\"joint\"/>","\n")
								cat("\t\t\t<trait name=\"coordinates\" tag=\"coordinates\">","\n")
								cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
								cat("\t\t\t</trait>","\n")
								cat("\t\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
								cat("\t\t\t<trait name=\"rate\" tag=\"coordinates.rate\">","\n")
								cat(paste0("\t\t\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
								cat("\t\t\t</trait>","\n")
								cat("\t\t</logTree>","\n")
							}
					}
			}
	}
sink(NULL)

# 6. Running BEAST and building the maximum clade consensus (MCC) tree

source("MCC_tree_extraction.r"); analyses = c(); runningNewAnalyses = FALSE
wd = getwd(); setwd(paste0(wd,"/Phylogeographic_runs/"))
for (i in 1:length(clusters2))
	{
		if ((dim(clusters2[[i]])[1] >= 3)&(sum(!is.na(clusters2[[i]][,"longitude"])) >= 3)) analyses = c(analyses, paste0("Clade_",i))
	}
if (runningNewAnalyses)
	{
		system("java -jar beast_1104.jar All_clades.xml", ignore.stdout=T, ignore.stderr=F)
		for (i in 1:length(analyses))
			{
system(paste0("BEAST_1104/bin/treeannotator -burninTrees 501 -heights keep ",analyses[i],".trees ",analyses[i],".tree"), ignore.stdout=F, ignore.stderr=F)
			}
	}
setwd(wd)

# 7. Extracting spatio-temporal information embedded in MCC and posterior trees

wd = getwd(); setwd(paste0(wd,"/Phylogeographic_runs/"))
for (i in 1:length(analyses))
	{
		if (!file.exists(paste0(analyses[i],".csv")))
			{
				index = as.numeric(unlist(strsplit(analyses[i],"_"))[2])
				mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collection_date"])
				mcc_tre = readAnnotatedNexus(paste0(analyses[i],".tree")); dates = c()
				mcc_tab = MCC_tree_extraction(mcc_tre, mostRecentSamplingDatum)
				write.csv(mcc_tab, paste0(analyses[i],".csv"), row.names=F, quote=F)
			}
	}
for (i in 1:length(analyses))
	{
		mcc_tab = read.csv(paste0(analyses[i],".csv"), head=T)
		if (!"tipLabel"%in%colnames(mcc_tab))
			{
				tipLabels = matrix(nrow=dim(mcc_tab)[1], ncol=1); colnames(tipLabels) = "tipLabel"
				for (j in 1:dim(mcc_tab)[1])
					{
						if (!mcc_tab[j,"node2"]%in%mcc_tab[,"node1"])
							{
								index = which((round(samplingData[,"longitude"])==round(mcc_tab[j,"endLon"]))&(round(samplingData[,"latitude"])==round(mcc_tab[j,"endLat"])))
								if (length(index) != 1)
									{
										print(c(i,j))
									}	else	{
										tipLabels[j,1] = samplingData[index,"sequence_ID"]
									}
							}	
					}
				mcc_tab = cbind(mcc_tab, tipLabels)
				write.csv(mcc_tab, paste0(analyses[i],".csv"), row.names=F, quote=F)
			}
	}
nberOfTreesToSample = 1000; randomSampling = FALSE; coordinateAttributeName = "coordinates"; nberOfCores = 10
for (i in 1:length(analyses))
	{
		localTreesDirectory = paste0(analyses[i],"_ext")
		if (!file.exists(paste0(localTreesDirectory,"/TreeExtractions_1.csv")))
			{
				index = as.numeric(unlist(strsplit(analyses[i],"_"))[2]); burnIn = 501
				mostRecentSamplingDatum = max(clusters2[[index]][which(!is.na(clusters2[[index]][,"longitude"])),"collection_date"])
				allTrees = scan(file=paste0(analyses[i],".trees"), what="", sep="\n", quiet=T, blank.lines.skip=F)
				treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName, nberOfCores)
			}
	}
for (i in 1:length(analyses))
	{
		for (j in 1:nberOfTreesToSample)
			{
				tab = read.csv(paste0(analyses[i],"_ext/TreeExtractions_",j,".csv"), head=T)
				if (!"tipLabel"%in%colnames(tab))
					{
						tipLabels = matrix(nrow=dim(tab)[1], ncol=1); colnames(tipLabels) = "tipLabel"
						for (k in 1:dim(tab)[1])
							{
								if (!tab[k,"node2"]%in%tab[,"node1"])
									{
										index = which((round(samplingData[,"longitude"])==round(tab[k,"endLon"]))&(round(samplingData[,"latitude"])==round(tab[k,"endLat"])))
										if (length(index) != 1)
											{
												print(c(i,j,k))
											}	else	{
												tipLabels[k,1] = samplingData[index,"sequence_ID"]
											}
									}	
							}
						tab = cbind(tab, tipLabels)
						write.csv(tab, paste0(analyses[i],"_ext/TreeExtractions_",j,".csv"), row.names=F, quote=F)
					}
			}
	}
cladesToExclude = c()
for (i in 1:length(analyses))
	{
		tab = read.csv(paste0(analyses[i],".csv"), head=T)
		if (i == 1)
			{
				all = tab
			}	else	{
				if (!analyses[i]%in%cladesToExclude)
					{
						maxNodeID = max(all[,c("node1","node2")])
						tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
						all = rbind(all, tab)
					}
			}
	}
write.csv(all, "All_clades.csv", row.names=F, quote=F)
dir.create(file.path("All_clades_ext1"), showWarnings=F)
dir.create(file.path("All_clades_ext2"), showWarnings=F)
dir.create(file.path("1st_wave_ext2"), showWarnings=F)
dir.create(file.path("2nd_wave_ext2"), showWarnings=F)
dir.create(file.path("All_Liege_ext1"), showWarnings=F)
dir.create(file.path("All_Liege_ext2"), showWarnings=F)
dir.create(file.path("1st_wave_Liege"), showWarnings=F)
dir.create(file.path("2nd_wave_Liege"), showWarnings=F)
nberOfExtractionFiles = nberOfTreesToSample
for (i in 1:nberOfExtractionFiles)
	{
		for (j in 1:length(analyses))
			{
				tab = read.csv(paste0(analyses[j],"_ext/TreeExtractions_",i,".csv"), head=T)
				if (j == 1)
					{
						all = tab
					}	else	{
						if (!analyses[i]%in%cladesToExclude)
							{
								maxNodeID = max(all[,c("node1","node2")])
								tab[,c("node1","node2")] = tab[,c("node1","node2")]+maxNodeID
								all = rbind(all, tab)
							}
					}
			}
		write.csv(all, paste0("All_clades_ext1/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		vS1 = raster::extract(pop_liege, all[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege, all[,c("endLon","endLat")])
		lie = all[which((!is.na(vS1))&(!is.na(vS2))),]
		write.csv(lie, paste0("All_Liege_ext1/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		temp1 = all[,c("startLon","startLat")]; temp2 = all[,c("endLon","endLat")]
		coordinates(temp1) = ~ startLon + startLat; crs(temp1) = crs(communes)
		coordinates(temp2) = ~ endLon + endLat; crs(temp2) = crs(communes)
		temp1 = spTransform(temp1, CRS("+init=epsg:4326"))@coords
		temp2 = spTransform(temp2, CRS("+init=epsg:4326"))@coords
		all[,c("startLon","startLat")] = temp1; all[,c("endLon","endLat")] = temp2
		write.csv(all, paste0("All_clades_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		bf1 = all[which(all[,"endYear"]<decimal_date(dmy("01-06-2020"))),]
		af1 = all[which(all[,"endYear"]>=decimal_date(dmy("31-08-2020"))),]
		write.csv(bf1, paste0("1st_wave_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		write.csv(af1, paste0("2nd_wave_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		vS1 = raster::extract(pop_liege_WGS, all[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege_WGS, all[,c("endLon","endLat")])
		lie = all[which((!is.na(vS1))&(!is.na(vS2))),]
		write.csv(lie, paste0("All_Liege_ext2/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		vS1 = raster::extract(pop_liege_WGS, bf1[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege_WGS, bf1[,c("endLon","endLat")])
		bf2 = bf1[which((!is.na(vS1))&(!is.na(vS2))),]
		vS1 = raster::extract(pop_liege_WGS, af1[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege_WGS, af1[,c("endLon","endLat")])
		af2 = af1[which((!is.na(vS1))&(!is.na(vS2))),]	
		write.csv(bf2, paste0("1st_wave_Liege/TreeExtractions_",i,".csv"), row.names=F, quote=F)
		write.csv(af2, paste0("2nd_wave_Liege/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}
setwd(wd)

# 8. Generating the dispersal history graphs (mapped MCC trees, 80% HPD polygons)

localTreesDirectory = paste0("Phylogeographic_runs/All_clades_ext1"); nberOfExtractionFiles = 1000
percentage = 80; prob = percentage/100; precision = 1/(365/7); croppingPolygons = TRUE
mcc = read.csv("Phylogeographic_runs/All_clades.csv", head=T); startDatum = min(mcc[,"startYear"])
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
if (showingPlots)
	{
		colourScale = rev(colorRampPalette(brewer.pal(11,"BrBG"))(141)[16:116])
		colourScale = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(141)[16:116])
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		startYears_indices = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		startYears_colours = colourScale[startYears_indices]
		endYears_colours = colourScale[endYears_indices]
		polygons_colours = rep(NA, length(polygons))
		cexNode = 0.6; LWD = 1.0
		for (i in 1:length(polygons))
			{
				date = as.numeric(names(polygons[[i]]))
				polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
				polygons_colours[i] = paste0(colourScale[polygon_index],"40")
			}
		firstTimePeriod = TRUE; secondTimePeriod = FALSE
		firstTimePeriod = FALSE; secondTimePeriod = TRUE
		firstTimePeriod = FALSE; secondTimePeriod = FALSE
		pdf("Figure_2_NEW.pdf", width=7.3, height=6) # dev.new(width=7.3, height=6)
		par(oma=c(0,0,0,0), mar=c(1,1,1,1), lwd=0.2, col="gray30")
		cols = c(colorRampPalette(brewer.pal(9,"YlGnBu"))(161)[1:101])
		cols = c(colorRampPalette(brewer.pal(9,"Greys"))(201)[1:101])
		if ((firstTimePeriod == TRUE)|(secondTimePeriod == TRUE))
			{
				# plot(pop, col=cols, axes=F, ann=F, box=F, legend=F)
				# plot(provinces, border="white", col=NA, add=T, lwd=LWD)
				plot(provinces, border=NA, col="gray95", lwd=LWD)
			}	else	{
				plot(provinces, border=NA, col="gray95", lwd=LWD)
			}
		if ((firstTimePeriod == FALSE)&(secondTimePeriod == FALSE))
			{
				for (i in 1:length(polygons))
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								polygons[[i]]@polygons[[j]] = maptools::checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]; crs(pol) = crs(belgium)
						if (croppingPolygons == TRUE) pol = crop(pol, belgium)
						plot(pol, axes=F, col=polygons_colours[i], add=T, border=NA)
					}
			}
		plot(provinces, border="white", col=NA, add=T, lwd=LWD)
		plot(belgium, border="gray30", col=NA, add=T, lwd=0.4); croppingPolygons = TRUE
		selectedBranches = 1:dim(mcc)[1]
		if (firstTimePeriod) selectedBranches = which(mcc[,"endYear"]<decimal_date(ymd("2020-06-01")))
		if (secondTimePeriod) selectedBranches = which(mcc[,"startYear"]>decimal_date(ymd("2020-08-31")))
		for (i in selectedBranches)
			{
				curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						    arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		if ((firstTimePeriod == FALSE)&(secondTimePeriod == FALSE))
			{
				for (i in 1:length(clusters2))
					{
						if (sum(!is.na(clusters2[[i]][,"longitude"])) < 3)
							{
								if (sum(!is.na(clusters2[[i]][,"longitude"])) == 2)
									{
										indices = which(!is.na(clusters2[[i]][,"longitude"]))
										if (firstTimePeriod)
											{
												indices = which((!is.na(clusters2[[i]][,"longitude"]))&(clusters2[[i]][,"collection_date"]<decimal_date(ymd("2020-03-18"))))
											}
										if (secondTimePeriod)
											{
												indices = which((!is.na(clusters2[[i]][,"longitude"]))&(clusters2[[i]][,"collection_date"]>decimal_date(ymd("2020-03-18"))))
											}
										if (length(indices) == 2)
											{
												curvedarrow(cbind(clusters2[[i]][indices[1],"longitude"],clusters2[[i]][indices[1],"latitude"]),
															cbind(clusters2[[i]][indices[2],"longitude"],clusters2[[i]][indices[2],"latitude"]),
															arr.length=0, arr.width=0, lwd=0.2, lty=2, lcol="gray10", arr.col=NA, arr.pos=F,
															curve=0.1, dr=NA, endhead=F)
											}
									}
							}
					}
			}
		for (i in 1:length(clusters2))
			{
				if (sum(!is.na(clusters2[[i]][,"longitude"])) < 3)
					{
						for (j in 1:dim(clusters2[[i]])[1])
							{
								if (!is.na(clusters2[[i]][j,"longitude"]))
									{
										plotTheNode = TRUE
										if ((firstTimePeriod==TRUE)&(clusters2[[i]][j,"collection_date"]>decimal_date(ymd("2020-06-01")))) plotTheNode = FALSE
										if ((secondTimePeriod==TRUE)&(clusters2[[i]][j,"collection_date"]<decimal_date(ymd("2020-08-31")))) plotTheNode = FALSE
										if (plotTheNode)
											{
												index = (((clusters2[[i]][j,"collection_date"]-minYear)/(maxYear-minYear))*100)+1
												points(clusters2[[i]][j,"longitude"], clusters2[[i]][j,"latitude"], pch=16, col=colourScale[index], cex=cexNode)
												points(clusters2[[i]][j,"longitude"], clusters2[[i]][j,"latitude"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
											}
									}
							}
					}
			}
		for (i in rev(selectedBranches))
			{
				if (!mcc[i,"node1"]%in%mcc[selectedBranches,"node2"])
					{
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[i], cex=cexNode)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
					}
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
				points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
			}
		selectedDates = decimal_date(ymd(c("2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01")))
		selectedLabels = c("01/03","01/04","01/05","01/06","01/07","01/08","01/09","01/10","01/11")
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
		plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.05,0.55,0.100,0.112),
			 legend.args=list(text="", cex=0.6, line=0.3, col="gray30"), horizontal=T,
		     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.6, col.axis="gray30", line=0, mgp=c(0,0.00,0),
		     at=selectedDates, labels=selectedLabels))
		dev.off()
	}
if (showingPlots)
	{
		polIndex1 = which(provinces@data[,"NAME_2"]=="Liège")
		maxArea = 0; polIndex2 = 0; cexNode = 0.7
		for (i in 1:length(provinces@polygons[[polIndex1]]@Polygons))
			{
				if (maxArea < provinces@polygons[[polIndex1]]@Polygons[[i]]@area)
					{
						maxArea = provinces@polygons[[polIndex1]]@Polygons[[i]]@area; polIndex2 = i
					}
			}
		pol = provinces@polygons[[polIndex1]]@Polygons[[polIndex2]]
		p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pol = sps; proj4string(pol) = communes@proj4string
		pop_liege = raster::mask(crop(pop,pol),pol)
		vS1 = raster::extract(pop_liege, mcc[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege, mcc[,c("endLon","endLat")])
		sub = mcc[which((!is.na(vS1))&(!is.na(vS2))),]		
		colourScale = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(141)[21:121])
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		startYears_indices = (((sub[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices = (((sub[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		startYears_colours = colourScale[startYears_indices]
		endYears_colours = colourScale[endYears_indices]
		pdf("Figure_2_NEW.pdf", width=11, height=8) # dev.new(width=11, height=8)
		par(mfrow=c(2,2), oma=c(0,2,1,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		for (i in 1:4)
			{
				cols = c(colorRampPalette(brewer.pal(9,"YlGnBu"))(161)[1:101])
				cols = c(colorRampPalette(brewer.pal(9,"Greys"))(201)[1:101])
				plot(pop_liege, col=cols, axes=F, ann=F, box=F, legend=F)
				plot(pol, border="gray30", col=NA, add=T, lwd=0.4)
				if (i == 1) selectedBranches = 1:dim(sub)[1]
				if (i == 2) selectedBranches = which(sub[,"endYear"]<decimal_date(dmy("01-06-2020")))
				if (i == 3) selectedBranches = which((sub[,"endYear"]>=decimal_date(dmy("01-06-2020")))&(sub[,"endYear"]<=decimal_date(dmy("31-08-2020"))))
				if (i == 4) selectedBranches = which(sub[,"endYear"]>decimal_date(dmy("31-08-2020")))
				for (j in selectedBranches)
					{
						curvedarrow(cbind(sub[j,"startLon"],sub[j,"startLat"]), cbind(sub[j,"endLon"],sub[j,"endLat"]), arr.length=0,
						    		arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
				for (j in rev(selectedBranches))
					{
						if (!sub[j,"node1"]%in%sub[selectedBranches,"node2"])
							{
								points(sub[j,"startLon"], sub[j,"startLat"], pch=16, col=startYears_colours[j], cex=cexNode)
								points(sub[j,"startLon"], sub[j,"startLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
							}
						points(sub[j,"endLon"], sub[j,"endLat"], pch=16, col=endYears_colours[j], cex=cexNode)
						points(sub[j,"endLon"], sub[j,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
					}
				for (j in 1:length(clusters2))
					{
						if (sum(!is.na(clusters2[[j]][,"longitude"])) < 3)
							{
								if (length(which(!is.na(clusters2[[j]][,"longitude"]))) > 2)
									{
										print(c(i,j))
										sub3 = clusters2[[j]][which(!is.na(clusters2[[j]][,"longitude"])),]
										vS3 = raster::extract(pop_liege, sub3[1,c("longitude","latitude")])
										vS4 = raster::extract(pop_liege, sub3[2,c("longitude","latitude")])
										if ((!is.na(vS3))&(!is.na(vS4))&(sum(!is.na(clusters2[[j]][,"longitude"]))==2))
											{
												indices = which(!is.na(clusters2[[i]][,"longitude"]))
												if (length(indices) == 2)
													{
														curvedarrow(cbind(clusters2[[i]][indices[1],"longitude"],clusters2[[i]][indices[1],"latitude"]),
																	cbind(clusters2[[i]][indices[2],"longitude"],clusters2[[i]][indices[2],"latitude"]),
																	arr.length=0, arr.width=0, lwd=0.2, lty=2, lcol="gray10", arr.col=NA, arr.pos=F,
																	curve=0.1, dr=NA, endhead=F)
													}
											}
									}
								for (k in 1:dim(clusters2[[j]])[1])
									{
										index = which(mutations_data[,"seqName"]==row.names(clusters2[[j]])[k])
										vS5 = raster::extract(pop_liege, cbind(clusters2[[j]][k,c("longitude")],clusters2[[j]][k,c("latitude")]))
										if ((!is.na(vS5))&(!is.na(clusters2[[j]][k,"longitude"])))
											{
												plotTheNode = TRUE
												if ((firstTimePeriod==TRUE)&(clusters2[[j]][k,"collection_date"]>decimal_date(ymd("2020-06-01")))) plotTheNode = FALSE
												if ((secondTimePeriod==TRUE)&(clusters2[[j]][k,"collection_date"]<decimal_date(ymd("2020-08-31")))) plotTheNode = FALSE
												if (plotTheNode)
													{
														index = (((clusters2[[j]][k,"collection_date"]-minYear)/(maxYear-minYear))*100)+1
														points(clusters2[[j]][k,"longitude"], clusters2[[j]][k,"latitude"], pch=16, col=colourScale[index], cex=cexNode)
														points(clusters2[[j]][k,"longitude"], clusters2[[j]][k,"latitude"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
													}
											}
									}
							}
					}
			}
		dev.off()
	}

# 9. Exploring the spatio-temporal distribution of specific mutations

mutations_data = read.csv("Belgian_mutations.csv", head=T)
mutations = c("D614G","S98F","A222V","S477N") # A222V (--> 20A.EU1) and S477N (--> 20A.EU2)
ULG_sequences = mutations_data[which(grepl("ULG-",mutations_data[,"seqName" ])),]
mcc = read.csv("Phylogeographic_runs/All_clades.csv", head=T)
for (i in 1:length(mutations))
	{
		cat(round((sum(grepl(mutations[i],ULG_sequences[,"spikeMutations" ]))/dim(ULG_sequences)[1])*100,1),"\t")
	}
if (showingPlots)
	{
		polIndex1 = which(provinces@data[,"NAME_2"]=="Liège")
		maxArea = 0; polIndex2 = 0; cexNode = 0.7
		for (i in 1:length(provinces@polygons[[polIndex1]]@Polygons))
			{
				if (maxArea < provinces@polygons[[polIndex1]]@Polygons[[i]]@area)
					{
						maxArea = provinces@polygons[[polIndex1]]@Polygons[[i]]@area; polIndex2 = i
					}
			}
		pol = provinces@polygons[[polIndex1]]@Polygons[[polIndex2]]
		p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pol = sps; proj4string(pol) = communes@proj4string
		pop_liege = raster::mask(crop(pop,pol),pol)
		vS1 = raster::extract(pop_liege, mcc[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege, mcc[,c("endLon","endLat")])
		colourScale = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(141)[21:121])
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		sub1 = mcc[which((!is.na(vS1))&(!is.na(vS2))),]		
		pdf("Figure_3_NEW.pdf", width=11, height=8) # dev.new(width=11, height=8)
		par(mfrow=c(2,2), oma=c(0,2,1,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		for (i in 1:length(mutations))
			{
				cols = c(colorRampPalette(brewer.pal(9,"YlGnBu"))(161)[1:101])
				cols = c(colorRampPalette(brewer.pal(9,"Greys"))(201)[1:101])
				plot(pop_liege, col=cols, axes=F, ann=F, box=F, legend=F)
				plot(pol, border="gray30", col=NA, add=T, lwd=0.4)
				selected_sequences = mutations_data[which(grepl(mutations[i],mutations_data[,"spikeMutations"])),"seqName"]
				sub2 = sub1[which(sub1[,"tipLabel"]%in%selected_sequences),]	
				nodesToExplore = sub1[which(sub1[,"node2"]%in%sub2[,"node1"]),"node2"]
				buffer1 = nodesToExplore			
				while (length(buffer1) != 0)
					{
						buffer1 = c(); buffer2 = c()
						for (j in 1:length(nodesToExplore))
							{
								if (sum(sub2[,"node1"]==nodesToExplore[j]) == 2)
									{
										sub2 = rbind(sub2, sub1[which(sub1["node2"]==nodesToExplore[j]),])
										buffer1 = c(buffer1, sub1[which(sub1["node2"]==nodesToExplore[j]),"node1"])
									}	else		{
										buffer2 = c(buffer2, nodesToExplore[j])
									}	
							}
						nodesToExplore = c(buffer1, buffer2)
					}			
				startYears_indices = (((sub2[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
				endYears_indices = (((sub2[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
				startYears_colours = colourScale[startYears_indices]
				endYears_colours = colourScale[endYears_indices]
				for (j in 1:dim(sub2)[1])
					{
						if (sum(sub2[,"node1"]==sub2[j,"node1"]) == 2)
							{
								curvedarrow(cbind(sub2[j,"startLon"],sub2[j,"startLat"]), cbind(sub2[j,"endLon"],sub2[j,"endLat"]), arr.length=0,
						    				arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
						    }
					}
				for (j in 1:dim(sub2)[1])
					{
						if ((!sub2[j,"node1"]%in%sub2[,"node2"])&(sum(sub2[,"node1"]==sub2[j,"node1"])==2))
							{
								points(sub2[j,"startLon"], sub2[j,"startLat"], pch=16, col=startYears_colours[j], cex=cexNode)
								points(sub2[j,"startLon"], sub2[j,"startLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
							}
						points(sub2[j,"endLon"], sub2[j,"endLat"], pch=16, col=endYears_colours[j], cex=cexNode)
						points(sub2[j,"endLon"], sub2[j,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
					}
				for (j in 1:length(clusters2))
					{
						if (sum(!is.na(clusters2[[j]][,"longitude"])) < 3)
							{
								if (length(which(!is.na(clusters2[[j]][,"longitude"]))) > 2)
									{
										sub3 = clusters2[[j]][which(!is.na(clusters2[[j]][,"longitude"])),]
										index1 = which(mutations_data[,"seqName"]==row.names(sub3)[1])
										index2 = which(mutations_data[,"seqName"]==row.names(sub3)[2])
										vS3 = raster::extract(pop_liege, sub3[1,c("longitude","latitude")])
										vS4 = raster::extract(pop_liege, sub3[2,c("longitude","latitude")])
										if ((grepl(mutations[i],mutations_data[index1,"spikeMutations"]))&(grepl(mutations[i],mutations_data[index2,"spikeMutations"]))
											&(!is.na(vS3))&(!is.na(vS4))&(sum(!is.na(clusters2[[j]][,"longitude"]))==2))
											{
												indices = which(!is.na(clusters2[[i]][,"longitude"]))
												if (length(indices) == 2)
													{
														curvedarrow(cbind(clusters2[[i]][indices[1],"longitude"],clusters2[[i]][indices[1],"latitude"]),
																	cbind(clusters2[[i]][indices[2],"longitude"],clusters2[[i]][indices[2],"latitude"]),
																	arr.length=0, arr.width=0, lwd=0.2, lty=2, lcol="gray10", arr.col=NA, arr.pos=F,
																	curve=0.1, dr=NA, endhead=F)
													}
											}
									}
								for (k in 1:dim(clusters2[[j]])[1])
									{
										index = which(mutations_data[,"seqName"]==row.names(clusters2[[j]])[k])
										vS5 = raster::extract(pop_liege, cbind(clusters2[[j]][k,c("longitude")],clusters2[[j]][k,c("latitude")]))
										if ((grepl(mutations[i],mutations_data[index,"spikeMutations"]))&(!is.na(vS5))&(!is.na(clusters2[[j]][k,"longitude"])))
											{
												index = (((clusters2[[j]][k,"collection_date"]-minYear)/(maxYear-minYear))*100)+1
												points(clusters2[[j]][k,"longitude"], clusters2[[j]][k,"latitude"], pch=16, col=colourScale[index], cex=cexNode)
												points(clusters2[[j]][k,"longitude"], clusters2[[j]][k,"latitude"], pch=1, col="gray10", cex=cexNode, lwd=0.4)
											}
									}
							}
					}
			}
		dev.off()
	}
if (showingPlots)
	{
		mcc = read.csv("Phylogeographic_runs/All_clades.csv", head=T); cols1 = list(); cols2 = list()
		cols1[[1]] = rgb(150,150,150,255,maxColorValue=255); cols2[[1]] = rgb(150,150,150,100,maxColorValue=255) # light grey
		cols1[[2]] = rgb(250,165,33,255,maxColorValue=255); cols2[[2]] = rgb(250,165,33,100,maxColorValue=255) # orange
		cols1[[3]] = rgb(222,67,39,255,maxColorValue=255); cols2[[3]] = rgb(222,67,39,100,maxColorValue=255) # red
		cols1[[4]] = rgb(70,118,187,255,maxColorValue=255); cols2[[4]] = rgb(70,118,187,100,maxColorValue=255) # blue
		cols1[[5]] = rgb(76,76,76,255,maxColorValue=255); cols2[[5]] = rgb(60,60,60,100,maxColorValue=255) # dark grey
		cols1[[6]] = rgb(77,77,77,255,maxColorValue=255); cols2[[6]] = rgb(0,0,0,0,maxColorValue=255) # black (and transparent)
		polIndex1 = which(provinces@data[,"NAME_2"]=="Liège")
		maxArea = 0; polIndex2 = 0; cexNode = 0.7
		for (i in 1:length(provinces@polygons[[polIndex1]]@Polygons))
			{
				if (maxArea < provinces@polygons[[polIndex1]]@Polygons[[i]]@area)
					{
						maxArea = provinces@polygons[[polIndex1]]@Polygons[[i]]@area; polIndex2 = i
					}
			}
		pol = provinces@polygons[[polIndex1]]@Polygons[[polIndex2]]
		p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		pol = sps; proj4string(pol) = communes@proj4string
		pop_liege = raster::mask(crop(pop,pol),pol)
		vS1 = raster::extract(pop_liege, mcc[,c("startLon","startLat")])
		vS2 = raster::extract(pop_liege, mcc[,c("endLon","endLat")])
		colourScale = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(141)[21:121])
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		sub1 = mcc[which((!is.na(vS1))&(!is.na(vS2))),]; mutations_years = list()
		MRCAs = mcc[which(!sub1[,"node1"]%in%sub1[,"node2"]),"startYear"]
		for (i in 1:length(mutations))
			{
				selected_sequences = mutations_data[which(grepl(mutations[i],mutations_data[,"spikeMutations"])),"seqName"]
				sub2 = sub1[which(sub1[,"tipLabel"]%in%selected_sequences),]	
				nodesToExplore = sub1[which(sub1[,"node2"]%in%sub2[,"node1"]),"node2"]
				buffer1 = nodesToExplore	; mutation_years = c()		
				while (length(buffer1) != 0)
					{
						buffer1 = c(); buffer2 = c()
						for (j in 1:length(nodesToExplore))
							{
								if (sum(sub2[,"node1"]==nodesToExplore[j]) == 2)
									{
										sub2 = rbind(sub2, sub1[which(sub1["node2"]==nodesToExplore[j]),])
										buffer1 = c(buffer1, sub1[which(sub1["node2"]==nodesToExplore[j]),"node1"])
									}	else		{
										buffer2 = c(buffer2, nodesToExplore[j])
									}	
							}
						nodesToExplore = c(buffer1, buffer2)
					}
				for (j in 1:dim(sub2)[1])
					{
						if ((!sub2[j,"node1"]%in%sub2[,"node2"])&(sum(sub2[,"node1"]==sub2[j,"node1"])==2))
							{
								mutation_years = c(mutation_years, sub2[j,"startYear"])
							}
						mutation_years = c(mutation_years, sub2[j,"endYear"])
					}
				mutations_years[[i]]	 = mutation_years	
			}
		ats = c(decimal_date(ymd("2020-02-01")), decimal_date(ymd("2020-03-01")), decimal_date(ymd("2020-04-01")), decimal_date(ymd("2020-05-01")),
				decimal_date(ymd("2020-06-01")), decimal_date(ymd("2020-07-01")), decimal_date(ymd("2020-08-01")), decimal_date(ymd("2020-09-01")),
				decimal_date(ymd("2020-10-01")), decimal_date(ymd("2020-11-01")), decimal_date(ymd("2020-12-01")))
		minYear = min(mcc[,"startYear"]); maxYear = decimal_date(ymd("2020-12-01"))
		dev.new(width=7, height=4); par(mgp=c(0,0,0), oma=c(1,1,0.5,0.5), mar=c(1.5,1.5,1,1))
		plot(density(mutations_years[[1]]), col=NA, ann=F, axes=F, xlim=c(minYear,maxYear), ylim=c(0,17))
		for (i in 1:length(mutations))
			{
				xx_l = c(density(mutations_years[[i]])$x,rev(density(mutations_years[[i]])$x))
				yy_l = c(rep(0,length(density(mutations_years[[i]])$x)),rev(density(mutations_years[[i]])$y))
				polygon(xx_l, yy_l, col=cols2[[i]], border=0)
			}
		for (i in 1:length(mutations))
			{
				lines(density(mutations_years[[i]])$x, density(mutations_years[[i]])$y, col=cols1[[i]], lwd=0.75)
			}
		# lines(density(MRCAs)$x, density(MRCAs)$y, col=cols1[[5]], lwd=0.75) # introduction events
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=ats, labels=c("01-02","01-03","01-04","01-05","01-06","01-07","01-08","01-09","01-10","01-11","01-12"))
		legend(minYear, 15, mutations, col=unlist(cols2[1:4]), text.col="gray30", pch=16, pt.cex=1.5, box.lty=0, cex=0.7, y.intersp=1.3)
		dev.new(width=7, height=4); par(mgp=c(0,0,0), oma=c(1,1,0.5,0.5), mar=c(1.5,1.5,1,1))
		timeSlices = seq(minYear,maxYear,(maxYear-minYear)/50)
		barplot_tab = matrix(nrow=length(mutations), ncol=length(timeSlices)-1)
		row.names(barplot_tab) = mutations; colnames(barplot_tab) = timeSlices[2:length(timeSlices)]-((maxYear-minYear)/100)
		for (i in 2:length(timeSlices))
			{
				for (j in 1:length(mutations))
					{
						indices = which((mutations_years[[j]]>timeSlices[i-1])&(mutations_years[[j]]<=timeSlices[i]))
						barplot_tab[j,i-1] = length(indices)
					}
			}
		plot(cbind(as.numeric(colnames(barplot_tab)),barplot_tab[1,]), col=NA, ann=F, axes=F, xlim=c(minYear,maxYear))
		for (i in 1:length(mutations))
			{
				xx_l = c(as.numeric(colnames(barplot_tab)),rev(as.numeric(colnames(barplot_tab))))
				yy_l = c(rep(0,length(as.numeric(colnames(barplot_tab)))),rev(barplot_tab[i,]))
				polygon(xx_l, yy_l, col=cols2[[i]], border=0)
			}
		for (i in 1:length(mutations))
			{
				lines(cbind(as.numeric(colnames(barplot_tab)),barplot_tab[i,]), col=cols1[[i]], lwd=0.75)
			}		
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.20,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.2, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=ats, labels=c("01-02","01-03","01-04","01-05","01-06","01-07","01-08","01-09","01-10","01-11","01-12"))
		legend(minYear, 15, mutations, col=unlist(cols2[1:4]), text.col="gray30", pch=16, pt.cex=1.5, box.lty=0, cex=0.7, y.intersp=1.3)
		dev.new(width=7, height=4); par(mgp=c(0,0,0), oma=c(1,1,0.5,0.5), mar=c(1.5,1.5,1,1))
		barplot(barplot_tab, col=c(cols2[[1]],unlist(cols1)[2:4]), border="gray30", axes=F, axisnames=F, ann=F, ylim=c(-10,280))
		selected_dates = c(decimal_date(ymd("2020-03-01")), decimal_date(ymd("2020-04-01")), decimal_date(ymd("2020-05-01")),decimal_date(ymd("2020-06-01")),
						   decimal_date(ymd("2020-07-01")), decimal_date(ymd("2020-08-01")), decimal_date(ymd("2020-09-01")),
						   decimal_date(ymd("2020-10-01")), decimal_date(ymd("2020-11-01")), decimal_date(ymd("2020-12-01")))
		selected_dates_mod = ((selected_dates-minYear)/(maxYear-minYear))*60
		selected_dates_name = c("01-03","01-04","01-05","01-06","01-07","01-08","01-09","01-10","01-11","01-12")
		axis(side=1, lwd.tick=0.4, cex.axis=0.8, mgp=c(0,0.30,0), lwd=0.2, tck=-0.02, col.tick="gray30", col.axis="gray30", col="gray30", at=selected_dates_mod, labels=selected_dates_name)
		axis(side=2, lwd.tick=0.4, cex.axis=0.8, mgp=c(0,0.30,0), lwd=0.2, tck=-0.02, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,300,50), pos=0)
		title(ylab="frequency", cex.lab=1.1, mgp=c(0.5,0,0), col.lab="gray30")
		legend(5, 250, mutations, col=c(cols2[[1]],unlist(cols1)[2:4]), text.col="gray30", pch=16, pt.cex=1.5, box.lty=0, cex=0.7, y.intersp=1.3)
	}

