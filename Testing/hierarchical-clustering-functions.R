source("Code//load_count_data.R")
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(bioDist)
library(edgeR)
library(stats)
library(gplots)

# set colors
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)


#Function for creating and plotting heatmap
create.heatmap = function(data, type = "samples", clustering.method = "complete", 
                          distance.measure = "euclidian", label_with = "ad"){
  # set color scheme for heatmap
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  if(type == "samples"){
    data = t(data)
  }
  if(distance.measure == "cor"){
    dist = cor.dist(data, abs = TRUE)
  }
  else{
    dist = dist(data, method = distance.measure)
  }
  
  dist.matrix = as.matrix(dist)
  if(type == "samples"){
    colnames(dist.matrix) = groups_ad
    if( label_with == "probiotic"){
      colnames(dist.matrix) = groups_p
    }
    if( label_with == "maternal atopy"){
      colnames(dist.matrix) = groups_mad
    }
  }
  else{ 
    colnames(dist.matrix) = NULL
  }
  rownames(dist.matrix) = NULL
  heatmap = pheatmap(dist.matrix,
                     clustering_distance_rows = dist,
                     clustering_distance_cols = dist,
                     col=colors, scale = "none",  clustering_method = clustering.method)
  plot(heatmap)
}

create.heatmap(log.cpm, type = "miRNA")