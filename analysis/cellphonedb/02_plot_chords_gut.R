library(circlize)
setwd('E:/AAA_Labwork/Tcell tissues/v2')
celltype_palette <- c('T Cells' = '#0351A8',
                      'NK' = '#8CB0E0',
                      'ILC' = '#D56D11',
                      'B Cells' = '#FFBB78',
                      'Plasma Cells' = '#234E08',
                      'Monocytes' = '#53CB8B',
                      'Macrophages' = '#D30083',
                      'Dendritic Cells' = '#CB788D',
                      'Endothelial Cells' = '#4E195A',
                      'Telocytes' = '#C58CCF',
                      'Fibroblastic Reticular Cells' = '#AA290F',
                      'Fibroblast' = '#B03FD1',
                      'Smooth Muscle Cells' = '#E8BCCF',
                      'Intestinal Epithelial Cells' = '#64605F',
                      'Enteric Glial Cells' = '#B2AD9A')
mat = read.csv('gut_3_cpdb_chord.csv',row.names = 1)
mat <- as.matrix(mat)
mat[mat<5] = 0
colnames(mat) = row.names(mat)
panel.fun <- function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
  circos.axis(h = "top", labels.cex = 0.5,  sector.index = sector.name, track.index = 2)
}

# set global parameters
circos.par(cell.padding = c(0,0,0,0), gap.degree = 2, start.degree = 90, points.overflow.warning = FALSE)
par(mar=c(5,6,4,1)+.1)

chordDiagram(mat, annotationTrack = "grid",self.link = 1, preAllocateTracks = 1, grid.col = celltype_palette, transparency = 0.3,keep.diagonal =FALSE, symmetric = TRUE)
circos.trackPlotRegion(track.index = 1, panel.fun = panel.fun, bg.border = NA)
dev.copy(jpeg,'vis/gut_3_chord.png', width=14, height=14, units="in", res=500)
dev.off()

mat = read.csv('gut_4_cpdb_chord.csv',row.names = 1)
mat <- as.matrix(mat)
mat[mat<5] = 0
colnames(mat) = row.names(mat)
panel.fun <- function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + 0.1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1)
  circos.axis(h = "top", labels.cex = 0.5,  sector.index = sector.name, track.index = 2)
}

# set global parameters
circos.par(cell.padding = c(0,0,0,0), gap.degree = 2, start.degree = 90, points.overflow.warning = FALSE)
par(mar=c(5,6,4,1)+.1)

chordDiagram(mat, annotationTrack = "grid",self.link = 1, preAllocateTracks = 1, grid.col = celltype_palette, transparency = 0.3,keep.diagonal =FALSE, symmetric = TRUE)
circos.trackPlotRegion(track.index = 1, panel.fun = panel.fun, bg.border = NA)


dev.copy(jpeg,'vis/gut_4_chord.png', width=14, height=14, units="in", res=500)
dev.off()