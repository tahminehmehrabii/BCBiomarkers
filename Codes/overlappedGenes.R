library(data.table)
library(ggvenn)

setwd("C:/Users/MHR/Desktop/PDAC")

degs <- fread("updown.txt", header = FALSE)$V1
modules <- fread("modulesGenes.txt", header = FALSE)$V1

dir.create("overlappedGenes", showWarnings = FALSE)

png("overlappedGenes/DEGsModules.png", height = 1600, width = 1600, res = 300)
ggvenn(list(DEGs = degs, M1 = modules), fill_color = c("blue", "red"), fill_alpha = 0.3, 
       text_size = 6, set_name_size = 5.5, text_color = "black", stroke_alpha = 0.5, stroke_size = 0.2)
dev.off()

write.table(intersect(degs, modules), "overlappedGenes/DEGsModules.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
