#####################################
###  Version 2023-01-06
#three required arguments to the function
# "Exons.txt" Avoir la liste des exons avec les positions.
# "SNPs.txt" Avoir la liste des SNPs avec leur position sur le B37/38.
# figure  Le nom de la figure


#function 
pwp_syn_nonsyn = function(snps_file = "input/SNPs.txt", exons_file = "input/Exons.txt" , figure = "SERPINA1-PWP.png") {


snps <- read.table(snps_file, header = T, sep = "\t")
exons <- read.table(exons_file, header = T, sep = "\t")

snps.order <- snps[order(snps$SNP.Position, decreasing = T) , ]

figure_name = sub('.png',paste0('_',Sys.Date(),'.png'),figure)
png(figure_name, height = 3840, width = 5120, pointsize = 75)

par(mar = c(6,0,0,0))
plot(x=0,y=0, xlim = c(max(exons$Exon.Start), min(exons$Exon.End)), ylim = c(0,15), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", cex.lab = 1.5)

box(col = "white", lwd = 10)

## Ligne de base
rect(min(exons$Exon.End), 1, max(exons$Exon.Start), 2, col = "black", border = NA)


## Axis 1
positions = seq(to = max(exons$Exon.Start)+100,from= min(exons$Exon.End)-100,length.out=15)
axis(1, at = positions, labels = round(positions/1000), cex.axis = 1.1, lwd = 15, line = 1)

mtext("Chromosome 14 (kb)", 1, cex = 1.5, line = 4)


## Exons #darkgrey = non codant
for(i in seq_along(exons$Exon.Start))
{
   #if(i ==1) rect(exons$Exon.End[i], 3, exons$Exon.Start[i], 0, col = "darkgrey", border = NA)
   rect(exons$Exon.End[i], 3, exons$Exon.Start[i], 0, col = "black", border = NA)

   spacer = 0 
   text((((exons$Exon.End[i]-exons$Exon.Start[i])/2)+exons$Exon.Start[i])+spacer, -0.3, exons$Exon.Name[i], cex = 1.5)
}

## SNPs-Position sur ligne de base
for (i in 1:nrow(snps.order))
{
	segments(snps.order$SNP.Position[i]-4, 2, snps.order$SNP.Position[i]+4, 4, col = "black",lwd =2 )
}



## SNP rs
dif = ((max(exons$Exon.Start) -  min(exons$Exon.End)) / nrow(snps)) # 

for (i in 1:nrow(snps.order))
{
	segments(snps.order$SNP.Position[i], 4, (max(exons$Exon.Start) - (i * dif)), 6, lwd =  2)
  segments(max(exons$Exon.Start) - (i * dif), 6, max(exons$Exon.Start) - (i * dif), 7,lwd= 2)
	text(x = max(exons$Exon.Start) - (i * dif), y = 7.1, labels = snps.order$SNP[i], col = snps.order$SNP.Couleur[i], srt = 90, cex = 0.8, pos = 4,offset=0)
	}

dev.off()

message(paste0('Done preparing figure ',figure,' --- Time is: ',Sys.time()))

}

#Run the function below...
pwp_syn_nonsyn(snps_file = "C:/Users/renseb01/Desktop/SNPs_ABCA3.txt", exons_file = "C:/Users/renseb01//Desktop/Exons_ABCA3.txt",figure = 'test.png')



