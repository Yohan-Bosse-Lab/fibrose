#####################################
###  Version 2023-01-06
# four arguments are required to run the function
# "Exons.txt" Fichier (et path) d'exons avec les positions.
# "SNPs.txt" Fichier (et path) de SNPs avec leur position sur le B37/38.
# figure  Le nom (et path) de la figure à produire
# title: titre (optionnel) de la figure 


#function definition
pwp_syn_nonsyn = function(snps_file = "input/SNPs.txt", exons_file = "input/Exons.txt" , figure = "SERPINA1-PWP.png",title = '') {


snps <- read.table(snps_file, header = T, sep = "\t")
exons <- read.table(exons_file, header = T, sep = "\t")

snps.order <- snps[order(snps$SNP.Position, decreasing = T) , ]

figure_name = sub('.png',paste0('_',Sys.Date(),'.png'),figure)
png(figure_name, height = 3840, width = 5120, pointsize = 75)

par(mar = c(6,1,3,1))
plot(x=0,y=0, xlim = c(max(exons$Exon.Start), min(exons$Exon.End)), ylim = c(0,8), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = title, cex.lab = 1.5)

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
   if(exons$code[i] == 'non codant') rect(exons$Exon.End[i], 3, exons$Exon.Start[i], 0, col = "darkgrey", border = NA)
   if(exons$code[i] == 'codant')     rect(exons$Exon.End[i], 3, exons$Exon.Start[i], 0, col = "black", border = NA)

   #rect(exons$Exon.End[i], 3, exons$Exon.Start[i], 0, col = "black", border = NA)
  
   spacer = 0 
   text((((exons$Exon.End[i]-exons$Exon.Start[i])/2)+exons$Exon.Start[i])+spacer, -0.3, exons$Exon.Name[i], cex = 1.5,xpd =T)
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
  labels =  paste0(snps.order$SNP[i],ifelse(snps.order$Changement.AA[i] == ' ','','\n'),snps.order$Changement.AA[i],'\n',snps.order$N.Patients[i])
	segments(snps.order$SNP.Position[i], 4, (max(exons$Exon.Start) - (i * dif)), 6, lwd =  2)
  segments(max(exons$Exon.Start) - (i * dif), 6, max(exons$Exon.Start) - (i * dif), 7,lwd= 2)
	text(x = max(exons$Exon.Start) - (i * dif), y = 7.1, labels = labels, col = snps.order$SNP.Couleur[i], srt = 90, cex = 0.8, pos = 4,offset=0)
	}

dev.off()

message(paste0('Done preparing figure ',figure,' --- Time is: ',Sys.time()))

}

#Example:
pwp_syn_nonsyn(snps_file = "sandbox/SNPs_ABCA3.txt", exons_file = "sandbox/Exons_ABCA3(003).txt",figure = 'test.png',title = 'this is a super plot')



