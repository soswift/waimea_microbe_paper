# this picks up where "rich.r" left off, using rdata file for portability
library(ggplot2)
load("ggdata_trophic_grouped.rdata")


	# i forget what this is for, but it works
	ggdata2 <- ggdata
	ggdata2$Trophic <- "All"
	ggdata2 <- rbind(ggdata, ggdata2)
	ggdata2$Trophic <- factor(ggdata2$Trophic, levels=c("Environmental", "PrimaryProducer", "Consumer", "All"))
	ggdata2$observed <- log10(ggdata2$observed)

	# remove "all" from ggdata2
	ggdata2 <- ggdata2[ggdata2$Trophic != "All", ]
	ggdata2$Trophic <- factor(ggdata2$Trophic, levels=c("Environmental", "PrimaryProducer", "Consumer" ))

	ddg <- 0.6
	p <- ggplot(ggdata2, aes(x=Trophic, y=observed, fill=Habitat))
	p <- p +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	p <- p + geom_violin(trim = TRUE, position = position_dodge(ddg))
	p <- p + geom_point(pch = ".", position = position_jitterdodge(dodge.width=ddg, jitter.width=0.20), size=0.7)
	p <- p + stat_summary(position = position_dodge(ddg), fun.data=function(x){mean_sdl(x, 1)})
	#p <- p + geom_boxplot(width = 0.15, position = position_dodge(ddg), outlier.shape=NA, coef = 0, alpha=0, color="white")
	p <- p + theme(legend.position = c(0.2, 0.15))
	p <- p + theme(legend.background=element_rect(fill=NA, color=NA))
	p <- p + theme(legend.title = element_blank())
	p <- p + theme(aspect.ratio=0.6)
	# order:                             marine     river      terr
	p <- p + scale_fill_manual(values=c("#1f78b4", "#a6cee3", "#fed9a6"))
	p <- p + theme(axis.title.x=element_blank(), axis.text=element_text(size=8))
	p <- p + labs(y="Log10 Observed Richness")
	# lines for means
	trs <- c("Environmental", "PrimaryProducer", "Consumer" )
	tr_cols <- c("#fc8d62",  "#66c2a5", "#756bb1")
	for(i in 1:length(trs)){
		tr_mu <- mean(ggdata2$observed[ggdata2$Trophic == trs[i]])
		tr_sd <- sd(ggdata2$observed[ggdata2$Trophic == trs[i]])
		p <- p + geom_segment(x=i-0.5, xend=i+0.5, y=tr_mu, yend=tr_mu, 
			color=tr_cols[i], size=1.0)
	}


	pdf("log10richness_by_trophic_grouped.pdf")
	p
	dev.off()

# statistics



	# do an anova
	library(car)
	get_eta2 <- function(res){
		res <- as.data.frame(res)
		res$eta2 <- res$"Sum Sq" / sum(res$"Sum Sq")
		return(res)
	}
	m1 <- aov(ggdata2$observed ~ ggdata2$Habitat * ggdata2$Trophic)
	res <- car::Anova(m1)
	res <- get_eta2(res)
	capture.output(res, file="anova_results.txt")
	pdf("anova_varexp_piechart.pdf")
	pie(x=res$eta2, labels=c("Habitat(***)", "Trophic(***)", "Hab*Tro(***)", "Residuals"))
	dev.off()


