nestdata <- readRDS("H2NODF.RDS")
library(gamlss)


# functions to summarize site
zscore <- function(emp, sim){ (emp - mean(sim))/sd(sim) }

pval_norm <- function(emp, sim){
	normfit <- gamlss(sim~1, family="NO")
	mu_est <- normfit$mu.coefficients[1]
	sd_est <- exp(normfit$sigma.coefficients[1])
	pnorm(q=emp, mean=mu_est, sd=sd_est, lower.tail=emp <= mu_est)
}

sum_site <- function(x){
	# compute Z-score and P-value for H2
	H2_Z <- zscore(emp=x$H2_emp, sim=x$H2_sim)
	H2_P <- pval_norm(emp=x$H2_emp, sim=x$H2_sim)
	# same for NODF
	NODF_Z <- zscore(emp=x$NODF_emp, sim=x$NODF_sim)
	NODF_P <- pval_norm(emp=x$NODF_emp, sim=x$NODF_sim)
	# output a row of a data.frame
	return(data.frame(
		site=x$site, habitat=x$habitat,
		filename=x$filename,
		H2_Zscore=H2_Z, H2_Pval=H2_P,
		H2_emp=x$H2_emp, H2_sim_mean=mean(x$H2_sim), h2_sim_sd=sd(x$H2_sim),
		NODF_Zscore=NODF_Z, NODF_Pval=NODF_P,
		NODF_emp=x$NODF_emp, NODF_sim_mean=mean(x$NODF_sim), NODF_sim_sd=sd(x$NODF_sim)

	))
}

output_rows <- lapply(X=nestdata, FUN=sum_site)

output_df <- do.call("rbind", output_rows)

write.table(file="results.csv", x=output_df, sep=",", quote=F, row.names=F)
saveRDS(output_df, file="results.RDS")

