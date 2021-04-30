## Local ASVs ------------------------------------------
library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)


# This analysis compares trends in elevation locally to trends in latitude globally.
# ASVs present in each habitat were determined using the script ASV_Habitat.R 
# Latitude and elevation ranges are based on ranges from EM_Range.R script.

# This analysis considers only ASVs that occured in at least 1 Waimea sample.
# ASVs were identified as 'local' (occuring only in Waimea) vs. 'non-local' (occuring in both EMP and waimea samples)
# Since this analysis is concerned with latitude, samples in EMP without latitude_deg data were dropped.
# As a results, some ASVs that occured in Waimea and the EMP dataset are not considered due to lack of latitude data.

# Total ASVs that occured in at least 1 waimea samples: 309467
# Of those ASVS, total that also occured in at least 1 EMP sample: 136432

# Total ASVs that occured in at least 1 waimea sample after removing controls etc.: 274176
# Of those ASVs, total that occured in at least 1 EMP samples that had latitude info: 116507

## Identify data files --------------
# ASVs per variable (includes ASVs by habitat, trophic, site name)
var_asvs_file <- "../data/processed/ASVs_by_variable.rds"

# All ASV ranges (latitude globally, elevation locally)
# Waimea marine sites were excluded from this analysis
# ranges_w_marine_file <- "../data/processed/empo_ranges/empo_3_sample_range_correlation.rds"
ranges_no_marine_file <- "../data/processed/empo_ranges/no_marine/empo_3_sample_range_correlation_no_marine.rds"

# All ASV occupancies and abundances (Includes occupancy and abundance in marine sites)
occ_abund_file <- "../data/processed/all_ASV_occupancy_abundance.csv"

# Sample metadata for waimea samples
# sample_metadata_file has the sample data submitted to qiita along with elevation and rainfall
sample_metadata_file <- "../data/sample/mm_16s_hiseqs_full_run_metadata_table.csv"

# Read in ------------------
var_asvs  <- readRDS(var_asvs_file)
empo_3    <- readRDS(ranges_no_marine_file )

occ_abund <-fread(occ_abund_file)
setnames(occ_abund, "otu", "OTU_ID")

wai_meta <- fread(sample_metadata_file)
site_meta <- wai_meta[ , .N, by = .(site_name, elevation)]


## Localized ASVs -------------
# get_local_asvs() identifies which ASVs were found only in Waimea
# looks for columns "coord_min" and "coord_max"
get_local_asvs <-function(range_table, lat_range = c(21.59, 21.65)){
  new_table <- copy(range_table)
  # Set local to logical FALSE
  new_table$local <- F
  # If the ASV has a latitude range greater than Waimea vallye, set local to logical TRUE
  endemic_ASVs  <- new_table[coord_min >= min(lat_range) &
                             coord_max <= max(lat_range),
                             local := T]
  return(new_table)
}

# Identify all local and cosmopolitan ASVs in the data set
# ASVs can be looked up in this table to determine whether they occured exclusively in Waimea
all_local <- get_local_asvs(range_table = empo_3,
                            lat_range = c(21.59, 21.65))

# Add sample occupancy and abundance information
all_local <- merge(all_local,
                   occ_abund,
                   by = "OTU_ID",
                   all.x = T)

fwrite(all_local, "../data/processed/all_ASV_ranges_and_localness.csv")


## Organize ASVS --------------------------------------------------
# make data tables of ASV ranges/group counts for each habitat, trophic, site

# determine ranges for all variable groups
ranges_list <- lapply( names(var_asvs) ,
                       function(var_group) {
  local_by_group <- all_local[OTU_ID %in% var_asvs[[var_group]]]
  local_by_group$group <- var_group
  return(local_by_group)
})

# Convert list to data.table for plotting
ranges_dt <-rbindlist(ranges_list)

# Calculate proportion local vs. total
local_ratios <- lapply(unique(ranges_dt$group),
                       function(var_group) {
  ranges_dt[group %in% var_group & local, .N] /
  ranges_dt[group %in% var_group, .N]
})

names(local_ratios) <- unique(ranges_dt$group)

print("Ratios of local ASVS by habitat, site, trophic")
print(local_ratios)


# create table of ASVs that are not exclusive Waimea
non_local <- all_local[local == F ]

# add log transformed data columns
non_local[ , c("log_occupancy","log_elev_max"):= 
          list( log(occupancy), log(elev_max + 1))]

## Range Boxplots --------------------------------------------------
# Generate boxplots of 
# get_comparisons takes a vector of factor and generates all possible pairwise comparisons 
# returns result list of integer vectors of length 2 that indicate pair positions
get_comparisons <- function(group_vec) {
  comp_mat <- combn(seq(1, length(unique(group_vec))), m = 2)
  comp_list <-
    lapply(seq_len(ncol(comp_mat)), function(i) comp_mat[, i])
  return(comp_list)
}

# Summarize ranges for each metadata variable (habitat, site, trophic)
# Separate out by variable
var_list <- list( 
  habitat   = unique(wai_meta$habitat),
  site_name = unique(wai_meta[order(longitude_deg) ,site_name]),
  trophic    = unique(wai_meta$trophic),
  habitroph = unique(wai_meta[order(habitat,trophic), paste(habitat, trophic)])
)

ranges_out <- lapply(names(var_list), function(var_name) {
  
  # get ASV ranges for a given variable
  var_ranges_dt <- ranges_dt[group %in% var_list[[var_name]]]
  var_ranges_dt[, group := factor(group, levels = var_list[[var_name]])]
  
  fwrite(var_ranges_dt, paste0("../data/processed/",var_name,"_ASV_lat_ranges.csv"))

  # set up comparisons based on the number of unique factors in the variable
  comps <- get_comparisons(var_ranges_dt$group)
  
  # boxplot
  p <-
    ggplot(var_ranges_dt, aes(x = group, y = range_coord, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    #geom_jitter(width = .2) +
    theme_pubr() +
    ylab("ASV Latitudinal Range") +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.title = element_blank()
    ) +
    labs(title = paste0("ASV ranges by ", var_name))
  
  p
  
  print(paste("saving", var_name, "plot"))
  ggsave(
    paste0("../outputs/figures/ASV_ranges_by_", var_name, ".png"),
    plot = p,
    width = 10,
    height = 5
  )
  
  
  #Pull out range and group count mean/sd for each habitat
  range_summary <- var_ranges_dt[, .(
                                      mean_range = mean(range_coord),
                                      sd_range = sd(range_coord),
                                      mean_types = mean(group_count),
                                      sd_types = sd(group_count)
                                    ),
                                    by = group]
  return(list(
    plot = p,
    summary = range_summary,
    ranges = var_ranges_dt
  ))
})

names(ranges_out) <- names(var_list)

# Gradient Site Range Regression ------------------------------
# Generate a regression of ASV ranges by sampling site going up the gradient
# get ASV ranges for a given variable
site_ranges <- ranges_out$site_name$ranges
setnames(site_ranges, "group", "site_name")

# drop marine sites
site_ranges <- site_ranges[site_name %in% wai_meta[habitat != "Marine", site_name] ]
site_ranges[ , site_name := droplevels(site_name)]
site_ranges <- merge(site_ranges, site_meta)


mean_ranges <-site_ranges[ , mean(range_coord), by = elevation]

# model
mod <- lm(range_coord ~ elevation, data = site_ranges)
mod

summary(mod)


# Plot mean ranges
# adding all ranges for each site is too messy, so we just plot the means
# potentially, excluding 'local' ASVs would help and let us make vioin plots

p_el <- ggplot(mean_ranges, aes(x = elevation, y = V1)) +
  geom_point() +
  geom_abline(slope = mod$coefficients[2],
              intercept = mod$coefficients[1]) +
  labs(y = "Latitude Range", caption = "1 degree of latitude for every 124 m elevation") +
  theme_pubr()
  
  
ggsave("../outputs/figures/ASV_lat_range_by_site_elevation_model.png", plot = p_el)


# Gradient ASV range correlation ----------
# generate a regression of ASV elevational ranges vs ASV latitudinal ranges

# plot all data points as violins
# Range elevation vs. Range latitude
ggplot(data = non_local,
       aes(x = as.factor(round(range_elev)), y =range_coord))+
  geom_violin()+
  labs(x = "ASV Elevation Range",
       y = "ASV Abs(Latitude) Range")+
  stat_summary(geom = "point", fun = mean)

ggsave("../outputs/figures/ASV_Range_Latitude_vs_Range_Elevation.png",
       width = 5,
       height = 4)

ggplot(data = non_local,
       aes(x = as.factor(round(elev_max)), y =coord_max))+
  geom_violin()+
  labs(x = "ASV Elevation Max",
       y = "ASV Abs(Latitude) Max")+
  stat_summary(geom = "point", fun = mean)

ggsave("../outputs/figures/ASV_Max_Latitude_vs_Range_Elevation.png",
       width = 5,
       height = 4)



## Plot means as points -----------------
# There are too many points to see any trends in a scatter plot.
# Instead, plot the means at each level (e.g. elevational range value) and add histograms.
# 1. Latitudinal range vs. elevational range = do ASVs distribute across waimea like they do global?
# 2. Max latitude vs. max elevation = are elevation/latitude limiting for ASVs?
# 3. Sample occupancy vs. max elevation = if elevation limits, are tolerant ASVs more common?
# 4. Latitudinal range vs. ASV elevation observed = do waimea sites look like global sites?
# Ranges and maxima are attributes that can be tied to each ASVs
# We also look at ASVs within a site, which repeatedly counts ASVs that occur in multiple.

# function add_histo() adds histograms to a correlation plot
add_histos <- function(main_plot,y_col, x_col, plot_data) {
  
  # set X axis to full data range
  # set plot limits
  main_plot <- main_plot +
    xlim(range(plot_data[[x_col]]))
  
  # define plot limits
  x_min <- min(main_plot$data[[1]])
  x_max <- max(main_plot$data[[1]])
  y_min <- min(main_plot$data[[2]]) # assumes col 3 is standard deviation
  y_max <- max(main_plot$data[[2]])
  
  # Histogram plot for X axis
  xplot <- gghistogram(plot_data, x_col, fill = "gray") +
    clean_theme() +
    scale_y_reverse()
  
  # Histogram plot Y axis  
  # since y axis is dynamic, highlight plotted area for y histogram and inset
  yplot <- gghistogram(plot_data, y_col, fill = "gray", add = "mean")+
    annotate("rect",
             xmin = y_min,
             xmax = y_max,
             ymin = 0,
             ymax = Inf,
             alpha = .2,
             fill = "blue")+
    clean_theme()+
    rotate()
  
  # identify mean so it can be added to main plot
  mean_val <- mean(plot_data[[y_col]])
  
  # define offset for insert placement
  x_offset <- (x_max - x_min) * 0.15
  y_offset <- (y_max - y_min) * 0.60
  # Insert grob inside the scatter plot (x and y swap because of rotated plot)
  yplot_grob <- ggplotGrob(yplot)
  
  ann_plot <- main_plot +
               annotation_custom(grob = yplot_grob,
                                  xmin = x_min,
                                  xmax = x_min + x_offset, 
                                  ymin = y_min + y_offset,
                                  ymax = y_max)+
              geom_line(y = mean_val,
                        linetype = "longdash")   # Add a matching mean line to the main plot

              
  # Arrange alongside the plot
  plot_arr <- ggarrange(
    ann_plot,
    xplot,
    nrow = 2,
    align = "hv",
    heights = c(3,1),
    widths = 3,
    common.legend = TRUE
  )
  plot_arr
  return(plot_arr)
}

# 1. Mean latitudinal range vs elevational range 
# Does elevational range of an ASV predict latitudinal range?

r_lat_by_r_el <- non_local[ , mean(range_coord), by = range_elev]

# linear model
mod.1 <- lm(range_coord ~ range_elev, data = non_local)

# plot
p_r_el <- ggplot(r_lat_by_r_el,
                 aes(x = range_elev,
                     y = V1)) +
                geom_point(size = 2) +
                geom_abline(slope     = mod.1$coefficients[2],
                       intercept = mod.1$coefficients[1],
                       color = "darkblue") +
                labs(y = "Mean(ASV latitudinal range)",
                     x = "ASV elevational range") +
                theme_pubr()

r_el_arr <- add_histos(main_plot = p_r_el,
                      x_col = "range_elev",
                      y_col = "range_coord",
                      plot_data = non_local)
r_el_arr

ggsave("../outputs/figures/ASV_lat_range_by_elev_range.png",
       plot = r_el_arr, width = 5, height = 4)


# 2. Maximum latitude vs maximum elevation 
# Does maximum elevation of an ASV predict maximum latitude?

m_lat_by_m_el <- non_local[ , mean(coord_max), by = elev_max]

# linear model
mod.2 <- lm(coord_max ~ elev_max, data = non_local)

# plot
p_m_el <- ggplot(m_lat_by_m_el,
                 aes(x = elev_max,
                     y = V1)) +
          geom_point(size = 2) +
          geom_abline(slope     = mod.2$coefficients[2],
                      intercept = mod.2$coefficients[1],
                      color = "darkblue")+
          labs(y = "Mean (ASV max. latitude)",
               x= "ASV max. elevation") +
          theme_pubr()

m_el_arr <- add_histos(main_plot = p_m_el,
                      x_col = "elev_max",
                      y_col = "coord_max",
                      plot_data = non_local)
m_el_arr

ggsave("../outputs/figures/Max_ASV_lat_by_max_elev.png", plot = m_el_arr, width = 5, height = 4)


# 3. Sample occupancy vs maximum elevation
# Does max elevation of ASV predict it's spread across samples?
# log transform occupancy for scale and normality
m_occ_by_m_el <- non_local[ , mean(log_occupancy), 
                            by =   elev_max]

# linear model
mod.3 <- lm(log_occupancy ~ elev_max, data = non_local)

#plot
p_occ_el <- ggplot(m_occ_by_m_el, aes(x = elev_max, y = V1)) +
                  geom_point(size = 2) +
                  geom_abline(slope     = mod.3$coefficients[2],
                              intercept = mod.3$coefficients[1],
                              color = "darkblue")+
                  labs(y = "Mean( log( ASV sample occupancy))",
                       x= "ASV max. elevation") +
                  theme_pubr()
p_occ_el

occ_el_arr <- add_histos(main_plot = p_occ_el,
                        x_col = "elev_max",
                        y_col = "log_occupancy",
                        plot_data = non_local)
occ_el_arr



ggsave("../outputs/figures/ASV_lat_range_by_elev_range.png",
       plot = r_el_arr, width = 5, height = 4)

# 4. Latitude range vs elevation occured
# Does the elevation an ASV is observed at predict latitudinal range?
# ASVs by non-marine site
non_local_site <- site_ranges[local == F]

r_lat_by_s_el <- non_local_site[ , mean(range_coord), by = elevation]

# linear model
mod.4 <- lm(range_coord ~ elevation, data = non_local_site)

# plot
p_lat_s_el <- ggplot(r_lat_by_s_el,
                     aes(x = elevation,
                         y = V1)) +
              geom_point(size = 2) +
              geom_abline(slope     = mod.4$coefficients[2],
                          intercept = mod.4$coefficients[1],
                          color = "darkblue")+

              labs(y = "Mean (ASV latitudinal range)", x= "ASV elevation occured") +
              theme_pubr()
p_lat_s_el

lat_s_el_arr <- add_histos(main_plot = p_lat_s_el,
                         x_col = "elevation",
                         y_col = "range_coord",
                         plot_data = non_local_site)
lat_s_el_arr
ggsave("../outputs/figures/ASV_lat_range_by_elev_occured.png",
       plot = lat_s_el_arr, width = 5, height = 4)


# Facet all plots --------------
grobs <-lapply(list( r_el_arr,
                     lat_s_el_arr,
                     m_el_arr,
                     occ_el_arr),
               ggplotGrob)

g <- ggarrange(plotlist=grobs,
               labels = as.list(letters[1:4]),
               nrow = 2,
               ncol = 2)
g
ggsave(paste0("../outputs/figures/ASV_elev_range_facet.pdf"),
       plot = g, 
       width = 8.5,
       height = 11)


## Occupancy (presence in samples) vs. Abundance (total sequences) for ASVs ----------------------

# Occupancy increases with elevation 
ggplot(non_local, aes(as.factor( x =round(elev_max)), y = range_coord))+
  geom_violin()+
  labs(y = "Latitude Range",
       x = "Maximum Elevation Occured (m)")


site_dt$group <- factor(site_dt$group, levels = grad_sites)

# determine number of sites each OTU occurs in
site_occ <- site_dt[ , .(site_count = .N) , by = OTU_ID]

site_dt <- merge(site_dt, site_occ, by = "OTU_ID")

# Crazy plot showing occupancy in global dataset by max elevation by site

ggplot(site_dt, aes(as.factor( x =round(elev_max)), y = log(occupancy), fill = group ))+
  geom_violin(scale = "count")+
  labs(y = "Log(Occupancy)",
       x = "Maximum Elevation Occured (m)")

ggplot( site_dt,
        aes(x = as.factor(x = round(elev_max)),  fill = group))+
  geom_bar()

ggplot( site_dt,
        aes(x = as.factor(x = round(elev_min)),  fill = group))+
  geom_bar()


# Plot site occupancy by site
ggplot(site_dt, aes(x = as.factor(elev_max), y = site_count))+
  geom_violin(scale = "count")


ggplot(non_local, aes(x = occupancy, y = coord_max))+
  geom_point()

