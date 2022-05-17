library(stringr)
library(tidyverse)
library("optparse")
library("matrixStats")
library(foreach)
library(doParallel)

### Parse input options
option_list = list(
            make_option(c("-d", "--count-directory"),
                        type="character",
                        dest = "count_dir",
                        help="full path directory of tsv count files generated
                        by slamdunk"),
            make_option(c("-t", "--cpm--threshold"),
                        type="double",
                        dest = "cpm_cutoff",
                        help="expression threshold, CPM cutoff"),
            make_option(c("-b", "--background"),
                        type="integer",
                        dest = "bg",
                        help="user specify if ctrl samples without 4SU are
                        present for analysis (0: no, 1:yes)"),
            make_option(c("-o", "--output-directory"),
                        type="character",
                        dest = "out_dir",
                        help="full path of output directory")
                        )

arguments <- parse_args(OptionParser(option_list = option_list))

setwd(arguments$count_dir)

# setwd("/mnt/sharc/fastdata/bo1cv/a549_slam/slam/count")
# setwd("/mnt/sharc/shared/sudlab1/General/projects/SLAMseq_CHO_Mikayla/cho_slam/slam_picr/count/")
# arguments <- data.frame(cpm_cutoff = 1,
#                         bg = 1)


# Create Filter ConversionRate table by CPM cut off -----------------------------------------------------------------
###Load data and apply CPM cut off
data_names_files = list.files(pattern="*h-R*.*_tcount.tsv")
data_names <- str_remove_all(data_names_files, pattern = "trimmed-")
data_names <- str_remove(data_names, pattern = "_processed_slamdunk_mapped_filtered_dedup_tcount.tsv")
data_samples <- lapply(data_names_files, read.table, header = T, colClasses=c(rep("character", 6), rep("numeric",10)), skip = 2)
names(data_samples) <- data_names
damples_ConvRate <- lapply(data_samples, subset, select = c("Chromosome","Start","End","Name","Length","Strand", "ConversionRate")) %>%
  purrr::reduce(full_join, by = c("Chromosome","Start","End","Name","Length","Strand"))
names(damples_ConvRate) <- c("Chromosome","Start","End","Name","Length","Strand",data_names)
damples_ConvRate %>%  unite("Coordinate", 1:6, sep = "_") -> damples_ConvRate

filter_data <- lapply(data_samples, subset, select = c("Chromosome","Start","End","Name","Length","Strand", "ReadsCPM")) %>%
  purrr::reduce(full_join, by = c("Chromosome","Start","End","Name","Length","Strand"))
names(filter_data) <- c("Chromosome","Start","End","Name","Length","Strand",data_names)
filter_data %>%  unite("Coordinate", 1:6, sep = "_") -> filter_data
row.names(filter_data) <- filter_data$Coordinate
filter_data$Coordinate <- NULL

###Keep only one with at least two rep with CPM > 1
replicates = str_extract(colnames(filter_data), regex("(?<=-).+(?=h)")) %>% as.numeric() %>% unique()
filter_cpm <- vector()
for (y in 1:length(replicates)) {
  col_vec = select(filter_data, contains(paste0(replicates[y], "h")) ) %>% colnames()

  for (i in 1:nrow(filter_data)) {
    if ( (sum(filter_data[i,col_vec] > 1) < 2) ) {
      filter_cpm <- append(filter_cpm, row.names(filter_data)[i] ) %>% unique()
    } } }

filter_data_ConvRate <- filter(damples_ConvRate, !(Coordinate %in% filter_cpm))

head(filter_data_ConvRate)
row.names(filter_data_ConvRate) <- filter_data_ConvRate$Coordinate
filter_data_ConvRate$Coordinate <- NULL
head(filter_data_ConvRate)

############"

if (arguments$bg == 1) {
  data_bg = list.files(pattern= regex("*no4su-R..*_tcount.tsv"))
  bg_convrate <- lapply(data_bg,read.table, header = T, colClasses=c(rep("character", 6), rep("numeric",10)), skip = 2)
  bg_convrate <- lapply(bg_convrate, subset, select = c("Chromosome","Start","End","Name","Length","Strand", "ConversionRate")) %>%
    purrr::reduce(full_join, by = c("Chromosome","Start","End","Name","Length","Strand"))
  bg_convrate %>% unite("Coordinate",1:6,sep = "_") -> bg_convrate
  bg_convrate <- bg_convrate[!(bg_convrate$Coordinate  %in% filter_cpm),]
  row.names(bg_convrate) <- bg_convrate$Coordinate
  bg_convrate$Coordinate <- NULL
  if (ncol(bg_convrate) > 2) {avg_bg <- apply(bg_convrate[,2:ncol(bg_convrate)], 1, mean)
  filter_data_ConvRate <- sweep(filter_data_ConvRate, 1, avg_bg, '-')
  } else {filter_data_ConvRate <- sweep(filter_data_ConvRate, 1, bg_convrate[, "ConversionRate"], '-') }
}

#Normalize on T0
avg_0h <- dplyr::select(filter_data_ConvRate, matches("-0h-")) %>% apply(1, mean)
ConvRate_rmbg_norm <- apply(filter_data_ConvRate, 2,function(x) x/avg_0h)
ConvRate_rmbg_norm <- na.omit(ConvRate_rmbg_norm)

setwd(arguments$out_dir)
write.table(ConvRate_rmbg_norm, file = "ConvRate_processed.tsv", quote = F,
row.names = T, sep = "\t")

#Timepoints
ts = str_extract(colnames(ConvRate_rmbg_norm), regex("(?<=-).+(?=h)")) %>% as.numeric()



# Model fitting -----------------------------------------------------------
fit_model1 <- function (expressions, times) {
  moptim1 <- function(a, expression){
    mRNA_guess = exp(-a * times)
    err = expression - mRNA_guess
    sq_err = err^2
    sumsq = sum(sq_err)
    return (sumsq)
  }

  outm <- optim(1, moptim1, expression=expressions)
  minm <- outm$value
  a_m <- outm$par[1]
  half1_m <- log(2)/a_m
  mod1half = half1_m
  #RMSE <- sqrt(minm/n)
  return (c(halflife=half1_m, decayrate = a_m))
}

model1_fits <- apply(ConvRate_rmbg_norm, 1, fit_model1, times=ts)
model1_fits <- t(model1_fits)

write.table(model1_fits, file = "halflife_unfiltered.tsv", row.names = T, quote = F, sep = "\t")


# Bootstrap Conf intervals ------------------------------------------------
#Set core number parralel
cores=detectCores()
core_loop <- makeCluster(cores[1]-1)

#Start parralel
registerDoParallel(core_loop)

#Function for to loop over
bootstrapping_hl <- function(dataframe_ConvRate, timepoints) {
  tmp_bs <- timepoints %>% unique()
  random_data <- matrix(nrow = nrow(dataframe_ConvRate), ncol = 0, byrow = T)
  rownames(random_data) <- rownames(dataframe_ConvRate)
  for (i in tmp_bs) {
    #Getting stats on data separate timepoints
    std <- dplyr::select(dataframe_ConvRate, matches(paste0("-",i,"h-"))) %>% apply(1,sd)
    mean <- dplyr::select(dataframe_ConvRate, matches(paste0("-",i,"h-"))) %>% apply(1,mean)
    col_num = ncol(dplyr::select(dataframe_ConvRate, matches(paste0("-",i,"h-"))))
    tmp_stats <- data.frame(rep(col_num, length(std)),mean, std)
    #Generating random data
    random_tmp <- apply(tmp_stats, 1, function(x) rnorm(x[1], mean=x[2], sd=x[3])) %>% t() %>% as.data.frame()
    colnames(random_tmp) <- paste0(rep(i, col_num), "h")
    random_data <- merge(random_data, random_tmp, by = 0)
    rownames(random_data) <- random_data$Row.names
    random_data$Row.names <- NULL
  }
  random_data <- as.matrix.data.frame(random_data)

  #Model 1
  model1_fits <- apply(random_data, 1, fit_model1, times=ts)
  model1_fits <- t(model1_fits)
  return(model1_fits)
}

#Need ConvRate data as a data frane to generate random data
dataframe_boostrap <- as.data.frame(ConvRate_rmbg_norm)

#Parralel looping
bootstrapped_halflifes <- foreach(i=1:1000, .combine=cbind) %dopar% {
  #Give libraries for job
  library(tidyverse)
  library("matrixStats")

  bootstrap_hl <- bootstrapping_hl(dataframe_ConvRate = dataframe_boostrap,
                                   timepoints = ts)

  bootstrap_hl #cbinding
}

head(bootstrapped_halflifes)
colnames(bootstrapped_halflifes) <- paste0(c("bs_halflife_", "bs_decayrate_"), rep(seq(1, 1000), each = 2))

#stop parallel
stopCluster(core_loop)


write.table(bootstrapped_halflifes, file = "bootstrap_halflife.tsv", quote = F, sep = "\t", row.names = T)

# Calculate Percentile and Empirical/Parametric Bootstrap Confidence Intervals  -----------------------------------------
bootstrapped_halflifes <- as.data.frame(bootstrapped_halflifes)
bs_halflife <- bootstrapped_halflifes %>%
  select(contains("id") | contains("bs_halflife"))

### Percentile Bootstrap Confidence Intervals
#Find  the  0.05  and  0.95  quantile
bs_quantiles = apply(bs_halflife, 1, quantile, c(0.05,0.95))  %>%t()
halflife_bsCI <- merge(model1_fits, bs_quantiles, by = 0)
colnames(halflife_bsCI)[4] <- "Percentile5%"
colnames(halflife_bsCI)[5]  = "Percentile95%"


### Empirical Bootstrap Confidence Intervals
#deltastar = bsmean - samplemean
deltastar <- as.data.frame(model1_fits) %>% select(halflife) %>% merge(bs_halflife, by = 0)
deltastar[,3:ncol(deltastar)] <- sweep(deltastar[,3:ncol(deltastar)], 1, deltastar[, "halflife"] , '-')
row.names(deltastar) <- deltastar$Row.names
deltastar$Row.names <- NULL
deltastar$halflife <- NULL

#Find  the  0.05  and  0.95  quantile  for  deltastar
dd = apply(deltastar, 1, quantile, c(0.05,0.95))  %>%t()

# samplemen - deltastar.05, samplemean - deltastar.95
halflife_bsCI$`Param5%` <-  halflife_bsCI[,"halflife"]  - dd[,2]
halflife_bsCI$`Param95%` <-  halflife_bsCI[,"halflife"]  - dd[,1]

#Add sats of bs
halflife_bsCI$bs_mean <- apply(bs_halflife, 1, mean)
halflife_bsCI$bs_median <- apply(bs_halflife, 1, median)
halflife_bsCI$bs_sd <- apply(bs_halflife, 1, sd)
halflife_bsCI$bs_Cofvar <- (halflife_bsCI$bs_sd / halflife_bsCI$bs_mean )*100
names(halflife_bsCI)[1] <- "Coordinate"

write.table(halflife_bsCI, "halflife_percentileCIs.tsv", quote = F, sep = "\t", row.names = F)

# Filtering ---------------------------------------------------------------

out_percentileCI <- halflife_bsCI[!(halflife_bsCI$`Percentile5%` < halflife_bsCI$halflife & halflife_bsCI$`Percentile95%` > halflife_bsCI$halflife) ,]
out_paramCI <- halflife_bsCI[!(halflife_bsCI$`Param5%` < halflife_bsCI$halflife & halflife_bsCI$`Param95%` > halflife_bsCI$halflife) ,]

shared <- length(out_percentileCI$Coordinate %in% out_paramCI$Coordinate)

above24 <- halflife_bsCI[halflife_bsCI$halflife > 24,] %>% nrow()

below0 <- halflife_bsCI[halflife_bsCI$halflife < 0,] %>% nrow()

f1_halflife <- halflife_bsCI[(halflife_bsCI$halflife < 24 & halflife_bsCI$halflife > 0) ,]

out_percentileCI_ok <- f1_halflife[!(f1_halflife$`Percentile5%` < f1_halflife$halflife & f1_halflife$`Percentile95%` > f1_halflife$halflife) ,]
out_paramCI_ok <- f1_halflife[!(f1_halflife$`Param5%` < f1_halflife$halflife & f1_halflife$`Param95%` > f1_halflife$halflife) ,]
out <- nrow(out_percentileCI_ok) + nrow(out_paramCI_ok)

f2_halflife <- f1_halflife[!(f1_halflife$Coordinate %in% out_paramCI$Coordinate) & !(f1_halflife$Coordinate %in% out_percentileCI$Coordinate)  ,]

filtering <- rbind(total = nrow(halflife_bsCI),
                   filtering = nrow(halflife_bsCI) - nrow(f2_halflife),
                   outofPercentileCI = nrow(out_percentileCI),
                   outofParamCI = nrow(out_paramCI),
                   shared_between_CIs_methods = shared,
                   above24h = above24,
                   below0h = below0)
f2_halflife <-f2_halflife %>% separate(Coordinate, into = c("chr", "start", "end", "transcript_id", "length", "strand"), sep = "_", convert = T)
write.table(f2_halflife, file = "halflife_filtered.tsv", quote = F, sep = "\t")
write.table(filtering, file = "halflife_filtered.log", quote = F, sep = "\t", col.names = F, row.names = T)


# STREME ------------------------------------------------------------------

#
# ###Generate bed for streme
# #Order them and select for bed file
# best_half <- best_half[order(best_half$halflife),1:6]
# #Take the 10% most stable and less stable
# howmany <- (nrow(best_half) * 0.1) %>% round()
# most_stable <- head(best_half, n = howmany)
# less_stable <- tail(best_half, n = howmany)
# write.table(most_stable, file = "most_stable.bed", quote = F, row.names = F, sep = "\t", col.names = F)
# write.table(less_stable, file = "less_stable.bed", quote = F, row.names = F, sep = "\t", col.names = F)
