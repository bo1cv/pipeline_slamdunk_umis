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

# setwd("/mnt/sharc/shared/sudlab1/General/projects/SynthUTR_hepG2_a549/a549_slam/slam/count/")
# setwd("/mnt/sharc/shared/sudlab1/General/projects/SLAMseq_CHO_Mikayla/cho_slam/slam_picr/count/")
# arguments <- data.frame(cpm_cutoff = 1,
#                         bg = 0)


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
  bg_convrate <- bg_convrate[bg_convrate$Coordinate  %in% filter_cpm,]
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
  mod1AIC = AIC1 = 2 + length(times)*log(minm)
  #RMSE <- sqrt(minm/n)

  return (c(thalf_m1=half1_m, AIC_m1=mod1AIC, a_m1 = a_m))
}

model1_fits <- apply(ConvRate_rmbg_norm, 1, fit_model1, times=ts)
model1_fits <- t(model1_fits)


fit_model2 <- function(expressions2, times2){
  moptim2 <- function(x,expression) {
    a = x[1]
    c = x[2]
    mRNA_guess = (1-c)*exp(-a*times2) + c
    diff = mRNA_guess - expression
    sumsqdiff= sum(diff^2)
    return(sumsqdiff)
  }
  outm2 <- optim(c(1,0), fn=moptim2, expression = expressions2, times2)
  a_m2 <- outm2$par[1]
  c_m2 <- outm2$par[2]

  model2_halfm <- function(x, a, c) {
    mRNA_half <- (1-c) * exp(-a*x) + c
    diff = mRNA_half - 0.5
    return (sum(diff^2))
  }
  out2_2m <- optim(1, fn = model2_halfm, a = a_m2, c = c_m2)
  mod2half <- out2_2m$par[1]

  mod2AIC = AIC2 = 4 + length(times2)*log(outm2$value)
  #RMSE <- sqrt(out2_2m$value/n)
  return (c(thalf_m2=mod2half, AIC_m2=mod2AIC, a_m2=a_m2, c_m2=c_m2))
}

model2_fits <- apply(ConvRate_rmbg_norm, 1, fit_model2, times2=ts)
model2_fits <- t(model2_fits)

fit_model3 <- function (times3, expressions3){
  moptim3 <- function(x, expression){
    a = x[1]
    b = x[2]
    c = x[3]
    mRNA_guess = (1-c)*exp(-b*times3) + c*exp(-a*times3)
    diff = mRNA_guess - expression
    sumqdiff2 = sum(diff^2)
    return(sumqdiff2)
  }
  out3m <- optim(c(1,1,0.1), fn=moptim3, expression = expressions3, times3)
  a_m3 <- out3m$par[1]
  b_m3 <- out3m$par[2]
  c_m3 <- out3m$par[3]

  model3_halfm <- function(x, a,b,c){
    mRNA_halfm <- c * exp(-a*x)+(1-c)*exp(-b*x)
    diff = mRNA_halfm - 0.5
    return(sum(diff^2))
  }
  out3_3m <- optim(1, fn = model3_halfm, a = a_m3, b = b_m3, c = c_m3)
  mod3half3 <- out3_3m$par[1]

  mod3AIC = AIC3 = 6 + length(times3)*log(out3m$value)
  return(c(thalf_m3=mod3half3, AIC_m3 = mod3AIC, a_m3 = a_m3, b_m3 = b_m3, c_m3 = c_m3))
}

model3_fits <- apply(ConvRate_rmbg_norm, 1, fit_model3, times3=ts)
model3_fits <- t(model3_fits)

#Models merge
models_aic <- merge(model1_fits, model2_fits, by = "row.names")
rownames(models_aic) <- models_aic[,1]
models_aic[,1] <- NULL
models_aic <- merge(models_aic, model3_fits, by = "row.names")
names(models_aic)[1] <- "id"
write.table(models_aic, file = "models_halflife_decay_aic.tsv", quote = F, row.names = F, sep = "\t")

# Model selection ----------------------------------------------------------

best_half <- matrix(nrow = nrow(models_aic), ncol = 3)
for (i in 1:nrow(models_aic)) {
  best_half[i,] <- c(models_aic[i,"id"], models_aic[i,paste0("thalf_m",which.min(models_aic[i, c("AIC_m1","AIC_m2","AIC_m3")]) )], paste0("model",which.min(models_aic[i, c("AIC_m1","AIC_m2","AIC_m3")]) ))
}
best_half <- data.frame(best_half, stringsAsFactors = F)
best_half <- type.convert(best_half, as.is = T)
colnames(best_half) <- c("id", "halflife", "model")

best_half2 <- separate(best_half, id, sep = "_", into = c("Chromosome","Start", "End", "Name", "Length", "Strand"),remove = FALSE)
best_half2$X1 <- NULL

write.table(best_half, file = "halflife_unfiltered.tsv", quote = F, row.names = F, sep = "\t")

###Create stats of half-lifes
per_model1 = ((sum(best_half2[,"model"] == "model1" ) / nrow(best_half2)))*100
per_model2 = ((sum(best_half2[,"model"] == "model2" ) / nrow(best_half2)))*100
per_model3 = ((sum(best_half2[,"model"] == "model3" ) / nrow(best_half2)))*100
per0_model1 = sum(best_half2[,"model"] == "model1" & best_half2[,"halflife"] < 0)
per0_model2 = sum(best_half2[,"model"] == "model2" & best_half2[,"halflife"] < 0)
per0_model3 = sum(best_half2[,"model"] == "model3" & best_half2[,"halflife"] < 0)
pertmp_model1 = sum(best_half2[,"model"] == "model1" & best_half2[,"halflife"] > max(ts))
pertmp_model2 = sum(best_half2[,"model"] == "model2" & best_half2[,"halflife"] > max(ts))
pertmp_model3 = sum(best_half2[,"model"] == "model3" & best_half2[,"halflife"] > max(ts))
per1_model1 = sum(best_half2[,"model"] == "model1" & best_half2[,"halflife"] < -1)
per1_model2 = sum(best_half2[,"model"] == "model2" & best_half2[,"halflife"] < -1)
per1_model3 = sum(best_half2[,"model"] == "model3" & best_half2[,"halflife"] < -1)
per_stat <- data.frame(model = c("model1","model2","model3"),
                       Percent_model = c(per_model1,per_model2,per_model3),
                       below0 = c(per0_model1,per0_model2,per0_model3),
                       above_tmp = c(pertmp_model1,pertmp_model2,pertmp_model3),
                       below_minus1 = c(per1_model1, per1_model2, per1_model3))
write.csv(per_stat, file = "model_fitting_summary.txt", row.names = F, quote = F)


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
  #Model 2
  model2_fits <- apply(random_data, 1, fit_model2, times2=ts)
  model2_fits <- t(model2_fits)

  #Model 3
  model3_fits <- apply(random_data, 1, fit_model3, times3=ts)
  model3_fits <- t(model3_fits)

  #Models merge
  models_aic <- merge(model1_fits, model2_fits, by = "row.names")
  rownames(models_aic) <- models_aic[,1]
  models_aic[,1] <- NULL
  models_aic <- merge(models_aic, model3_fits, by = "row.names")
  head(models_aic)
  names(models_aic)[1] <- "id"

  best_half <- matrix(nrow = nrow(models_aic), ncol = 3)
  for (i in 1:nrow(models_aic)) {
    best_half[i,] <- c(models_aic[i,"id"], models_aic[i,paste0("thalf_m",which.min(models_aic[i, c("AIC_m1","AIC_m2","AIC_m3")]) )], paste0("model",which.min(models_aic[i, c("AIC_m1","AIC_m2","AIC_m3")]) ))
  }
  best_half <- as.data.frame(best_half)
  rownames(best_half) <- best_half[,1]
  best_half[,1] <- NULL


  return(best_half)
}

#Need ConvRate data as a data frane to generate random data
dataframe_boostrap <- as.data.frame(ConvRate_rmbg_norm)

#Parralel looping
bootstrapped_halflifes <- foreach(i=1:10000, .combine=cbind) %dopar% {
  #Give libraries for job
  library(tidyverse)
  library("matrixStats")

  bootstrap_hl <- bootstrapping_hl(dataframe_ConvRate = dataframe_boostrap,
                                   timepoints = ts)

  bootstrap_hl #cbinding
}

bootstrapped_halflifes <- type.convert(bootstrapped_halflifes, as.is = T)
colnames(bootstrapped_halflifes) <- paste0(c("bs_halflife_", "bs_model_"), rep(seq(1, 1000), each = 2))
head(bootstrapped_halflifes)

#stop parallel
stopCluster(core_loop)


write.table(bootstrapped_halflifes, file = "bootstrap_halflife.tsv", quote = F, sep = "\t")



# Percentile Bootstrap Confidence Intervals  -----------------------------------------
bootstrapped_halflifes <- as.data.frame(bootstrapped_halflifes)
bs_halflife <- bootstrapped_halflifes %>%
  select(contains("id") | contains("bs_halflife"))

#Find  the  0.05  and  0.95  quantile
bs_quantiles = apply(bs_halflife, 1, quantile, c(0.05,0.95))  %>%t()
halflife_PbsCI <- merge(best_half, bs_quantiles, by.x = "id",by.y = 0)

#Add means of bs
halflife_PbsCI$bs_mean <- apply(bs_halflife, 1, mean)
halflife_PbsCI$bs_median <- apply(bs_halflife, 1, median)
halflife_PbsCI$bs_sd <- apply(bs_halflife, 1, sd)
halflife_PbsCI$bs_Cofvar <- (halflife_PbsCI$bs_sd / halflife_PbsCI$bs_mean )*100
names(halflife_PbsCI)[1] <- "Coordinate"

write.table(halflife_PbsCI, "halflife_percentileCIs.tsv", quote = F, sep = "\t")


# Filtering ---------------------------------------------------------------

out_halflife <- halflife_PbsCI[!(halflife_PbsCI$`5%` < halflife_PbsCI$halflife & halflife_PbsCI$`95%` > halflife_PbsCI$halflife) ,] %>% nrow()

f1_halflife <- halflife_PbsCI[(halflife_PbsCI$`5%` < halflife_PbsCI$halflife & halflife_PbsCI$`95%` > halflife_PbsCI$halflife) ,]

above24 <- halflife_PbsCI[halflife_PbsCI$halflife > 24,] %>% nrow()

f2_halflife <- f1_halflife[(f1_halflife$halflife < 24 ) ,]

below0 <- halflife_PbsCI[halflife_PbsCI$halflife < 0,] %>% nrow()

f3_halflife <- f2_halflife[(f2_halflife$halflife < 24 ) ,]

filtering <- rbind(total = nrow(halflife_PbsCI),
                   fileting = (out_halflife + above24 + below0),
                   outofCI = out_halflife,
                   above24h = above24,
                   below0h = below0)
write.table(f3_halflife, file = "halflife_filtered.tsv", quote = F, sep = "\t")
write.table(filtering, file = "halflife_filtered.log", quote = F, sep = "\t", col.names = F)


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
