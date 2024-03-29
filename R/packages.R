# ## library() calls go here
# ###load packages and functions
# here::here()
# here::i_am("packages.R")
#   if(!require("pacman")) install.packages("pacman")
#   library(pacman)
#   package.list <- c("conflicted", "dotenv", "drake","data.table","gtools","rlist",
#                     "RCurl","plyr","ggpubr","tidyverse","furrr", "fnmate", "moments","fuzzySim",
#                     "dflow","tictoc","chron","lubridate","httr","TTR",
#                     "grid","gridExtra", "ggridges", "MuMIn", "here",
#                     "viridis", "broom","bbmle","ggthemes", "ggeffects", "betareg",
#                     "igraph","ggraph","magick","cowplot","rriskDistributions",
#                     "rstan", "brms", "tidybayes", "parallel", "hillR", "RInSp","rsample",
#                     "emmeans")
#   # package.list <- c('conflicted', 'dotenv', 'drake', 'abind', 'brms','cowplot',
#   #                   'dplyr', 'egg', 'fnmate', 'furrr', 'ggplot2', 'ggpubr', 'grid', 'gridExtra',
#   #                   'loo', 'lubridate', 'plyr', 'purrr', 'RInSp', 'rlist',
#   #                   'rriskDistributions', 'rstan', 'stringr', 'tibble', 'tidyr', 'vegan',
#   #                   'viridis')
#   p_load(char = package.list, install = FALSE, character.only = TRUE)
#   remotes::install_github("jimjunker1/junkR", upgrade = "never")
#   library(junkR)
#   devtools::install_github("rmcelreath/rethinking", upgrade = "never")
#   conflict_prefer('count', 'dplyr')
#   conflict_prefer('mutate', 'dplyr')
#   conflict_prefer('group_by', 'dplyr')
#   conflict_prefer('traceplot', 'coda')
#   remotes::install_github("milesmcbain/dflow", upgrade = "never")
#   library(dflow)
#   rm("package.list" )
#   
#   source("https://gist.githubusercontent.com/jimjunker1/0ec857c43b1e3a781363c1b1ea7e12ad/raw/4dd2d1078a00500963822d29b2e58ebf39831fb3/geom_flat_violin.R")
#   cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#   load(here::here("data/ocecolors.rda"))
#   overkt_to_C <- function(a){1/(a*(8.61733*10^-5)) - 273.15}
#   theme_mod = function(){theme_bw() %+replace% theme(panel.grid = element_blank())}
#   theme_black = function() {theme_bw() %+replace% theme(panel.background = element_rect(fill = 'transparent', colour = NA),panel.grid = element_blank(), axis.ticks = element_line(color = 'white'),
#                                                           axis.title = element_text(color = 'white'), axis.text = element_text(color = 'white'), plot.background =element_rect(fill = 'transparent', colour = NA),
#                                                           panel.border = element_rect(fill = NA, colour = 'white'))}
#   
#   multiply_prod <- function(x, NPE,...) x/NPE
#   daily_prod <- function(x) x/as.numeric(.$days)
#   
nboot = 1e3
# set.seed(1312)
# boot_samp = paste0("X",sample.int(1e3, size = nboot, replace = FALSE))


'%ni%' <- Negate('%in%')
#   
  # # ! ++++ Plotting aesthetics ++++ ! #
  # oce_temp_disc = c("#E5FA6A","#CF696C","#544685","#072C44","#082A40","#0B222E")#color codes
  oce_temp_pos <- c(256,212,168,124,80,1)#color positions in 'temperature' list of ocecolors
  stream_order <- factor(c("hver", "st6","st9", "st7","oh2","st14"))#stream ordering
  stream_order_list <- stream_order  |>  as.list() |> setNames(object = _,nm = stream_order) #stream ordering for lists
  stream_temp_labels <- c("27.2","17.6","11.2","5.8","5.5","5.0")#stream annual temperature labels
  names(stream_temp_labels) = stream_order #setting named character vector of stream labels
#   
# source("./R/fluxweb_mod_function.R")
# source("./R/Lorenz.R")
# source("./R/Evenness.R")
# quiet <- function(x) { 
#   sink(tempfile()) 
#   on.exit(sink()) 
#   invisible(force(x)) 
# } 
# theme_set(theme_mod())
# options(mc.cores = parallel::detectCores()-1)
# rstan_options(auto_write = TRUE)
# library(knitr)
# library(rmarkdown)

