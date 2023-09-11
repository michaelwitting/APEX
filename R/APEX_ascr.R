library(igraph)
library(dplyr)
library(readr)
library(Spectra)
library(MsBackendMgf)
library(MetaboCoreUtils)


### define variables ----------------------------------------------------------
register(bpstart(SnowParam(2, progressbar = TRUE)))


# get working directory
wd <- getwd()

# set path to desired dataset
path <- paste0(wd, "/data/MSV000086293/")


# merged network file
net_fl <- paste0(path, 
                 "/merged_networks_ascr.tsv")

# ids file (manual annotations)
ids_fl <- paste0(path, 
                 "/ids.csv")

ids <- read.csv(ids_fl)[,1:7]

# align column names 

colnames(ids) <- c("Row","m.z", "RT"  ,"Annotation.ID","formula" ,"adduct" , "evidence")

# ms2 files
ms2_fl <- paste0(path, 
                 "ms2_spectra_combined.mgf")

# optional cosine filter 
# can be used if a higher threshold should be used for APEX if spectral 
# similarity edges to non-GPNAE species arise
cosine_filt <- 0.7



# get all manual ascr species that will be used as seed node set (SNS)
ascr <- ids[grepl("ascr", ids$Annotation.ID) |
              grepl("icas", ids$Annotation.ID),]$Row



# 1) Load merged networks ------------------------------------------------------
net_df <- read_tsv(net_fl)

# separate the network in order to sort the `from` and `to` column
#
# order nodes:
## homologous series
homol <- net_df[net_df$interaction == "HomolSeries", 
                c("from", "to", "HSIDs","mzincrement", "RTincrement")]

homol$from2 <- with(homol, pmin(from, to))
homol$to2 <- with(homol, pmax(from, to))
homol$from <- homol$from2
homol$to <- homol$to2
homol <- homol %>% 
  select(c("from", "to", "HSIDs","mzincrement", "RTincrement"))

## spectral similarity
spec <- net_df[net_df$interaction == "Similarity", 
               c("from", "to", "gnps_cosine", "gnps_mass")]

spec$from2 <- with(spec, pmin(from, to))
spec$to2 <- with(spec, pmax(from, to))
spec$from <- spec$from2
spec$to <- spec$to2
spec <- spec %>% select(c(c("from", "to", "gnps_cosine", "gnps_mass")))

## mass-difference
diff <- net_df[net_df$interaction == "MassDifference", 
               c("from", "to", "group_mass", "formula_mass", "mass_mass")]

diff$from2 <- with(diff, pmin(from, to))
diff$to2 <- with(diff, pmax(from, to))
diff$from <- diff$from2
diff$to <- diff$to2
diff <- diff %>% 
  select(c("from", "to", "group_mass", "formula_mass", "mass_mass"))


# merge again to have an ordered network
net_df <- merge(spec, diff, by = c("from", "to"), all = TRUE)
net_df <- merge(net_df, homol, by = c("from", "to"), all = TRUE)

write_csv(net_df, 
          file = paste0(path, "results/net_df.csv"))

# 2) Load corresponding MS2 data -----------------------------------------------

ms2_spectra <- Spectra(ms2_fl,
                       source = MsBackendMgf(),
                       backend = MsBackendDataFrame()
)


# 3) Align Nodes/Edges ---------------------------------------------------------

## Remove Features that are not in node table 'ids'
net_df <- net_df[which(net_df$from %in% ids$Row &
                         net_df$to %in% ids$Row), ]

## Remove MS2 spectra that are not in node table 'ids'
ms2_spectra <- ms2_spectra[which(ms2_spectra$CLUSTER_ID %in% ids$Row)]

## add node attribute if MS2 is available
ids$ms2_available <- ifelse(ids$Row %in% unique(ms2_spectra$CLUSTER_ID),
                            "Yes", "No"
)


## Filtering the network if needed 
net_df <- net_df |> filter(is.na(gnps_cosine) | gnps_cosine >= cosine_filt, )


# nodes that are involved in homol series: 
homolNodes <- unique(c(homol$from, homol$to))
ids$homol <- ifelse(ids$Row %in% homolNodes, "yes", "no")

# change to character in order to later use `strsplit` in the APEX workflow
net_df$mzincrement <- as.character(net_df$mzincrement)
net_df$mass_mass <- as.character(net_df$mass_mass)


# 4) Create igraph network -------------------------------------------------
net <- graph_from_data_frame(d = net_df, vertices = ids, directed = TRUE)


# 5) Find all neighbors of ascrs ----------------------------------------------
## graph of all neighbors of ascrs
ascr_net <- induced_subgraph(net, neighbors(net, ascr, mode = "total"))


# 6) Annotation ----------------------------------------------------------------


layers <- list(
  "1_layer_homol" = c("HSIDs"),
  "1_layer_mass" = c("mass_mass"),
  "1_layer_gnps" = c("gnps_cosine"),
  "2_layers_mass_homol" = c("mass_mass", "HSIDs"),
  "2_layers_gnps_homol" = c("gnps_cosine", "HSIDs"),
  "2_layers_gnps_mass" = c("gnps_cosine", "mass_mass"),
  "3_layers" = c("gnps_cosine", "mass_mass", "HSIDs")
)


## Make a new column/ node attribute which is called "network_annotation" and
## "network_comment"
ids$network_annotation <- NA
ids$network_formula <- NA
ids$network_inputNode <- NA
ids$network_layer <- NA
ids$network_comment <- NA
ids$key_fragments <- NA
ids$mz_diff <- NA

results <- data.frame(matrix(ncol = 16, nrow = 0))
colnames(results) <- colnames(ids)

#i = 8
for (i in 1:length(ascr)) {
  df_i <- net_df[(net_df$from == ascr[i] | net_df$to == ascr[i]), ]
  #df_i <- net_df[c(net_df$from == ascr[i], net_df$to == ascr[i]), ]
 # df_i <- net_df[(net_df$from == "Cluster_2815"| net_df$to == "Cluster_2815"), ]
  
  
  
  for (s in 1:nrow(df_i)) {
    sub <- df_i[s, ]
    
    ### check if ascr is in from or to ###
    col <- ifelse(sub$to == ascr[i], "from", "to")
    row <- which(ids$Row == sub[[col]])
    row_ascr <- which(ids$Row == ascr[i])
    
    ascr_mzdiff <- sign(ids[row_ascr, ]$m.z -
                           ids[row, ]$m.z)
    ## test if substraction of mass diff of seed node happens
    ## and if yes test if mass difference formula is contained in seed sum formula
    if (ascr_mzdiff == 1 & !is.na(sub$formula_mass)){
      if (!containsElements(ids[row_ascr, ]$formula, 
                            sub$formula_mass)) { ## if not contained it will not be annotated
        next
      }} else {

        
        if (any(complete.cases(sub[, layers[["3_layers"]]]))) {
          
          ids[row, ]$network_annotation <-
            paste0(
              ids[row_ascr, ]$Annotation.ID,
              ifelse(ascr_mzdiff == -1, "+", "-"), sub$formula_mass
            )
          ## layers that have been used
          ids[row, ]$network_layer <- "3_layers"
          
          
        } else if (any(complete.cases(sub[, layers[["2_layers_gnps_mass"]]]))) {
          ids[row, ]$network_annotation <-
            paste0(ids[row_ascr, ]$Annotation.ID,
                   ifelse(ascr_mzdiff == -1, "+", "-"), sub$formula_mass)
          
          ## layers that have been used
          ids[row, ]$network_layer <- "2_layers_gnps_mass"
          
          
        } else if (any(complete.cases(sub[, layers[["2_layers_gnps_homol"]]]))) {
          ids[row, ]$network_annotation <-
            paste0(ids[row_ascr, ]$Annotation.ID,
                   ifelse(ascr_mzdiff == -1, "+", "-"), sub$gnps_mass)
          
          ## layers that have been used
          ids[row, ]$network_layer <- "2_layers_gnps_homol"
          
          
        } else if (any(complete.cases(sub[, layers[["1_layer_gnps"]]]))) {
          ids[row, ]$network_annotation <-
            paste0(ids[row_ascr, ]$Annotation.ID,
                   ifelse(ascr_mzdiff == -1, "+", "-"), sub$gnps_mass)
          ## layers that have been used
          ids[row, ]$network_layer <- "1_layer_gnps"
          
          
        }  else if (any(complete.cases(sub[, layers[["2_layers_mass_homol"]]]))) {


              ids[row, ]$network_annotation <-
                paste0(
                  ids[row_ascr, ]$Annotation.ID,
                  ifelse(ascr_mzdiff == -1, "+", "-"), sub$formula_mass
                )
              ## layers that have been used
              ids[row, ]$network_layer <- "2_layers_mass_homol"
              
       
        }else if (any(complete.cases(sub[, layers[["1_layer_mass"]]]))) {

              ids[row, ]$network_annotation <-
                paste0(
                  ids[row_ascr, ]$Annotation.ID,
                  ifelse(ascr_mzdiff == -1, "+", "-"), sub$formula_mass
                )
              ## layers that have been used
              ids[row, ]$network_layer <- "1_layer_mass"
              
          
          
          
          
          
        } else if (any(complete.cases(sub[, layers[["1_layer_homol"]]]))) {

            ids[row, ]$network_annotation <-
              paste0(
                ids[row_ascr, ]$Annotation.ID,
                ifelse(ascr_mzdiff == -1, "+", "-"), sub$formula_mass
              )
            
          
            ## layers that have been used
            ids[row, ]$network_layer <- "1_layer_homol"
            
          
          
          
        }else
          next
        
        ## input node information
        ids[row, ]$network_inputNode <- paste0(ascr[i]) 
        
        ids[row, ]$mz_diff <- ids[row_ascr, ]$m.z -
          ids[row, ]$m.z
        
        ## "network_comment" = values of layers
        ids[row, ]$network_comment <-
          paste0(
            "simil:", sub$gnps_cosine,
            "; mzdiff:", sub$group_mass,
            "; homol mz:", sub$mzincrement, "; homol RT:", sub$RTincrement
          )
        
        if (!grepl("/",sub$formula_mass) & 
            ids[row_ascr, ]$formula != "" &
            !is.na(sub$formula_mass)) {
          ids[row, ]$network_formula <- paste0(
            ifelse(ascr_mzdiff == -1, 
                   addElements(ids[row_ascr, ]$formula, sub$formula_mass),
                   subtractElements(ids[row_ascr, ]$formula, sub$formula_mass)))
        }else if (grepl("/",sub$formula_mass) & ## if multiple formula are in mass diff match
                  ids[row_ascr, ]$formula != "" &
                  !is.na(sub$formula_mass)) {
          md <- unlist(strsplit(# "CH4/NH2" ,#
            sub$formula_mass, 
            split = "/"))
          ids[row, ]$network_formula <- 
            ifelse(ascr_mzdiff == -1, 
                   paste0(sapply(md, function(x) addElements(x, ids[row_ascr, ]$formula)), collapse = "|"),
                   paste0(sapply(md, function(x) subtractElements(x, ids[row_ascr, ]$formula)), collapse = "|")
            )
          ### NOTE: If substraction makes no sense, result will be NA (e.g. XYN2 - N3)
          ###
          
        } else 
          next 
        
        
        
      }
  }
  # combine all results from all neighbors to results master file
  results <- rbind(results, ids[row & complete.cases(ids$network_annotation),])
  
}



# 7) unite results ------------------------------------------------------
# Remove duplicated results
results <- results[!duplicated(results), ]

# number of seed nodes used to annotate neighbor 
seedNo <- table(results$Row) %>% as.data.frame()
colnames(seedNo) <- c("Row", "seedNo")

# aggregate multiple results 
results_aggr <- aggregate(results,
                          list(results$Row),
                          function(x) {
                            paste0(unique(x), collapse ="|")
                          })

results_aggr <- merge(results_aggr, seedNo, by = "Row")

## keep only one entry per annotation:
## the one with the highest number of layers involved and/or 
## the smallest mass-difference
results_top <- results[order(c(results$network_layer, -results$mz_diff),  
                             decreasing = TRUE),] 

# remove duplicates
results_top <- results_top[ !duplicated(results_top$Row), ] 

# add number of seed nodes involved 
results_top <- merge(results_top, seedNo, by = "Row")
# 8) Export -------------------------------------------------------------------


write_csv(results_top, 
          file = paste0(path, "results/ascr_results_top.csv"))

