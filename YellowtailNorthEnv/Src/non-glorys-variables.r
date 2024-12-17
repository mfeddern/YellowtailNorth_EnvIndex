source("_00_yellowtail-header.r")


pdo = data.frame(read.csv( paste0(data_dir,'cciea_OC_PDO.csv')))


process_envir <- function(dfile=pdo, index_name = 'PDO'){
  dfile = dfile %>% rename(date = time) %>% 
    mutate(date = substring(date,1,10))
  dfile = dfile[dfile$date !="UTC",]
  dfile$year = year(dfile$date)
  dfile$month = month(dfile$date)
  dfile$index = as.numeric(dfile[,index_name])
  dfile$season = 'Spr'
  dfile$season[dfile$month %in% c(7:9)] <- "Sum"
  dfile$season[dfile$month %in% c(10:12)] <- "Fal"
  dfile$season[dfile$month %in% c(1:3)] <- "Win"
  dfile = dfile %>% group_by(year, season) %>% 
    reframe(index = mean(index)) %>%
    mutate(index_name = paste0(index_name,"_",season))
  dfile = dfile %>% dplyr::select(!season)
  dfilex = pivot_wider(dfile, values_from = index, names_from = index_name) 
  return(dfilex) # file to spit out
}

PDO = process_envir(pdo,"PDO")
head(PDO)




