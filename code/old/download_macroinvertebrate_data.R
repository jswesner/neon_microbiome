library(neonUtilities)
library(tidyverse)
library(lubridate)

# A) Download data ---------------------------------------------------------

# these are large files, and may take a while to fully download dependent on web traffic and connection speed. 

# download all macroinvertebrate colllection data from January 2017 to december 2019
macro <- loadByProduct(dpID = "DP1.20120.001",
                       startdate = "2017-01",
                       enddate = "2019-12",
                       check.size = FALSE,
                       nCores = 3)
#lakes to remove
lakes <- readRDS("~/GitHub/NEON_nsfproposal/size_spectra_proposal/data/lakes.rds")

sample_size <- macro$inv_fieldData %>% select(sampleID, benthicArea)

macro_df <- macro$inv_taxonomyProcessed %>% 
  left_join(sample_size) %>% 
  filter(!siteID %in% lakes) %>% 
  mutate(no_m2 = individualCount/benthicArea) %>% 
  # calculate total count / per m2 for each size class
  group_by(siteID, collectDate, genus, scientificName, family, order) %>%
  summarise(mean = mean(no_m2)) %>% 
  mutate(month = month(collectDate),
         year = year(collectDate))

site_lat_lon <- as_tibble(macro$inv_fieldData) %>% 
  distinct(siteID, decimalLatitude, decimalLongitude)

macro_df_morethan_15 <- macro_df  %>% 
  group_by(scientificName, year) %>% 
  tally() %>%
  filter(n >= 15) %>% 
  filter(!grepl("sp.", scientificName)) %>% 
  arrange(-n)




# sites with more than 15 and the number of taxa they contain
sites_with_morethan_15 <- macro_df %>% filter(scientificName %in% unique(macro_df_morethan_15$scientificName)) %>%  
  left_join(site_lat_lon) %>% 
  mutate(year = year(collectDate)) %>% 
  group_by(siteID, year) %>% 
  tally() %>% 
  arrange(-n) %>% 
  View()


world <- map_data("world")
states <- map_data("state")


site_label <- sites_with_morethan_10 %>% 
  select(decimalLatitude, decimalLongitude, siteID) %>% 
  distinct(siteID, .keep_all = T)

ggplot() + 
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey10", color = "white") + 
  geom_polygon(data = states, aes(x = long, y = lat, group =group), color = "white")+
  # coord_quickmap()+
  coord_fixed(1.2) +
  geom_label_repel(data = site_label, aes(x = decimalLongitude, y = decimalLatitude, label = siteID),
                   size = 3) +
  geom_point(data = sites_with_morethan_15, aes(x = decimalLongitude, y = decimalLatitude, label = siteID), size = 3,
             fill = "yellow", shape = 21) +
  ylim(c(10,75))+
  xlim(c(-170,-65)) +
  theme_void() +
  facet_wrap(~year)

