# devtools::install_github("daijiang/neonDivData")
library(neonDivData)
library(tidyverse)
library(janitor)
library(taxize)
library(lubridate)


# Algae -------------------------------------------------------------------


neon_site_info <- read.csv("neon_site_info.csv") %>% clean_names() %>% separate(domain_site, c("domain", "siteID", by = " - ")) %>%
  select(-" - ")
neon_data <- data_algae
dunbeck_data <- read.csv("dunbeck_density.csv") %>% rename(chl_a = total_chl)
chl_a_neon <- read.csv("chla_temp_table.csv") %>% clean_names() %>% rename(siteID = site)


combined_data <- data_algae %>%
  group_by(unit, siteID, observation_datetime, unique_sample_id) %>%
  summarize(value = sum(value)) %>%
  left_join(chl_a_neon) %>%
  left_join(neon_site_info) %>%
  filter(unit == "cells/mL") %>%
  mutate(source = "neon") %>%
  select(siteID, value, chl_a, site_type) %>%
  mutate(source = "neon")  %>%
  bind_rows(dunbeck_data %>% select(site, abundance, chl_a) %>% rename(siteID = site, value = abundance) %>%
              mutate(siteID = as.character(siteID),
                     source = "dunbeck",
                     site_type = "River"))


combined_data %>%
ggplot(aes(x = value, y = ..scaled.., color = interaction(site_type, source))) +
  geom_density() +
  scale_x_log10()




data_algae %>%
  left_join(neon_site_info) %>%
  group_by(unit, siteID, observation_datetime, site_type) %>%
  mutate(value = sum(value)) %>%
  filter(unit == "cells/mL") %>%
  mutate(date = ymd(as.Date(observation_datetime)),
         month = month(date),
         year = year(date)) %>%
  ggplot(aes(x = month, y = value, color = siteID)) +
  geom_point() +
  # geom_boxplot(aes(group = month)) +
  geom_line(aes(group = interaction(siteID,year))) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(2, 6, 10),
                     labels = c("2", "6", "10")) +
  labs(y = "cells/ml") +
  facet_grid(~site_type) +
  # coord_cartesian(ylim = c(0, 100000)) +
  NULL



# Macroinvertebrates ------------------------------------------------------

neon_inverts <- data_macroinvertebrate %>%
  separate(taxon_name, c('genus', 'species'), remove = F)

# invert_taxa <- neon_inverts %>% distinct(taxon_name) %>%
#   distinct(genus)

# invert_taxize <- classification(invert_taxa$genus, db = "ncbi")
# saveRDS(invert_taxize, file = "data/invert_taxize.rds")

invert_taxize <- readRDS(file = "data/invert_taxize.rds")

taxized_inverts <- rbind(invert_taxize) %>% as_tibble() %>%
  # select(-id) %>%
  distinct(name, rank, query) %>%
  filter(rank %in% c("order", "family", "class", "genus")) %>%
  pivot_wider(names_from = rank, values_from = name) %>%
  rename(genus_taxized = genus,
         genus = query)

neon_inverts_taxized <- left_join(neon_inverts, taxized_inverts)




# filter for NEON samples
sites <- c("LEWI", "MART", "BLUE", "BLDE", "LECO", "BIGC", "POSE", "HOPB", "MCRA", "MCDI")
dates <- c("2019-03-19", "2019-09-30", "2019-10-09", "2019-09-04", "2019-10-14", "2019-10-01", "2019-11-04",
           "2019-10-03", "2019-09-25", "2020-06-05")
sample_type <- c("core",
                 "surber",
                 "kicknet",
                 rep("surber", 7))

filter_data <- tibble(sites, dates, sample_type) %>%
  unite("merge_id", c(sites, dates, sample_type))




filtered_data <- neon_inverts_taxized %>%
  mutate(date = as.Date(observation_datetime)) %>%
  unite("merge_id", c(siteID, date, samplerType), remove = F) %>%
  right_join(filter_data)


filtered_data %>%
  filter(class == "Insecta") %>%
  group_by(siteID, merge_id, samplerType, unit, taxon_name, taxon_id, order, family) %>%
  tally() %>%
  group_by(siteID) %>%
  mutate(order_n = max(n)) %>%
  ggplot(aes(x = n, y = reorder(order, n), color = order)) +
  geom_point(position = position_jitter(width = 0.05),
             shape = 21) +
  facet_wrap(~siteID)




taxa_by_site <- filtered_data %>%
  # filter(siteID == "MCRA") %>%
  group_by(order, family, taxon_name, siteID, genus_taxized) %>%
  tally() %>%
  ungroup() %>%
  arrange(siteID, order, -n)

write.csv(taxa_by_site , file = "data/taxa_by_site.csv", row.names = F)





# Wolbachia NEON loaned samples

site_sample_info <- read_csv("data/site_sample_info.csv") %>%
  filter(!is.na(sample_id))

taxa_for_wolbachia <- read_csv("data/taxa_we_have.csv") %>%
  left_join(site_sample_info) %>%
  filter(!is.na(n_wolbachia))
