library(tidyverse)
library(readxl)
library(brms)
library(janitor)
library(tidybayes)
library(ggridges)
library(snakecase)
library(ggthemes)
library(cowplot)


sd_fun = function(model){
  as_draws_df(model) %>% 
    clean_names() %>% 
    select(starts_with(c("sd_"))) %>%
    mutate(.draw = row_number()) %>% 
    pivot_longer(cols = -.draw) %>% 
    separate(name, into = c(NA, "factor", NA)) %>% 
    mutate(factor = case_when(factor == "order" ~ "Insect Order",
                              factor == "locations" ~ "Stream",
                              factor == "functionalgroup" ~ "Functional Group",
                              factor == "genus" ~ "Genus",
                              factor == "family" ~ "Family")) 
} # posterior draws of the varying intercepts

# get data
sample_ids <- read_csv("data/sample_ids.csv") # might need for later

alpha_beta_diversity = read_excel("data/GLM MODELLING DATA.xlsx", 
                                  sheet = "mapping file") %>% 
  clean_names(case = "big_camel") %>% 
  mutate(BrayCurtis = case_when(BrayCurtis == 0 ~ 0.01,
                                BrayCurtis == 1 ~ 0.99,
                                 TRUE ~ BrayCurtis),
         FunctionalGroup = str_replace_all(FunctionalGroup, "[^[:alnum:]]", "")) %>% 
  rename(Functionalgroup = FunctionalGroup) %>% 
  select(-X12) %>% 
  mutate(GenusTaxized = case_when(GenusTaxized == "NA" ~ "unknown",
                                  TRUE ~ GenusTaxized))



# beta diversity ----------------------------------------------------------

# prior predictive
# sample priors
prior_bray <- brm(BrayCurtis ~ 1 + (1|Functionalgroup) + (1|Family) + (1|GenusTaxized) + (1|Locations),
                  family = Beta(link = "logit"),
                  data = alpha_beta_diversity,
                  prior = c(prior(normal(0, 1), class = "Intercept"),
                             prior(exponential(0.5), class = "sd")),
                  sample_prior = "only",
                  iter = 100, chains = 1,
                  file = "models/prior_bray.rds",
                  file_refit = "on_change")


saveRDS(prior_bray, file = "models/prior_bray.rds")


# fit model
fit_bray = update(prior_bray, sample_prior = "no",
                  newdata = alpha_beta_diversity,
                  iter = 2000, chains = 4,
                  file = "models/fit_bray.rds",
                  file_refit = "on_change")

pp_check(fit_bray, type = "stat_grouped", group = "Family")
pp_checkplot = pp_check(fit_bray, type = "boxplot")

ggsave(pp_checkplot, file = "plots/pp_checkplot.jpg", dpi = 600, width = 5, height = 2)


# Plot --------------------------------------------------------------------

# data to plot
order_groups = alpha_beta_diversity %>% distinct(Order, GenusTaxized, Family) %>% rename(Genus = GenusTaxized) %>% 
  mutate(order_color = Order) %>% 
  pivot_longer(cols = c(-order_color)) %>% 
  select(-name)

posts_to_plot <- fit_bray$data %>% 
  add_epred_draws(fit_bray, re_formula = NULL) %>% 
  rename(Genus = GenusTaxized,
         Stream = Locations) %>% 
  pivot_longer(cols = c("Functionalgroup", "Family", "Genus", "Stream")) %>% 
  left_join(order_groups) %>%  
  group_by(value, order_color) %>% 
  mutate(median = median(.epred)) %>% 
  mutate(taxonomic = case_when(name %in% c("Functionalgroup") ~ "Trophic",
                               name == "Stream" ~ "Geography",
                               TRUE ~"Taxonomy")) %>% 
  mutate(name = fct_relevel(name, "Family", "Genus", "Functionalgroup")) 

post_medians = posts_to_plot %>% distinct(name, value, median)

data_to_plot = fit_bray$data %>%
  rename(Genus = GenusTaxized,
         Stream = Locations) %>% 
  pivot_longer(cols = -BrayCurtis) %>% 
  left_join(post_medians) %>% 
  mutate(taxonomic = case_when(name %in% c("Functionalgroup") ~ "Trophic",
                               name == "Stream" ~ "Geography",
                               TRUE ~"Taxonomy")) %>% 
  mutate(name = fct_relevel(name, "Family", "Genus", "Functionalgroup")) %>% 
  left_join(order_groups)


sd_datatoplot = sd_fun(fit_bray)

# make plots
library(viridis)

taxa_grid = posts_to_plot %>%
  group_by(order_color) %>% 
  mutate(order_median = median(.epred)) %>% 
  filter(taxonomic == "Taxonomy") %>% 
  filter(value != "unknown") %>%
  mutate(name = fct_relevel(as_factor(name), "Family")) %>% 
  ggplot(aes(x = .epred, y = reorder(value, order_median), group = value, fill = reorder(order_color, median))) + 
  stat_slab(alpha = 0.8) +
  theme_default() +
  facet_wrap(~name, scales = "free_y") +
  scale_fill_viridis_d(name = "Order") +
  guides(alpha = "none") + 
  labs(y = "",
       x = "Beta Diversity (Bray Curtis Dissimilarity)") +
  geom_point(data = data_to_plot %>%
               group_by(order_color) %>% 
               mutate(order_median = median(BrayCurtis)) %>%
               filter(taxonomic == "Taxonomy") %>% distinct() %>% filter(value != "unknown"), aes(x = BrayCurtis),
             shape = "|")

nontaxa_grid = posts_to_plot %>%
  filter(taxonomic != "Taxonomy") %>%
  mutate(name = str_replace(name, "Functionalgroup", "Functional Group")) %>% 
  ggplot(aes(x = .epred, y = reorder(value, median))) + 
  stat_halfeye(color = "black") +
  theme_default() +
  facet_wrap(~name, scales = "free_y") +
  scale_fill_viridis_d() +
  guides(fill = "none") + 
  labs(y = "",
       x = "Beta Diversity (Bray Curtis Dissimilarity)") +
  geom_point(data = data_to_plot%>%
               mutate(name = str_replace(name, "Functionalgroup", "Functional Group")) %>%
               filter(taxonomic 
                      != "Taxonomy") %>% distinct(), aes(x = BrayCurtis),
             shape = "|")

bray_sd_plot = sd_datatoplot %>% 
  group_by(factor) %>% 
  mutate(median = median(value)) %>% 
  ggplot(aes(x = value, y = reorder(factor, median)), color = "black") + 
  stat_halfeye() +
  guides(fill = "none") + 
  theme_default() +
  labs(subtitle = "Variation attributable to each grouping",
       x = "Standard Deviation",
       y = "Grouping") + 
  coord_cartesian(xlim = c(NA, 3))



library(patchwork)

taxa_nontaxa_grid = taxa_grid/nontaxa_grid
ggsave(taxa_nontaxa_grid, file = "plots/taxa_nontaxa_grid.jpg", dpi = 600, width = 8, height = 8)
ggsave(bray_sd_plot, file = "plots/bray_sd_plot.jpg", dpi = 600, width = 5, height = 4)




# Prior versus Posterior --------------------------------------------------

# plot draws by group
priors_to_plot = prior_bray$data %>% 
  add_epred_draws(prior_bray, re_formula = NULL) %>% 
  rename(Genus = GenusTaxized,
         Stream = Locations) %>% 
  pivot_longer(cols = c("Functionalgroup", "Family", "Genus", "Stream")) %>% 
  left_join(order_groups) %>%  
  group_by(value, order_color) %>% 
  mutate(median = median(.epred)) %>% 
  mutate(taxonomic = case_when(name %in% c("Functionalgroup") ~ "Trophic",
                               name == "Stream" ~ "Geography",
                               TRUE ~"Taxonomy")) %>% 
  mutate(name = fct_relevel(name, "Family", "Genus", "Functionalgroup")) 

post_medians = posts_to_plot %>% distinct(name, value, median)

prior_v_post = egg::tag_facet(posts_to_plot %>% mutate(model = "posterior") %>% 
  bind_rows(priors_to_plot %>% mutate(model = "prior")) %>%
  # filter(value == "Tipulidae") %>%
  # filter(taxonomic == "Taxonomy") %>% 
  filter(value != "unknown") %>%
  mutate(name = fct_relevel(as_factor(name), "Family")) %>% 
  ggplot(aes(x = .epred, 
             group = interaction(model, value), 
             color = model)) + 
  # stat_slab() +
  geom_density() +
  theme_default() +
  facet_wrap(~name, scales = "free_y") +
  scale_color_colorblind() +
  guides(alpha = "none") + 
  labs(y = "Prior or Posterior Density",
       color = "Model",
       x = "Beta Diversity (Bray Curtis Dissimilarity)")
)
  

ggsave(prior_v_post, file = "plots/prior_v_post.jpg", dpi = 600, width = 6, height = 4)

prior_bray_draws %>% 
  ggplot(aes(x = level, y = value)) + 
  geom_violin() + 
  facet_wrap(~group, nrow = 1, scales = "free_x")

# Summarize ---------------------------------------------------------------

bray_means_95CI = posts_to_plot %>% 
  group_by(name, value) %>% 
  mean_qi(.epred) %>% 
  rename(mean = .epred)

sd_means_95CI = sd_datatoplot %>% 
  group_by(factor) %>% 
  mean_qi(value) %>% 
  arrange(-value) %>% 
  rename(mean = value)


write_csv(bray_means_95CI, file = "tables/bray_means_95CI.csv")
write_csv(sd_means_95CI, file = "tables/sd_means_95CI.csv")


as_draws_df(fit_bray) %>% as_tibble() %>% clean_names() %>% 
  mutate(b_intercept = inv_logit_scaled(b_intercept)) %>% 
  summarize(mean = mean(b_intercept),
            sd = sd(b_intercept))


sd_datatoplot %>% 
   pivot_wider(values_from = value, names_from = factor) %>% 
   clean_names() %>% 
   mutate(ffg_family = functional_group - family,
          ffg_genus = functional_group - genus,
          ffg_stream = functional_group - stream) %>% 
   pivot_longer(cols = contains("ffg")) %>% 
   group_by(name) %>% 
   summarize(prob_ffg_greater = sum(value > 0)/max(draw))

sd_datatoplot %>% 
  pivot_wider(values_from = value, names_from = factor) %>% 
  clean_names() %>% 
  mutate(ffg_family = functional_group - family,
         ffg_genus = functional_group - genus,
         ffg_stream = functional_group - stream) %>% 
  pivot_longer(cols = contains("ffg")) %>% 
  ggplot(aes(x = value, y = name)) + 
  stat_halfeye()


fit_bray$data %>% as_tibble() %>% clean_names() %>% 
  ungroup() %>% 
  # count(family) %>% 
  filter(grepl("Perli", family))
