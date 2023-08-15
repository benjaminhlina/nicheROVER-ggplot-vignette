# ---- Bring in Packages ---- 
{
  library(dplyr)
  library(ellipse)
  library(ggplot2)
  library(ggtext)
  library(here)
  library(nicheROVER) 
  library(purrr)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tidyr)
}

# ---- Bring in example data; will need to replace with your own ----

df <- fish %>% 
  janitor::clean_names() %>% 
  select(-d34s)


# set the number of posterior samples that will be taken 
nsample <- 1000


# ---- Use map from purrr and niw.post from nicheROVER to estimate niches ----- 
fish_par <- df %>% 
  split(.$species) %>% 
  map(~ select(., d15n, d13c)) %>% 
  map(~niw.post(nsample = nsample, X = .))


# ---- Extract mu from nicheROVER object ----- 
df_mu <- map(fish_par, pluck, 1) %>% 
  imap(~ as_tibble(.x) %>% 
         mutate( 
           metric = "mu", 
           species = .y
         )
  ) %>%
  bind_rows() %>% 
  mutate(
    species = factor(species, 
                     levels = c("ARCS", "BDWF", "LKWF", "LSCS"))
  ) %>% 
  group_by(species) %>% 
  mutate(
    sample_number = 1:1000
  ) %>% 
  ungroup()




df_mu_long <- df_mu %>% 
  pivot_longer(cols = -c(metric, species, sample_number), 
               names_to = "isotope", 
               values_to = "mu_est") %>% 
  mutate(
    element = case_when(
      isotope == "d15n" ~ "N",
      isotope == "d13c" ~ "C",
    ), 
    neutron = case_when(
      isotope == "d15n" ~ 15,
      isotope == "d13c" ~ 13,
    ) 
  )

# ---- Extract sigma from nicheROVER object ----- 

df_sigma <- map(fish_par, pluck, 2) %>%
  imap(~ as_tibble(.x) %>%
         mutate(
           metric = "sigma",
           id = c("d15n", "d13c"),
           species = .y
         )
  ) %>%
  bind_rows() %>%
  pivot_longer(cols = -c("id", "species", "metric"),
               names_to = "isotope",
               values_to = "post_sample"
  )  %>%
  separate(isotope, into = c("isotopes", "sample_number"), sep = "\\.")



df_sigma_cn <- df_sigma %>%
  filter(id != isotopes)

# ---- Plot posterior samples for mu estimates from nicheROver ----- 
posterior_plots <- df_mu_long %>%
  split(.$isotope) %>%
  imap(
    ~ ggplot(data = ., aes(x = mu_est)) +
      geom_density(aes(fill = species), alpha = 0.5) +
      scale_fill_viridis_d(begin = 0.25, end = 0.75,
                           option = "D", name = "Species") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title.x =  element_markdown(),
            axis.title.y =  element_markdown(),
            legend.position = "none"
      ) +
      labs(
        x = paste("\u00b5<sub>\U03B4</sub>", "<sub><sup>",
                  unique(.$neutron), "</sup></sub>",
                  "<sub>",unique(.$element), "</sub>", sep = ""),
        y = paste0("p(\u00b5 <sub>\U03B4</sub>","<sub><sup>",
                   unique(.$neutron), "</sub></sup>",
                   "<sub>",unique(.$element),"</sub>",
                   " | X)"), sep = "")
  )

posterior_plots$d15n +
  theme(legend.position = c(0.18, 0.84)) + 
  posterior_plots$d13c

# ---- Prepare sigma dataframes for plotting ----

df_sigma_cn <- df_sigma_cn %>%
  mutate(
    element_id = case_when(
      id == "d15n" ~ "N",
      id == "d13c" ~ "C",
    ),
    neutron_id = case_when(
      id == "d15n" ~ 15,
      id == "d13c" ~ 13,
    ),
    element_iso = case_when(
      isotopes == "d15n" ~ "N",
      isotopes == "d13c" ~ "C",
    ),
    neutron_iso = case_when(
      isotopes == "d15n" ~ 15,
      isotopes == "d13c" ~ 13,
    )
  )

# ---- Plot posterior samples for mu estimates from nicheROver ----- 
sigma_plots <- df_sigma_cn %>%
  group_split(id, isotopes) %>%
  imap(
    ~ ggplot(data = ., aes(x = post_sample)) +
      geom_density(aes(fill = species), alpha = 0.5) +
      scale_fill_viridis_d(begin = 0.25, end = 0.75,
                           option = "D", name = "Species") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.title.x =  element_markdown(),
            axis.title.y =  element_markdown(),
            legend.position = "none"
      ) +
      labs(
        x = paste("\U03A3","<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_id), "</sub></sup>",
                  "<sub>",unique(.$element_id),"</sub>"," ",
                  "<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_iso), "</sub></sup>",
                  "<sub>",unique(.$element_iso),"</sub>", sep = ""),
        y = paste("p(", "\U03A3","<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_id), "</sub></sup>",
                  "<sub>",unique(.$element_id),"</sub>"," ",
                  "<sub>\U03B4</sub>",
                  "<sub><sup>", unique(.$neutron_iso), "</sub></sup>",
                  "<sub>",unique(.$element_iso),"</sub>", " | X)", sep = ""),
      )
  )

sigma_plots[[1]] + 
  theme(legend.position = c(0.1, 0.82))


# ---- Manipulate sigma dataframes for ellipse loops ----- 
df_sigma_wide <- df_sigma %>%
  pivot_wider(names_from = id,
              values_from = post_sample)


# ---- Set the confidence interval for ellipses ----
p.ell <- 0.95

# ---- create blank list to dump results and category object to loop for 
species_name <- unique(df_sigma_wide$species)



all_ellipses <- list()


# ---- Loop to create ellipse ------ 
for (i in 1:length(species_name)) {
  
  sigma_species <- df_sigma_wide %>% 
    filter(species %in% species_name[i])
  
  mu_species <- df_mu %>% 
    filter(species %in% species_name[i])
  
  ell <- NULL
  post.id <- NULL
  
  for(j in 1:length(unique(sigma_species$sample_number))) {
    
    sigma_ind <- sigma_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(d15n, d13c) 
    
    Sigma <- as.matrix(sigma_ind, 2, 2)
    row.names(Sigma) <- c("d15n", "d13c")
    
    mu <- mu_species %>%
      filter(sample_number %in% sample_number[j]) %>% 
      dplyr::select(sample_number, d15n, d13c) %>% 
      pivot_longer(cols = -sample_number, 
                   names_to = "isotope", 
                   values_to = "mu") %>% 
      .$mu
    
    out <- ellipse::ellipse(Sigma, centre = mu, which = c(1, 2), level = p.ell)
    
    ell <- rbind(ell, out)
    post.id <- c(post.id, rep(j, nrow(out)))
  }
  ell <- as.data.frame(ell)
  ell$rep <- post.id
  all_ellipses[[i]] <- ell
}

# ---- Combine ellipse list into dataframe and add species names back in -----
ellipse_df <- bind_rows(all_ellipses, .id = "id") %>% 
  mutate(
    species = factor(
      case_when(
        id == "1" ~ "ARCS",
        id == "2" ~ "BDWF",
        id == "3" ~ "LKWF",
        id == "4" ~ "LSCS",
      ), level = c("ARCS", "BDWF", "LKWF", "LSCS")
    )
  ) %>% 
  as_tibble()

# ---- Randomly sample 10 ellipse for plotting ----- 

ellipse_df %>% 
  group_by(species, rep) %>% 
  nest() %>%
  group_by(species) %>% 
  slice_sample(n = 10, replace = TRUE) %>% 
  ungroup() %>% 
  unnest(cols = c(data)) -> random_ellipse 


# ---- Plot ellipse, biplot and each isotopes density ----- 
ellipse_plots <- ggplot() + 
  geom_polygon(data = random_ellipse,
               mapping = aes(x = d13c, y = d15n,
                             group = interaction(rep, species),
                             color = species),
               fill = NA,
               linewidth = 0.5) + 
  
  scale_colour_viridis_d(begin = 0.25, end = 0.75, 
                         option = "D", name = "species",
  ) + 
  scale_x_continuous(breaks = rev(seq(-20, -40, -2))) +
  scale_y_continuous(breaks = seq(6, 16, 2)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N")))



iso_long <- df %>%
  pivot_longer(cols = -species,
               names_to = "isotope", 
               values_to = "value") %>% 
  mutate(
    element = case_when(
      isotope == "d15n" ~ "N",
      isotope == "d13c" ~ "C",
    ), 
    neutron = case_when(
      isotope == "d15n" ~ 15,
      isotope == "d13c" ~ 13,
    )
  )



iso_density <- iso_long %>% 
  group_split(isotope) %>% 
  imap(
    ~ ggplot(data = .) + 
      geom_density(aes(x = value, 
                       fill = species), 
                   alpha = 0.35, 
                   linewidth = 0.8) +
      scale_fill_viridis_d(begin = 0.25, end = 0.75,
                           option = "D", name = "Species") +
      theme_bw(base_size = 10) +
      theme(axis.text = element_text(colour = "black"),
            panel.grid = element_blank(), 
            legend.position = c(0.15, 0.65), 
            legend.title.align = 0.5,
            legend.background = element_blank(), 
            axis.title.x = element_markdown(family = "sans")) + 
      labs(x =  paste("\U03B4",
                      "<sup>", unique(.$neutron), "</sup>",unique(.$element), 
                      sep = ""), 
           y = "Density")
  )

d13c_density <- iso_density[[1]] + 
  scale_x_continuous(breaks = rev(seq(-20, -34, -2)),
                     limits = rev(c(-20, -34)))

d15n_density <- iso_density[[2]] +
  scale_x_continuous(breaks = seq(5, 15, 2.5), 
                     limits = c(5, 15)) + 
  theme(
    legend.position = "none"
  )



iso_biplot <- ggplot() + 
  geom_point(data = df, aes(x = d13c, y = d15n,
                            fill = species),
             shape = 21, colour = "black", 
             stroke = 0.8,
             size = 3, alpha = 0.70) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "species") +
  scale_x_continuous(breaks = rev(seq(-20, -39, -1))) +
  scale_y_continuous(breaks = seq(5, 17, 1)) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(), 
        legend.position = "none", 
        legend.title.align = 0.5,
        legend.background = element_blank()) + 
  labs(x = expression(paste(delta ^ 13, "C")), 
       y = expression(paste(delta ^ 15, "N")))


d13c_density + ellipse_plots + iso_biplot + d15n_density +
  plot_annotation(tag_levels = "a", 
                  tag_suffix = ")")


# ---- Estimate niche similarities with nicheROVER ----- 
over_stat <- overlap(fish_par, nreps = nsample, nprob = 1000, 
                     alpha = 0.95)


# ---- Extract niche similarities and convert into a dataframe ---- 
over_stat_df <- over_stat %>% 
  as_tibble(rownames = "species_a") %>% 
  mutate(
    id = 1:nrow(.), 
    species_a = factor(species_a, 
                       level = c("ARCS", "BDWF", "LKWF", "LSCS"))
  ) %>% 
  pivot_longer(cols = -c(id, species_a), 
               names_to = "species_b", 
               values_to = "mc_nr")  %>% 
  separate(species_b, into = c("species_c", "sample_number"), 
           sep = "\\.") %>% 
  select(-id) %>% 
  rename(species_b = species_c) %>% 
  mutate(
    species_b =  factor(species_b, 
                        level = c("ARCS", "BDWF", "LKWF", "LSCS")
    ), 
    mc_nr_perc = mc_nr * 100
  )



over_sum <- over_stat_df %>% 
  group_by(species_a, species_b) %>% 
  summarise(
    mean_mc_nr = round(mean(mc_nr_perc), digits = 2),
    qual_2.5 = round(quantile(mc_nr_perc, probs = 0.025, na.rm = TRUE), digits = 2), 
    qual_97.5 = round(quantile(mc_nr_perc, probs = 0.975, na.rm = TRUE), digits = 2)
  ) %>% 
  ungroup() %>% 
  pivot_longer(cols = -c(species_a, species_b, mean_mc_nr), 
               names_to = "percentage", 
               values_to = "mc_nr_qual") %>% 
  mutate(
    percentage = as.numeric(str_remove(percentage, "qual_"))
  ) 


# ---- plot niche similarities ----- 
ggplot(data = over_stat_df, aes(x = mc_nr_perc)) + 
  geom_density(aes(fill = species_a)) + 
  geom_vline(data = over_sum, aes(xintercept = mean_mc_nr), 
             colour = "black", linewidth = 1) +
  geom_vline(data = over_sum, aes(xintercept = mc_nr_qual), 
             colour = "black", linewidth = 1, linetype = 6) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75,
                       option = "D", name = "Species", 
                       alpha = 0.35) + 
  ggh4x::facet_grid2(species_a ~ species_b, 
                     independent = "y",
                     scales = "free_y") + 
  theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    axis.text = element_text(colour = "black"), 
    legend.background = element_blank(),
    strip.background = element_blank()
  ) +
  labs(x = paste("Overlap Probability (%)", "\u2013", 
                 "Niche Region Size: 95%"), 
       y = "p(Percent Overlap | X)")


# ---- Estimate niche size ----- 
niche_size <- sapply(fish_par, function(spec) {
  apply(spec$Sigma, 3, niche.size)
})


# ---- Convert niche size into a dataframe ----- 
niche_size_df <- niche_size %>% 
  as_tibble() %>% 
  mutate(
    id = 1:nrow(.)
  ) %>% 
  pivot_longer(
    cols = -id, 
    names_to = "species", 
    values_to = "niche_size"
  ) %>% 
  mutate(
    id = 1:nrow(.), 
    species = factor(species,
                     level = c("ARCS", "BDWF", 
                               "LKWF", "LSCS"))
  )


# ---- Extract niche size mean, sd, and sem ----- 
niche_size_mean <- niche_size_df %>% 
  group_by(species) %>% 
  summarise(
    mean_niche = round(mean(niche_size), digits = 2), 
    sd_niche = round(sd(niche_size), digits = 2), 
    sem_niche = round(sd(niche_size) / sqrt(n()), digits = 2)
  ) %>% 
  ungroup()


# ---- Plot niche sizes as violin plots ------ 
ggplot(data = niche_size_df) + 
  geom_violin(
    aes(x = species, y = niche_size),
    width = 0.2) + 
  geom_point(data = niche_size_mean, aes(x = species, y = mean_niche)) +
  geom_errorbar(data = niche_size_mean, aes(x = species, 
                                            ymin = mean_niche  - sem_niche, 
                                            ymax = mean_niche  + sem_niche), 
                width = 0.05) + 
  theme_bw(base_size = 15) + 
  theme(panel.grid = element_blank(), 
        axis.text = element_text(colour = "black")) + 
  labs(x = "Species", 
       y = "Niche Size") 

