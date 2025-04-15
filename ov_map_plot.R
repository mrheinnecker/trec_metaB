# Install packages if not already installed
# install.packages(c("ggplot2", "rnaturalearth", "rnaturalearthdata", "sf"))

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(grid)  # for grobs
library(gridExtra)  # just in case
library(maps)

source("/g/schwab/Marco/repos/trec_metaB/fncts.R")


input_tables <- load_required_tables()

p_list=lapply(sort(input_tables$site_mapping$site %>% unique()), 
       get_stacked_barplot_per_site, 
       df_asvs=input_tables$df_asvs, 
       site_mapping=input_tables$site_mapping,
       mapping=input_tables$mapping)


# Example: filter for a specific city


all_sites <- c(sort(input_tables$site_mapping$site %>% unique()))
             #  'Roscoff', 'Tallinn', 'Kristineberg', 'Bilbao', 'Porto', 
             #   'Barcelona', 
            #   'Naples', 'Athens')  %>% unique()
data(world.cities)

cities <- subset(world.cities, name %in% all_sites) %>%
  select(name, lat, lon=long) %>%
  bind_rows(
    tibble(
      name = c("Roscoff", "Kristineberg",  "Villefranche"),
      lon = c(-3.9833, 11.4394,  7.3120),
      lat = c(48.7269, 58.2375, 43.7040)
    )
  ) %>%
  mutate(name=factor(name, levels=levels(input_tables$site_mapping$site))) %>%
  arrange(name)







# Load full world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create a data frame with city names and coordinates


# Plot Europe (approximate bounding box set in coord_sf only)
map <- ggplot(data = world) +
  geom_sf(fill = "antiquewhite") +
  geom_point(data = cities, aes(x = lon, y = lat), color = "red", size = 3) +
  geom_text(data = cities, aes(x = lon, y = lat, label = name), 
            nudge_y = 1, size = 3, color = "black") +
  coord_sf(
    xlim = c(-23, 32), ylim = c(35, 65)
           ) +  # Adjust view here
  theme_minimal() +
  labs(title = "Map of Europe with Selected Cities",
       x = NULL, y = NULL)



w <- 8
h <- 2

# Loop through cities and add each plot as a grob
for (i in seq_len(nrow(cities))) {
  #city <- cities$name[i]
  lon <- cities$lon[i]
  lat <- cities$lat[i]
  plot_grob <- ggplotGrob(p_list[[i]]+geom_col(show.legend = F)+
                            theme_void())
  
  map <- map + annotation_custom(
    grob = plot_grob,
    xmin = lon - w/2,
    xmax = lon + w/2,
    ymin = lat - h/2,
    ymax = lat + h/2
  )
}


legend_grob <- cowplot::get_legend(p_list[[1]]+geom_col()+theme(legend.title = element_blank()))


map+ annotation_custom(
  grob = legend_grob,
  xmin = -25, xmax = -15,  # Choose a good spot
  ymin = 48, ymax = 60
)

