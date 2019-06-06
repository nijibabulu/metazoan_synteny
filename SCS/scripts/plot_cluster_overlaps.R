library(UpSetR)
library(purrr)
library(devtools)
load_all()

longnames <- load_data_set('longnames')
sp.names <- load_data_set('sp.names')

clusters <- load_data_set('clusters')
block_species <- getClusters(clusters) %>%  map(species)

int.df <- sp.names  %>% map(function(sp)
  map_dbl(block_species, ~if(sp %in% .) 1 else 0) %>% unname()) %>% 
  set_names(longnames[sp.names]) %>% 
  do.call(data.frame, .) %>% 
  set_names(longnames)

# give them nicer names

upset(int.df, sets=unlist(longnames[sp.names]), order.by='freq', point.size=4,
      text.scale=1.5, mainbar.y.label = 'Clusters', keep.order = T)
