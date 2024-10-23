# key <- vroom::vroom('~/keyDNAextractions...') #CHANGE
ASVtbl <- vroom::vroom('~/Downloads/ITS_data6.fungi.seq.tax.rar.csv') #CHANGE
library(dplyr)
glimpse(ASVtbl); str(ASVtbl) #suggestions: ASVtbl$  #anyOK
# format species table as typical (samples x species), needed by vegan::diversity()
ASVtbl %>% select(asv_id, starts_with('X')) %>% #dim() #shows268samples/268col's
  
  # format as table.frame in base R while keeping rownames/IDs
  tibble::column_to_rownames('asv_id') %>% 
  t() %>% #.[1:3,1:3] #see
  as_tibble(rownames = NA) %>% #keeps rowname info, but invisible w/ '*' symbol
  
  # left_join(key) %>% #distingush samples
  # group_by(PLOT, SESSION, SUBPLOT, ...)
  
  #calculate while keeping rest of table
  mutate('diversity' = vegan::diversity(.)) %>% 
  
  # vegan::diversity() cant have ID columns, so moved to after calculation
  tibble::rownames_to_column('ASVs') %>% 

  ggplot(aes(diversity)) + geom_histogram()
  
  # ggplot(aes(SESSION, diversity)) + 
  # geom_point(color = 'gray') + 
  # stat_summary() +
  # facet_wrap(~ PLOT)