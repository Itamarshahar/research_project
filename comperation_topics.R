comparisons_specific <- list(
  'X15.k1' = list(c('HA', 'SuperAgers'), c('HA', 'Young CTRL')),
  'X15.k2' = list(c('HA', 'SuperAgers'), c('SuperAgers', 'Young CTRL')),
  'X15.k3' = list(c('HD', 'MCI', 'SuperAgers'), c('SuperAgers', 'Young CTRL'))
  ,'X15.k14' = list(c('HD', 'HA','SuperAgers'), c('SuperAgers', 'Young CTRL'))
)

genarate_comperation_specific <- function(colname, diagmosis1, diagnosis2){
  rv <- list()
  rv[[colname]] <- list(diagmosis1, diagnosis2)
  return (rv)
}
