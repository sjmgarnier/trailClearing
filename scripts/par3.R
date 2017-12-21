# load required libraries
if (!require(pacman)) install.packages("pacman")
p_load(readr, dplyr, ggplot2, trailClearing)

numCores <- detectCores()
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)

# load data
dat <- read_csv("../data/data.csv") %>%
  mutate(Treatment = ifelse(Treatment == "LS", "15P",
                            ifelse(Treatment == "LH", "15C",
                                   ifelse(Treatment == "SS", "5P", "5C"))))

dat_15 <- filter(dat, Treatment == "15P" | Treatment == "15C")
dat_5 <- filter(dat, Treatment == "5P" | Treatment == "5C")

# calculate percolation distribution for full grids
full_grid_15 <- makeFullGrid(nrow = 20, ncol = 15)
distr_full_15 <- allPaths(full_grid_15)$l

full_grid_5 <- makeFullGrid(nrow = 20, ncol = 5)
distr_full_5 <- allPaths(full_grid_5)$l

# calculate grids for all samples
grids_15 <- group_by(dat_15, Treatment, Rep, Day) %>%
  do(grid = makeGrid(., nrow = 20, ncol = 15))

grids_5 <- group_by(dat_5, Treatment, Rep, Day) %>%
  do(grid = makeGrid(., nrow = 20, ncol = 5))

grids <- rbind_list(grids_15, grids_5)

# calculate random grids stats
res <- data.frame()

for (j in 1:nrow(grids)) {
  s <- ifelse(grids$Treatment[j] == "5C" | grids$Treatment[j] == "5P",
              sum(full_grid_5), sum(full_grid_15))

  n <- sum(grids$grid[[j]]) - s

  tmp <- data.frame()

  for (i in 1:1000) {
    print(paste0("Processing grid #", j, ": ", i, "/1000"))

    if (grids$Treatment[j] == "5C" | grids$Treatment[j] == "5P") {
      rnd_grid <- full_grid_5
      distr_full <- distr_full_5
    } else {
      rnd_grid <- full_grid_15
      distr_full <- distr_full_15
    }

    idx <- sample(which(rnd_grid == 0), n)
    rnd_grid[idx] <- 1

    perc <- allPaths(rnd_grid)$l

    tmp <- rbind_list(tmp, data.frame(median = median(perc), mean = mean(perc),
                                      min = min(perc), max = max(perc),
                                      overlap = overlap(distr_full, perc)))
  }

  tmp$idx <- j

  res <- rbind_list(res, tmp)
}

stopCluster(cl)




