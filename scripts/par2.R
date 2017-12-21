# load required libraries
if (!require(pacman)) install.packages("pacman")
p_load(readr, dplyr, foreach, doSNOW, tcltk, parallel, trailClearing)

# load data
dat <- read_csv("../data/data.csv") %>%
  mutate(Treatment = ifelse(Treatment == "LS", "15P",
                            ifelse(Treatment == "LH", "15C",
                                   ifelse(Treatment == "SS", "5P", "5C"))))

dat_15 <- filter(dat, Treatment == "15P" | Treatment == "15C")
dat_5 <- filter(dat, Treatment == "5P" | Treatment == "5C")

# calculate percolation distribution for full grids
full_grid_15 <- makeFullGrid(nrow = 20, ncol = 15)
distr_full_15 <- percolate2(full_grid_15, n = 100000)

full_grid_5 <- makeFullGrid(nrow = 20, ncol = 5)
distr_full_5 <- percolate2(full_grid_5, n = 100000)

# calculate grids for all samples
grids_15 <- group_by(dat_15, Treatment, Rep, Day) %>%
  do(grid = makeGrid(., nrow = 20, ncol = 15))

grids_5 <- group_by(dat_5, Treatment, Rep, Day) %>%
  do(grid = makeGrid(., nrow = 20, ncol = 5))

grids <- rbind_list(grids_15, grids_5)

# calculate random grids stats
numCores <- detectCores()
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)

res <- data.frame()

pb <- txtProgressBar(max = 1000, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)

tmp <- foreach(i = 1:nrow(grids), .combine = rbind_list, .packages = "trailClearing", .options.snow = list(progress = progress)) %dopar% {

  s <- ifelse(grids$Treatment[i] == "5C" | grids$Treatment[i] == "5P",
              sum(full_grid_5), sum(full_grid_15))
  n <- sum(grids$grid[[i]]) - s

  if (grids$Treatment[i] == "5C" | grids$Treatment[i] == "5P") {
    rnd_grid <- full_grid_5
    distr_full <- distr_full_5
  } else {
    rnd_grid <- full_grid_15
    distr_full <- distr_full_15
  }

  idx <- sample(which(rnd_grid == 0), n)
  rnd_grid[idx] <- 1

  perc <- percolate2(rnd_grid, n = 100000)

  data.frame(median = median(perc), min = min(perc), max = max(perc),
             overlap = overlap(distr_full, perc))
}






for (j in 1:nrow(grids)) {
  s <- ifelse(grids$Treatment[j] == "5C" | grids$Treatment[j] == "5P",
              sum(full_grid_5), sum(full_grid_15))

  n <- sum(grids$grid[[j]]) - s


  tmp <- foreach (i = 1:1000, .combine = rbind_list, .packages = "trailClearing",
                  .options.snow = list(progress = progress)) %dopar% {
                    if (grids$Treatment[j] == "5C" | grids$Treatment[j] == "5P") {
                      rnd_grid <- full_grid_5
                      distr_full <- distr_full_5
                    } else {
                      rnd_grid <- full_grid_15
                      distr_full <- distr_full_15
                    }

                    idx <- sample(which(rnd_grid == 0), n)
                    rnd_grid[idx] <- 1

                    perc <- percolate2(rnd_grid, n = 100000)

                    data.frame(median = median(perc), min = min(perc), max = max(perc),
                               overlap = overlap(distr_full, perc))
                  }

  tmp$idx <- j

  res <- rbind_list(res, tmp)

  close(pb)
}

stopCluster(cl)