# load required libraries
if (!require(pacman)) install.packages("pacman")
p_load(readr, dplyr, foreach, doSNOW, tcltk, parallel, trailClearing)

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
distr_full_15 <- diag(bfSearch(full_grid_15, start = cbind(x = 1, y = 1:nrow(full_grid_15)),
                          end = cbind(x = ncol(full_grid_15), y = 1:nrow(full_grid_15))))

full_grid_5 <- makeFullGrid(nrow = 20, ncol = 5)
distr_full_5 <- diag(bfSearch(full_grid_5, start = cbind(x = 1, y = 1:nrow(full_grid_5)),
                         end = cbind(x = ncol(full_grid_5), y = 1:nrow(full_grid_5))))

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

  print(paste0("Processing grid #", j, ":"))

  pb <- txtProgressBar(max = 10000, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)

  tmp <- foreach (i = 1:10000, .combine = rbind_list, .packages = "trailClearing",
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

                    perc <- diag(bfSearch(rnd_grid, start = cbind(x = 1, y = 1:nrow(rnd_grid)),
                                     end = cbind(x = ncol(rnd_grid), y = 1:nrow(rnd_grid))))

                    data.frame(median = median(perc), mean = mean(perc),
                               min = min(perc), max = max(perc),
                               overlap = overlap(distr_full, perc))
                  }

  tmp$idx <- j

  res <- rbind_list(res, tmp)

  close(pb)
}

stopCluster(cl)








