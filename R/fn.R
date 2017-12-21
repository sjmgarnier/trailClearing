#' @export
makeGrid <- function(df, nrow = 20, ncol = 15) {
  myMat <- matrix(1, nrow = nrow * 2, ncol = ncol * 2 + 1)
  cols <- match(df$Col, letters)
  rows <- df$Row + (df$Row - 1) + (cols %% 2)
  cols <- cols + (cols - 1) + 1
  idx <- (cols - 1) * nrow * 2 + rows
  myMat[idx[df$Cleared == 0]] <- 0
  myMat
}


#' @export
makeFullGrid <- function(nrow = 20, ncol = 15) {
  myMat <- matrix(1, nrow = nrow * 2, ncol = ncol * 2 + 1)
  cols <- rep(1:ncol, each = nrow)
  rows <- 1:nrow + (1:nrow - 1) + (cols %% 2)
  cols <- cols + (cols - 1) + 1
  idx <- (cols - 1) * nrow * 2 + rows
  myMat[idx] <- 0
  myMat
}


#' @export
percolate <- function(grid, n = 1000) {
  nLanes <- nrow(grid)
  finish <- ncol(grid)
  r <- sample(1:nLanes, n, replace = TRUE)
  c <- rep(1, n)
  clock <- rep(0, n)

  done <- FALSE
  while (!done) {
    arrived <- c >= finish
    clock[!arrived] <- clock[!arrived] + 1

    idx <- c * nLanes + r
    clear <- grid[idx] == 1 & !arrived

    c[clear] <- c[clear] + 1
    r[!clear] <- r[!clear] + sample(c(-1, 1), sum(!clear), replace = TRUE)
    r[r > nLanes] <- nLanes
    r[r < 1] <- 1

    done <- sum(arrived) == n
  }

  clock
}


#' @export
percolate2 <- function(grid, n = 1000) {
  nLanes <- nrow(grid)
  finish <- ncol(grid)
  r <- sample(1:nLanes, n, replace = TRUE)
  c <- rep(1, n)
  clock <- rep(0, n)

  done <- FALSE
  while (!done) {
    arrived <- c >= finish
    clock[!arrived] <- clock[!arrived] + 1

    idx <- c * nLanes + r
    clear <- grid[idx] == 1 & !arrived

    mov <- sample(c(-1, 0, 1), sum(clear), replace = TRUE)
    r[clear] <- r[clear] + mov
    c[clear][mov == 0] <- c[clear][mov == 0] + 1

    mov <- sample(c(-1, 1), sum(!clear), replace = TRUE)
    r[!clear] <- r[!clear] + mov

    r[r > nLanes] <- nLanes - 1
    r[r < 1] <- 2

    done <- sum(arrived) == n
  }

  clock
}


#' @export
overlap <- function(d1, d2) {
  tab1 <- as.data.frame(table(d1))
  tab1$d1 <- as.numeric(as.character(tab1$d1))

  tab2 <- as.data.frame(table(d2))
  tab2$d2 <- as.numeric(as.character(tab2$d2))

  rg <- range(c(tab1$d, tab2$d))
  d <- rg[1]:rg[2]

  count1 <- count2 <- rep(0, length(d))
  count1[d %in% tab1$d1] <- tab1$Freq / length(d1)
  count2[d %in% tab2$d2] <- tab2$Freq / length(d2)

  1 - sum(abs(count1 - count2)) / 2
}


#' @export
lee <- function(m, start, end) {
  if (m[end[1], end[2]] == 0)
    stop("The endpoint is not available.")

  if (!any(m %in% c(0, 1)))
    stop("The matrix m should only have 0s and 1s.")

  path <- raster::raster(m)
  path[path == 1] <- NA
  path[path == 0] <- Inf
  path[start[1], start[2]] <- 1

  w <- matrix(c(NA, 1, NA, 1, NA, 1, NA, 1, NA), nrow = 3)
  done <- FALSE

  while (!done) {
    suppressWarnings(path <- raster::focal(path, w = w, na.rm = TRUE, NAonly = TRUE, pad = TRUE,
                                           function(...) { tmp <- min(...) + 1; ifelse(tmp == Inf, NA, tmp)}))
    done <- !is.na(path[end[1], end[2]])
  }

  as.numeric(path[end[1], end[2]])
}


#' @export
allPaths <- function(m) {
  require(foreach)

  # df <- expand.grid(start_r = 1:nrow(m), start_c = 1, end_r = 1:nrow(m), end_c = ncol(m))
  df <- data.frame(start_r = 1:nrow(m), start_c = 1, end_r = 1:nrow(m), end_c = ncol(m))

  pb <- txtProgressBar(max = nrow(df), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)

  res <- foreach(i = 1:nrow(df), .combine = c, .packages = "trailClearing",
                 .options.snow = list(progress = progress)) %dopar% {
                   lee(m, c(df$start_r[i], df$start_c[i]), c(df$end_r[i], df$end_c[i]))
                 }
  close(pb)

  df$l <- res
  df
}


sub2ind <- function(col, row, nrow) {
  (col - 1) * nrow + row;
}

#' @export
bfSearch <- function(m, start, end) {
  g <- expand.grid(y = 1:nrow(m), x = 1:ncol(m))
  g$id <- sub2ind(g$x, g$y, nrow(m))
  g$valid <- as.vector(m) == 1
  g$start <- g$id %in% sub2ind(start[, 1], start[, 2], nrow(m))
  g$end <- g$id %in% sub2ind(end[, 1], end[, 2], nrow(m))

  dd <- dist(g[, 1:2])
  ddm <- as.matrix(dd)

  dadj <- ddm == 1
  dadj[!g$valid, ] <- FALSE
  dadj[, !g$valid] <- FALSE

  starts <- which(g$start)
  ends <- which(g$end)

  dg <- igraph::graph.adjacency(dadj, mode = "undirected")

  igraph::distances(dg, starts, ends)
}
