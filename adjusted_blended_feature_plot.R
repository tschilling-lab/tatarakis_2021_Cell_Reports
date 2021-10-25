############Run full script to change blending colors for feature plot##################
############Find and replace colors - current are "dodgerblue3" and "darkgreen"########
#######################Must specify reduction and pt.size###############################
###################Final function call will be Feature_Adjusted()#######################

library(rlang)
library(cowplot)

SetQuantile <- function (cutoff, data) 
{
  if (grepl(pattern = "^q[0-9]{1,2}$", x = as.character(x = cutoff), 
            perl = TRUE)) {
    this.quantile <- as.numeric(x = sub(pattern = "q", replacement = "", 
                                        x = as.character(x = cutoff)))/100
    data <- unlist(x = data)
    data <- data[data > 0]
    cutoff <- quantile(x = data, probs = this.quantile)
  }
  return(as.numeric(x = cutoff))
}

RandomName <- function (length = 5L, ...) 
{
  CheckDots(..., fxns = "sample")
  return(paste(sample(x = letters, size = length, ...), collapse = ""))
}

CheckDots <- function (..., fxns = NULL) 
{
  args.names <- names(x = list(...))
  if (length(x = list(...)) == 0) {
    return(invisible(x = NULL))
  }
  if (is.null(x = args.names)) {
    stop("No named arguments passed")
  }
  if (length(x = fxns) == 1) {
    fxns <- list(fxns)
  }
  for (f in fxns) {
    if (!(is.character(x = f) || is.function(x = f))) {
      stop("CheckDots only works on characters or functions, not ", 
           class(x = f))
    }
  }
  fxn.args <- suppressWarnings(expr = sapply(X = fxns, FUN = function(x) {
    x <- tryCatch(expr = if (isS3stdGeneric(f = x)) {
      as.character(x = methods(generic.function = x))
    }
    else {
      x
    }, error = function(...) {
      return(x)
    })
    x <- if (is.character(x = x)) {
      sapply(X = x, FUN = argsAnywhere, simplify = FALSE, 
             USE.NAMES = TRUE)
    }
    else if (length(x = x) <= 1) {
      list(x)
    }
    return(sapply(X = x, FUN = function(f) {
      return(names(x = formals(fun = f)))
    }, simplify = FALSE, USE.NAMES = TRUE))
  }, simplify = FALSE, USE.NAMES = TRUE))
  fxn.args <- unlist(x = fxn.args, recursive = FALSE)
  fxn.null <- vapply(X = fxn.args, FUN = is.null, FUN.VALUE = logical(length = 1L))
  if (all(fxn.null) && !is.null(x = fxns)) {
    stop("None of the functions passed could be found")
  }
  else if (any(fxn.null)) {
    warning("The following functions passed could not be found: ", 
            paste(names(x = which(x = fxn.null)), collapse = ", "), 
            call. = FALSE, immediate. = TRUE)
    fxn.args <- Filter(f = Negate(f = is.null), x = fxn.args)
  }
  dfxns <- vector(mode = "logical", length = length(x = fxn.args))
  names(x = dfxns) <- names(x = fxn.args)
  for (i in 1:length(x = fxn.args)) {
    dfxns[i] <- any(grepl(pattern = "...", x = fxn.args[[i]], 
                          fixed = TRUE))
  }
  if (any(dfxns)) {
    dfxns <- names(x = which(x = dfxns))
    if (any(nchar(x = dfxns) > 0)) {
      fx <- vapply(X = Filter(f = nchar, x = dfxns), FUN = function(x) {
        if (isS3method(method = x)) {
          x <- unlist(x = strsplit(x = x, split = "\\."))
          x <- x[length(x = x) - 1L]
        }
        return(x)
      }, FUN.VALUE = character(length = 1L))
      message("The following functions and any applicable methods accept the dots: ", 
              paste(unique(x = fx), collapse = ", "))
      if (any(nchar(x = dfxns) < 1)) {
        message("In addition, there is/are ", length(x = Filter(f = Negate(f = nchar), 
                                                                x = dfxns)), " other function(s) that accept(s) the dots")
      }
    }
    else {
      message("There is/are ", length(x = dfxns), "function(s) that accept(s) the dots")
    }
  }
  else {
    unused <- Filter(f = function(x) {
      return(!x %in% unlist(x = fxn.args))
    }, x = args.names)
    if (length(x = unused) > 0) {
      msg <- paste0("The following arguments are not used: ", 
                    paste(unused, collapse = ", "))
      switch(EXPR = getOption(x = "Seurat.checkdots"), 
             warn = warning(msg, call. = FALSE, immediate. = TRUE), 
             stop = stop(msg), silent = NULL, stop("Invalid Seurat.checkdots option. Please choose one of warn, stop, silent"))
      unused.hints <- sapply(X = unused, FUN = OldParamHints)
      names(x = unused.hints) <- unused
      unused.hints <- na.omit(object = unused.hints)
      if (length(x = unused.hints) > 0) {
        message("Suggested parameter: ", paste(unused.hints, 
                                               "instead of", names(x = unused.hints), collapse = "; "), 
                "\n")
      }
    }
  }
}

BlendMatrix <- function (n = 10, col.threshold = 0.5, two.colors = c("dodgerblue3", "darkgreen"), negative.color = "white") 
{
  if (0 > col.threshold || col.threshold > 1) {
    stop("col.threshold must be between 0 and 1")
  }
  C0 <- colorRamp(colors = negative.color)(1)
  ramp <- colorRamp(colors = two.colors)
  C1 <- ramp(x = 0)
  C2 <- ramp(x = 1)
  merge.weight <- min(255/(C1 + C2 + C0 + 0.01))
  sigmoid <- function(x) {
    return(1/(1 + exp(-x)))
  }
  blend_color <- function(i, j, col.threshold, n, C0, C1, C2, 
                          merge.weight) {
    c.min <- sigmoid(5 * (1/n - col.threshold))
    c.max <- sigmoid(5 * (1 - col.threshold))
    c1_weight <- sigmoid(5 * (i/n - col.threshold))
    c2_weight <- sigmoid(5 * (j/n - col.threshold))
    c0_weight <- sigmoid(5 * ((i + j)/(2 * n) - col.threshold))
    c1_weight <- (c1_weight - c.min)/(c.max - c.min)
    c2_weight <- (c2_weight - c.min)/(c.max - c.min)
    c0_weight <- (c0_weight - c.min)/(c.max - c.min)
    C1_length <- sqrt(sum((C1 - C0)^2))
    C2_length <- sqrt(sum((C2 - C0)^2))
    C1_unit <- (C1 - C0)/C1_length
    C2_unit <- (C2 - C0)/C2_length
    C1_weight <- C1_unit * c1_weight
    C2_weight <- C2_unit * c2_weight
    C_blend <- C1_weight * (i - 1) * C1_length/(n - 1) + 
      C2_weight * (j - 1) * C2_length/(n - 1) + (i - 1) * 
      (j - 1) * c0_weight * C0/(n - 1)^2 + C0
    C_blend[C_blend > 255] <- 255
    C_blend[C_blend < 0] <- 0
    return(rgb(red = C_blend[, 1], green = C_blend[, 2], 
               blue = C_blend[, 3], alpha = 255, maxColorValue = 255))
  }
  blend_matrix <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      blend_matrix[i, j] <- blend_color(i = i, j = j, col.threshold = col.threshold, 
                                        n = n, C0 = C0, C1 = C1, C2 = C2, merge.weight = merge.weight)
    }
  }
  return(blend_matrix)
}

BlendExpression <- function (data) 
{
  if (ncol(x = data) != 2) {
    stop("'BlendExpression' only blends two features")
  }
  features <- colnames(x = data)
  data <- as.data.frame(x = apply(X = data, MARGIN = 2, FUN = function(x) {
    return(round(x = 9 * (x - min(x))/(max(x) - min(x))))
  }))
  data[, 3] <- data[, 1] + data[, 2] * 10
  colnames(x = data) <- c(features, paste(features, collapse = "_"))
  for (i in 1:ncol(x = data)) {
    data[, i] <- factor(x = data[, i])
  }
  return(data)
}

SingleDimPlot <- function (data, dims, col.by = NULL, cols = NULL, pt.size = NULL, 
          shape.by = NULL, order = NULL, label = FALSE, repel = FALSE, 
          label.size = 4, cells.highlight = NULL, cols.highlight = "#DE2D26", 
          sizes.highlight = 1, na.value = "grey50") 
{
  pt.size <- pt.size %||% AutoPointSize(data = data)
  if (length(x = dims) != 2) {
    stop("'dims' must be a two-length vector")
  }
  if (!is.data.frame(x = data)) {
    data <- as.data.frame(x = data)
  }
  if (is.character(x = dims) && !all(dims %in% colnames(x = data))) {
    stop("Cannot find dimensions to plot in data")
  }
  else if (is.numeric(x = dims)) {
    dims <- colnames(x = data)[dims]
  }
  if (!is.null(x = cells.highlight)) {
    highlight.info <- SetHighlight(cells.highlight = cells.highlight, 
                                   cells.all = rownames(x = data), sizes.highlight = sizes.highlight %||% 
                                     pt.size, cols.highlight = cols.highlight, col.base = cols[1] %||% 
                                     "#C3C3C3", pt.size = pt.size)
    order <- highlight.info$plot.order
    data$highlight <- highlight.info$highlight
    col.by <- "highlight"
    pt.size <- highlight.info$size
    cols <- highlight.info$color
  }
  if (!is.null(x = order) && !is.null(x = col.by)) {
    if (typeof(x = order) == "logical") {
      if (order) {
        data <- data[order(data[, col.by]), ]
      }
    }
    else {
      order <- rev(x = c(order, setdiff(x = unique(x = data[, 
                                                            col.by]), y = order)))
      data[, col.by] <- factor(x = data[, col.by], levels = order)
      new.order <- order(x = data[, col.by])
      data <- data[new.order, ]
      if (length(x = pt.size) == length(x = new.order)) {
        pt.size <- pt.size[new.order]
      }
    }
  }
  if (!is.null(x = col.by) && !col.by %in% colnames(x = data)) {
    warning("Cannot find ", col.by, " in plotting data, not coloring plot")
    col.by <- NULL
  }
  else {
    col.index <- match(x = col.by, table = colnames(x = data))
    if (grepl(pattern = "^\\d", x = col.by)) {
      col.by <- paste0("x", col.by)
    }
    else if (grepl(pattern = "-", x = col.by)) {
      col.by <- gsub(pattern = "-", replacement = ".", 
                     x = col.by)
    }
    colnames(x = data)[col.index] <- col.by
  }
  if (!is.null(x = shape.by) && !shape.by %in% colnames(x = data)) {
    warning("Cannot find ", shape.by, " in plotting data, not shaping plot")
  }
  plot <- ggplot(data = data) + geom_point(mapping = aes_string(x = dims[1], 
                                                                y = dims[2], color = paste0("`", col.by, "`"), shape = shape.by), 
                                           size = pt.size) + guides(color = guide_legend(override.aes = list(size = 3))) + 
    labs(color = NULL)
  if (label && !is.null(x = col.by)) {
    plot <- LabelClusters(plot = plot, id = col.by, repel = repel, 
                          size = label.size)
  }
  if (!is.null(x = cols)) {
    if (length(x = cols) == 1 && (is.numeric(x = cols) || 
                                  cols %in% rownames(x = brewer.pal.info))) {
      scale <- scale_color_brewer(palette = cols, na.value = na.value)
    }
    else if (length(x = cols) == 1 && (cols %in% c("alphabet", 
                                                   "alphabet2", "glasbey", "polychrome", "stepped"))) {
      colors <- DiscretePalette(length(unique(data[[col.by]])), 
                                palette = cols)
      scale <- scale_color_manual(values = colors, na.value = na.value)
    }
    else {
      scale <- scale_color_manual(values = cols, na.value = na.value)
    }
    plot <- plot + scale
  }
  plot <- plot + theme_cowplot()
  return(plot)
}

BlendMap <- function (color.matrix) 
{
  color.heat <- matrix(data = 1:prod(dim(x = color.matrix)) - 
                         1, nrow = nrow(x = color.matrix), ncol = ncol(x = color.matrix), 
                       dimnames = list(1:nrow(x = color.matrix), 1:ncol(x = color.matrix)))
  xbreaks <- seq.int(from = 0, to = nrow(x = color.matrix), 
                     by = 2)
  ybreaks <- seq.int(from = 0, to = ncol(x = color.matrix), 
                     by = 2)
  color.heat <- Melt(x = color.heat)
  color.heat$rows <- as.numeric(x = as.character(x = color.heat$rows))
  color.heat$cols <- as.numeric(x = as.character(x = color.heat$cols))
  color.heat$vals <- factor(x = color.heat$vals)
  plot <- ggplot(data = color.heat, mapping = aes_string(x = "rows", 
                                                         y = "cols", fill = "vals")) + geom_raster(show.legend = FALSE) + 
    theme(plot.margin = unit(x = rep.int(x = 0, times = 4), 
                             units = "cm")) + scale_x_continuous(breaks = xbreaks, 
                                                                 expand = c(0, 0), labels = xbreaks) + scale_y_continuous(breaks = ybreaks, 
                                                                                                                          expand = c(0, 0), labels = ybreaks) + scale_fill_manual(values = as.vector(x = color.matrix)) + 
    theme_cowplot()
  return(plot)
}

Melt <- function (x) 
{
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(rows = rep.int(x = rownames(x = x), times = ncol(x = x)), 
                    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, 
                                             times = nrow(x = x))), vals = unlist(x = x, use.names = FALSE)))
}


Feature_Adjusted <- function (object, features, dims = c(1, 2), cells = NULL, cols = if (blend) {
  c("gray90", "dodgerblue3", "darkgreen")
} else {
  c("lightgrey", "blue")
}, pt.size = NULL, order = FALSE, min.cutoff = NA, max.cutoff = NA, 
reduction = NULL, split.by = NULL, shape.by = NULL, slot = "data", 
blend = FALSE, blend.threshold = 0.5, label = FALSE, label.size = 4, 
repel = FALSE, ncol = NULL, combine = TRUE, coord.fixed = FALSE, 
by.col = TRUE, sort.cell = FALSE) 
{
  no.right <- theme(axis.line.y.right = element_blank(), axis.ticks.y.right = element_blank(), 
                    axis.text.y.right = element_blank(), axis.title.y.right = element_text(face = "bold", 
                                                                                           size = 14, margin = margin(r = 7)))
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (length(x = dims) != 2 || !is.numeric(x = dims)) {
    stop("'dims' must be a two-length integer vector")
  }
  if (blend && length(x = features) != 2) {
    stop("Blending feature plots only works with two features")
  }
  if (blend) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)), 
                   `0` = {
                     warning("No colors provided, using default colors", 
                             call. = FALSE, immediate. = TRUE)
                     default.colors
                   }, `1` = {
                     warning("Only one color provided, assuming specified is double-negative and augmenting with default colors", 
                             call. = FALSE, immediate. = TRUE)
                     c(cols, default.colors[2:3])
                   }, `2` = {
                     warning("Only two colors provided, assuming specified are for features and agumenting with '", 
                             default.colors[1], "' for double-negatives", 
                             call. = FALSE, immediate. = TRUE)
                     c(default.colors[1], cols)
                   }, `3` = cols, {
                     warning("More than three colors provided, using only first three", 
                             call. = FALSE, immediate. = TRUE)
                     cols[1:3]
                   })
  }
  if (blend && length(x = cols) != 3) {
    stop("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% colnames(x = object)
  data <- FetchData(object = object, vars = c(dims, "ident", 
                                              features), cells = cells, slot = slot)
  if (ncol(x = data) < 4) {
    stop("None of the requested features were found: ", paste(features, 
                                                              collapse = ", "), " in slot ", slot, call. = FALSE)
  }
  else if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found", call. = FALSE)
  }
  features <- colnames(x = data)[4:ncol(x = data)]
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[, 
                                                           feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff, 
                                              max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  brewer.gran <- ifelse(test = length(x = cols) == 1, yes = brewer.pal.info[cols, 
                                                                            ]$maxcolors, no = length(x = cols))
  data[, 4:ncol(x = data)] <- sapply(X = 4:ncol(x = data), 
                                     FUN = function(index) {
                                       data.feature <- as.vector(x = data[, index])
                                       min.use <- SetQuantile(cutoff = min.cutoff[index - 
                                                                                    3], data.feature)
                                       max.use <- SetQuantile(cutoff = max.cutoff[index - 
                                                                                    3], data.feature)
                                       data.feature[data.feature < min.use] <- min.use
                                       data.feature[data.feature > max.use] <- max.use
                                       if (brewer.gran == 2) {
                                         return(data.feature)
                                       }
                                       data.cut <- if (all(data.feature == 0)) {
                                         0
                                       }
                                       else {
                                         as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), 
                                                                          breaks = brewer.gran)))
                                       }
                                       return(data.cut)
                                     })
  colnames(x = data)[4:ncol(x = data)] <- features
  rownames(x = data) <- cells
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  }
  else {
    switch(EXPR = split.by, ident = Idents(object = object)[cells], 
           object[[split.by, drop = TRUE]][cells])
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list", length = ifelse(test = blend, 
                                                 yes = 4, no = length(x = features) * length(x = levels(x = data$split))))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, 
                                                                   dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, 
                                                               dims[2]])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(two.colors = cols[2:3], col.threshold = blend.threshold, 
                                negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ], 
                   as.vector(x = color.matrix))
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, 
                      , drop = FALSE]
    if (blend) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[, 
                                                       features]) == 0]
      if (length(x = no.expression) != 0) {
        stop("The following features have no value: ", 
             paste(no.expression, collapse = ", "), call. = FALSE)
      }
      data.plot <- cbind(data.plot[, c(dims, "ident")], 
                         BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (blend) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, 
                                                              feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      }
      else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, "ident", feature, 
                                   shape.by)]
      if (sort.cell) {
        data.single <- data.single[order(data.single[, 
                                                     feature]), ]
      }
      plot <- SingleDimPlot(data = data.single, dims = dims, 
                            col.by = feature, order = order, pt.size = pt.size, 
                            cols = cols.use, shape.by = shape.by, label = FALSE) + 
        scale_x_continuous(limits = xlims) + scale_y_continuous(limits = ylims) + 
        theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
      if (label) {
        plot <- LabelClusters(plot = plot, id = "ident", 
                              repel = repel, size = label.size)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, 
                                                         colour = "black"))
        plot <- plot + if (i == 1) {
          labs(title = feature)
        }
        else {
          labs(title = NULL)
        }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = ident)) + 
                             no.right)
        }
        if (j != 1) {
          plot <- plot + theme(axis.line.y = element_blank(), 
                               axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
                               axis.title.y.left = element_blank())
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(axis.line.x = element_blank(), 
                               axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
                               axis.title.x = element_blank())
        }
      }
      else {
        plot <- plot + labs(title = feature)
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        }
        else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warning("All cells have the same value (", 
                    unique.feature.exp, ") of ", feature, ".")
            if (unique.feature.exp == 0) {
              cols.grad <- cols[1]
            }
            else {
              cols.grad <- cols
            }
          }
          plot <- suppressMessages(expr = plot + scale_color_gradientn(colors = cols.grad, 
                                                                       guide = "colorbar"))
        }
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plot <- plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (blend) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(x = plots, 
                                              values = list(blend.legend + scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) > 
                                                                                                                                  1, yes = levels(x = data$split)[ii], no = "")), 
                                                                                              expand = c(0, 0)) + labs(x = features[1], y = features[2], 
                                                                                                                       title = if (ii == 1) {
                                                                                                                         paste("Color threshold:", blend.threshold)
                                                                                                                       } else {
                                                                                                                         NULL
                                                                                                                       }) + no.right), after = 4 * ii - 1))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (combine) {
    if (is.null(x = ncol)) {
      ncol <- 2
      if (length(x = features) == 1) {
        ncol <- 1
      }
      if (length(x = features) > 6) {
        ncol <- 3
      }
      if (length(x = features) > 9) {
        ncol <- 4
      }
    }
    ncol <- ifelse(test = is.null(x = split.by) || blend, 
                   yes = ncol, no = length(x = features))
    legend <- if (blend) {
      "none"
    }
    else {
      split.by
    }
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() + 
                                  ggtitle("") + scale_y_continuous(sec.axis = dup_axis(name = "")) + 
                                  no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 
                 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(plots[[i]] + scale_y_continuous(sec.axis = dup_axis(name = features[[idx]])) + 
                                         no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots)%%length(x = features) == 
                      1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) + 
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
      }
      plots <- plots[c(do.call(what = rbind, args = split(x = 1:length(x = plots), 
                                                          f = ceiling(x = seq_along(along.with = 1:length(x = plots))/length(x = features)))))]
      plots <- CombinePlots(plots = plots, ncol = nsplits, 
                            legend = legend)
    }
    else {
      plots <- CombinePlots(plots = plots, ncol = ncol, 
                            legend = legend, nrow = split.by)
    }
  }
  return(plots)
}