#' Load required libraries
#'
#' This script uses several plotting libraries including ggplot2, grid, signal, gridExtra, reshape, and cowplot.

library(ggplot2)
library(grid)
library(signal)
library(gridExtra)
library(reshape)
library(cowplot)


#' Format p-value for display
#'
#' Formats a p-value into a string with appropriate formatting. If the p-value is less than 0.001, it omits the equal sign.
#'
#' @param p Numeric, the p-value to format.
#' @return A character string representing the formatted p-value.
#' @examples
#' p.format(0.0005)
p.format <- function(p) {
  if (p < 0.001) {
    label <- paste0("p ", format.pval(pv = p, digits = 3, eps = 0.001))
  } else {
    label <- paste0("p = ", format.pval(pv = p, digits = 3, eps = 0.001))
  }
  return(label)
}


#' Arrange multiple ggplot objects in a grid
#'
#' This function accepts multiple ggplot objects (either via ... or a provided list) and displays them in a grid layout.
#'
#' @param ... ggplot objects passed as individual arguments.
#' @param plotlist An optional list of ggplot objects.
#' @param file (unused) parameter reserved for file output.
#' @param cols Number of columns in the layout.
#' @param layout A matrix specifying the layout. If provided, cols is ignored.
#' @param heights A numeric vector specifying the relative heights of rows.
#' @return None. The function prints the plots directly.
#' @examples
#' multiplot(p1, p2, cols = 2)
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL, heights = c(2, 1)) {
  require(grid)
  # Combine the ggplot objects
  plots <- c(list(...), plotlist)
  numPlots <- length(plots)

  # Create layout matrix if not provided
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                     ncol = cols, nrow = ceiling(numPlots / cols))
  }

  if (numPlots == 1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), heights = heights)))

    # Loop through each plot and place it in the grid
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#' Plot EDA data with optional accelerometer data
#'
#' This function creates a line plot for EDA data and, optionally, an additional plot for accelerometer data.
#'
#' @param data A data frame containing at least columns 'Timestamp' and 'EDA'. Also expected are 'X', 'Y', 'Z' if accelerometer is desired.
#' @param title Character string for the plot title. Default is 'EDA Data'.
#' @param show_acc Logical. If TRUE, an additional plot for acceleration is shown.
#' @param show_raw_acc Logical. If TRUE, raw accelerometer data (X, Y, Z) is plotted; otherwise, the magnitude (Movement) is computed and plotted.
#' @return A ggplot object or a cowplot object if accelerometer data is added.
#' @examples
#' plot.eda(mydata, title = 'Subject 1 EDA', show_acc = TRUE)
plot.eda <- function(data, title = "EDA Data", show_acc = FALSE, show_raw_acc = FALSE) {
  eda <- ggplot(data = data, aes(x = Timestamp, y = EDA)) + 
    geom_line(colour = "#2CC4FF") + 
    xlab("Time") + ylab("EDA") + ggtitle(title)

  if (show_acc) {
    if (show_raw_acc) {
      acc <- ggplot(data = data) + 
        geom_line(aes(x = Timestamp, y = X), colour = "red") + 
        geom_line(aes(x = Timestamp, y = Y), colour = "green") + 
        geom_line(aes(x = Timestamp, y = Z), colour = "blue") + 
        xlab("Time") + ylab("Acc")
    } else {
      data$Movement <- sqrt(data$X * data$X + data$Y * data$Y + data$Z * data$Z)
      acc <- ggplot(data = data) + 
        geom_line(aes(x = Timestamp, y = Movement), colour = "red") + 
        xlab("Time") + ylab("Motion")
    }
    return(plot_grid(eda, acc, ncol = 1, rel_heights = c(3, 1)))
  } else {
    return(eda)
  }
}


#' Plot dyadic data using ggplot2
#'
#' This function reshapes data for two variables (dyad) and plots them as lines with distinct colors.
#'
#' @param d A data frame with at least a 'Timestamp' column and two other variables to compare.
#' @param title Plot title.
#' @param ylabel Label for the y-axis. Default is 'EDA'.
#' @param legend_label Label for the legend. Default is 'Participant'.
#' @param grayscale Logical. If TRUE, the plot is rendered in grayscale.
#' @param ... Additional parameters passed to ggplot.
#' @return A ggplot object.
#' @examples
#' plotDyad(mydata, title = 'Dyad Plot')
plotDyad <- function(d, title = "", ylabel = "EDA", legend_label = "Participant", grayscale = FALSE, ...) {
  data <- melt(d, id.vars = c("Timestamp"), variable_name = legend_label)
  g <- ggplot(data, aes_string(x = "Timestamp", y = "value", col = legend_label)) + 
    geom_line(size = 1) + ylab(ylabel) + xlab("Time") + ggtitle(title)
  if (grayscale) {
    g <- g + scale_color_grey()
  }
  print(g)
  return(g)
}


#' Plot dyadic data using Plotly
#'
#' This function reshapes data for two variables and creates an interactive Plotly line plot.
#'
#' @param d A data frame with at least a 'Timestamp' column and two variables.
#' @param title Plot title.
#' @param ylabel Label for the y-axis. Default is 'EDA'.
#' @param legend_label Label for the legend. Default is 'Participant'.
#' @param grayscale Logical. If TRUE, applies a grayscale color scheme.
#' @param ... Additional parameters.
#' @return A Plotly plot object.
#' @examples
#' plotlyDyad(mydata, title = 'Dyad Interactive Plot')
plotlyDyad <- function(d, title = "", ylabel = "EDA", legend_label = "Participant", grayscale = FALSE, ...) {
  data <- melt(d, id.vars = c("Timestamp"), variable_name = legend_label)
  p <- plot_ly(data, x = ~Timestamp, y = ~value, color = ~get(legend_label))
  p <- add_lines(p)
  print(p)
  return(p)
}


#' Plot EDA by Condition
#'
#' This function reads EDA data and overlays condition codes as shaded rectangles on the plot.
#'
#' @param eda A source for EDA data; expected to be processed by read.eda.
#' @param codes A data frame containing condition codes with columns 'Start Time', 'End Time', and 'Condition'.
#' @param title Plot title.
#' @return A ggplot object representing the EDA plot with condition overlays.
#' @examples
#' plotEDAByCondition(eda_file, condition_codes, title = "EDA by Condition")
plotEDAByCondition <- function(eda, codes, title = "") {
  eda <- read.eda(eda)
  g1 <- ggplot(data = eda, aes(x = Timestamp, y = EDA)) + 
    geom_line(colour = "#2CC4FF") + xlab("Time") + ylab("EDA (uS)") + ggtitle(title)
  g1 <- g1 + geom_rect(data = codes, aes(xmin = `Start Time`, xmax = `End Time`, ymin = -1.0, ymax = 0, fill = Condition, col = NULL))
  print(g1)
  return(g1)
}


#' Create a histogram with density and normality annotation
#'
#' This function generates a density plot for a numeric vector. If groups are provided, it overlays densities by group.
#'
#' @param x Numeric vector to plot.
#' @param groups Optional grouping variable. If provided, densities are colored by group.
#' @return A ggplot object with the density plot and normality test annotation (if applicable).
#' @examples
#' gghist(mydata$variable)
gghist <- function(x, groups = NULL) {

  norm.test <- function(data) {
    w <- shapiro.test(data)
    if (w$p.value < 0.001) {
      p.label <- paste0("p ", format.pval(pv = w$p.value, digits = 3, eps = 0.001))
    } else {
      p.label <- paste0("p = ", format.pval(pv = w$p.value, digits = 3, eps = 0.001))
    }
    label <- paste0("Shapiro-Wilk Test of Normality\nW = ", format(w$statistic, digits = 3), ", ", p.label)
    return(label)
  }

  x.mu <- mean(x, na.rm = TRUE)
  x.sd <- sd(x, na.rm = TRUE)

  if (!is.null(groups)) {
    plt <- ggplot(data = data.frame(x = x, group = groups), aes(x = x, fill = group)) + 
      geom_density(alpha = 0.5)
  } else {
    plt <- ggplot(data = data.frame(x = x), aes(x = x)) + 
      geom_density(col = "darkgray", fill = "gray", lwd = 1) +
      stat_function(fun = dnorm, colour = "red", args = list(mean = x.mu, sd = x.sd))

    y.range <- ggplot_build(plt)$layout$panel_ranges[[1]]$y.range
    x.range <- ggplot_build(plt)$layout$panel_ranges[[1]]$x.range
    text.y <- y.range[1] + (y.range[2] - y.range[1]) * 0.1
    text.x <- x.range[1] + (x.range[2] - x.range[1]) * 0.5
    plt <- plt + annotate("label", x = text.x, y = text.y, label = try(norm.test(x)), size = 5, hjust = "center")
  }
  return(plt)
}


#' Add bracket annotations for statistical comparisons on a ggplot
#'
#' This function adds bracket annotations and p-value labels to a ggplot, based on t-tests between groups.
#'
#' @param plt A ggplot object with a grouping variable mapped to the x-axis.
#' @param pad Numeric, padding value for the annotations.
#' @param offset Numeric, unused in this implementation.
#' @param sig.cutoff Numeric, significance level threshold. Comparisons with p-values below this are annotated.
#' @param ci.func Function to compute confidence intervals, default is mean_cl_normal.
#' @return A ggplot object with added annotation segments and text.
#' @examples
#' ggbracket(my_plot)
ggbracket <- function(plt, pad = 0.1, offset = 0, sig.cutoff = 0.05, ci.func = mean_cl_normal) {
  f.test <- function(x, y) {
    mdl <- lm(y ~ x)
    s <- summary(aov(mdl))
    f_val <- round(s[[1]][[4]][1], digits = 3)
    p.val <- s[[1]][[5]][1]
    if (p.val < 0.001) {
      p.label <- paste0("p ", format.pval(pv = p.val, digits = 3, eps = 0.001))
    } else {
      p.label <- paste0("p = ", format.pval(pv = p.val, digits = 3, eps = 0.001))
    }
    df <- paste(s[[1]][[1]], collapse = ",")
    str <- paste0("F(", df, ") = ", f_val, ", ", p.label)
    return(str)
  }

  y.range <- ggplot_build(plt)$layout$panel_ranges[[1]]$y.range
  pad <- 0.04 * (y.range[2] - y.range[1])

  x_vals <- plt$data[[as.character(plt$mapping$x)]]
  y_vals <- plt$data[[as.character(plt$mapping$y)]]
  groups <- levels(as.factor(x_vals))

  x_vals <- simplify2array(x_vals)
  y_vals <- simplify2array(y_vals)

  n <- 0
  g <- plt
  comparisons <- combn(groups, m = 2, simplify = FALSE)

  g <- g + annotate("text", x = groups[2], y = max(y.range) + pad * 10 * length(comparisons),
                    label = f.test(x_vals, y_vals), hjust = "center")

  for (p in rev(comparisons)) {
    group.1 <- p[1]
    group.1.idx <- which(groups == group.1)
    group.2 <- p[2]
    group.2.idx <- which(groups == group.2)
    y.g1 <- subset(y_vals, x_vals == group.1)
    y.g2 <- subset(y_vals, x_vals == group.2)

    t_result <- t.test(y.g1, y.g2)

    if (t_result$p.value < sig.cutoff) {
      if (t_result$p.value < 0.001) {
        label <- paste0("p ", format.pval(pv = t_result$p.value, digits = 3, eps = 0.001))
      } else {
        label <- paste0("p = ", format.pval(pv = t_result$p.value, digits = 3, eps = 0.001))
      }
      y.1 <- ci.func(y.g1, conf.int = (1 - sig.cutoff))$ymax + pad
      y.2 <- ci.func(y.g2, conf.int = (1 - sig.cutoff))$ymax + pad
      y.max <- max(tapply(y_vals, x_vals, FUN = function(a) { return(ci.func(a)$ymax) })) + pad * (n + 1) * 4
      df <- data.frame(x1 = c(group.1, group.1, group.2),
                       x2 = c(group.1, group.2, group.2),
                       y1 = c(y.max - pad, y.max, y.max),
                       y2 = c(y.max, y.max, y.max - pad))
      g <- g + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = df, col = "black")
      x.label <- (max(c(group.1.idx, group.2.idx)) - min(c(group.1.idx, group.2.idx))) / 2.0 + min(c(group.1.idx, group.2.idx))
      g <- g + annotate("text", x = x.label, y = y.max + pad * 2, label = label)
      n <- n + 1
    }
  }
  return(g)
}


#' Plot a metric with grouping and error bars
#'
#' This function generates a summary plot for a given metric partitioned by a group, including mean points and error bars.
#'
#' @param data A data frame containing the data to plot.
#' @param metric Character string naming the metric column to plot.
#' @param group Character string naming the grouping variable.
#' @param sig.cutoff Numeric significance cutoff for brackets. Default is 0.05.
#' @return A ggplot object with the metric plot and statistical brackets.
#' @examples
#' plotMetric(mydata, 'response', 'group')
plotMetric <- function(data, metric, group, sig.cutoff = 0.05) {
  n_fun <- function(x) {
    return(data.frame(y = mean(x, na.rm = TRUE) + se(x), label = paste0("n = ", length(x))))
  }
  
  p <- ggplot(data = data, aes_string(x = group, y = metric, col = group)) + 
    stat_summary(fun.y = mean, geom = "point") + 
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar", width = 0.25, fun.args = list(conf.int = (1 - sig.cutoff)))
  p <- ggbracket(p, sig.cutoff = sig.cutoff, ci.func = mean_cl_normal)
  return(p)
}


#' Calculate a relative x-coordinate in a ggplot
#'
#' Given a ggplot object, this function returns an x-coordinate at a relative position between the plot limits.
#'
#' @param plt A ggplot object. Defaults to the last plot.
#' @param x A numeric value between 0 and 1 representing the relative position. Default is 0.5.
#' @return A numeric value representing the x-coordinate.
#' @examples
#' ggrel.x()
ggrel.x <- function(plt = last_plot(), x = 0.5) {
  x.range <- ggplot_build(plt)$layout$panel_params[[1]]$x.range
  coord.x <- min(x.range) + (max(x.range) - min(x.range)) * x
  return(coord.x)
}


#' Calculate a relative y-coordinate in a ggplot
#'
#' Given a ggplot object, this function returns a y-coordinate at a relative position between the plot limits.
#'
#' @param plt A ggplot object. Defaults to the last plot.
#' @param y A numeric value between 0 and 1 representing the relative position. Default is 0.5.
#' @return A numeric value representing the y-coordinate.
#' @examples
#' ggrel.y()
ggrel.y <- function(plt = last_plot(), y = 0.5) {
  y.range <- ggplot_build(plt)$layout$panel_params[[1]]$y.range
  coord.y <- min(y.range) + (max(y.range) - min(y.range)) * y
  return(coord.y)
}


#' Annotate a ggplot with relative coordinates
#'
#' This function adds an annotation (default text) at a relative position on the ggplot.
#'
#' @param plt A ggplot object. Defaults to the last plot.
#' @param geom Character string specifying the type of geom to add; defaults to 'text'.
#' @param x Relative x-position (0 to 1). Default is 0.5.
#' @param y Relative y-position (0 to 1). Default is 0.1.
#' @param ... Additional arguments passed to the annotate function.
#' @return The annotation layer added to the plot.
#' @examples
#' annotate.relative(plt, label = 'Hello')
annotate.relative <- function(plt = last_plot(), geom = "text", x = 0.5, y = 0.1, ...) {
  x.range <- ggplot_build(plt)$layout$panel_params[[1]]$x.range
  y.range <- ggplot_build(plt)$layout$panel_params[[1]]$y.range
  coord.x <- min(x.range) + (max(x.range) - min(x.range)) * x
  coord.y <- min(y.range) + (max(y.range) - min(y.range)) * y
  return(annotate(geom = geom, x = coord.x, y = coord.y, ...))
}


#' Create a scatter plot with a linear fit and correlation annotation
#'
#' This function generates a scatter plot for two variables, adds a linear regression line, and annotates the plot with the correlation coefficient and p-value.
#'
#' @param x Character string representing the name of the x variable in the data frame.
#' @param y Character string representing the name of the y variable in the data frame.
#' @param data A data frame containing the variables.
#' @return A ggplot object with the scatter plot and annotations.
#' @examples
#' scatter_plot('height', 'weight', mydata)
scatter_plot <- function(x, y, data) {
  cor_fun <- function(x, y, method = "pearson") {
    r <- NULL
    try((r <- cor.test(x, y, method = method)))
    if (is.null(r)) {
      return("r = 0, p = n.s.")
    } else {
      return(paste0("r = ", round(r$estimate, 2), ", ", p.format(r$p.value)))
    }
  }
  
  p <- ggplot(data, aes_string(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE)
  cor.all <- cor_fun(data[[x]], data[[y]])
  p <- p + annotate.relative(geom = "label", plt = p, y = 0.95, x = 0.9, label = cor.all, hjust = "right")
  return(p)
}
