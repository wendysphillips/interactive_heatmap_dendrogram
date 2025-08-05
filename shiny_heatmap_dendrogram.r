# Load necessary libraries
library(shiny)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggdendro)
library(patchwork)

# Load your data here - replace with your actual data loading
# Assumes data has the items to be clustered as rows and conditions as columns
# This example data comes from the publication: https://doi.org/10.1126/sciadv.abf5733
df <- read.table("tre_data.tsv", header = TRUE, row.names = 1)

# Generate colors for clusters
generate_cluster_colors <- function() {
  z <- colors()
  not_to_use <- c(
    "gray", "grey", "white", "azure", "snow", "cornsilk", "beige", "aliceblue",
    "mintcream", "bisque", "gainsboro", "honeydew", "ivory", "lavender",
    "lemonchiffon", "lightcyan", "linen", "mistyrose", "yellow", "oldlace",
    "papayawhip", "moccasin", "palegoldenrod", "peachpuff", "seashell",
    "wheat", "thistle"
  )

  for (cl in not_to_use) {
    z <- z[grep(cl, z, invert = TRUE)]
  }

  set.seed(111)
  sample(z, length(z), replace = FALSE)
}

random_colors <- generate_cluster_colors()

# Function to create a zoomable heatmap app
create_zoom_heatmap <- function() {
  # Define UI and server components
  ui <- fluidPage(
    titlePanel("Zoomable Heatmap with Dendrogram"),
    sidebarLayout(
      sidebarPanel(
        selectInput("distance_method", "Distance Method:",
          choices = c("manhattan", "euclidean", "maximum", "binary"),
          selected = "euclidean"
        ),
        selectInput("cluster_method", "Clustering Method:",
          choices = c("complete", "single", "average", "ward.D2"),
          selected = "ward.D2"
        ),
        br(),
        h4("Cluster Controls"),
        numericInput("cut_height", "Cut Height for Clusters:",
          value = 10, min = 0, max = 100, step = 1
        ),
        br(),
        h4("Y-axis Zoom Controls"),
        numericInput("y_min", "Y-axis minimum:",
          value = NULL,
          min = 1,
          step = 1
        ),
        numericInput("y_max", "Y-axis maximum:",
          value = NULL,
          min = 1,
          step = 1
        ),
        actionButton("apply_zoom", "Apply Y-axis Zoom", class = "btn-primary"),
        br(), br(),
        textInput("search_label", "Search for gene:",
          placeholder = "Enter gene name"
        ),
        actionButton("zoom_to_label", "Zoom to Gene", class = "btn-primary"),
        br(), br(),
        actionButton("reset_zoom", "Reset Zoom", class = "btn-secondary"),
        br(), br(),
        downloadButton("download_clusters", "Export Cluster Data", class = "btn-success"),
        br(), br(),
        downloadButton("download_image", "Export Image (PNG)", class = "btn-info"),
        width = 3
      ),
      mainPanel(
        plotOutput("heatmap_plot",
          height = "800px"
        ),
        width = 9
      )
    )
  )

  server <- function(input, output, session) {
    # Reactive values for zoom
    ranges <- reactiveValues(x = NULL, y = NULL)

    # Update numeric input limits based on data
    observe({
      plot_info <- plot_data()
      nrows <- plot_info$nrows

      updateNumericInput(session, "y_min",
        value = if (is.null(ranges$y)) 1 else ranges$y[1],
        min = 1, max = nrows
      )
      updateNumericInput(session, "y_max",
        value = if (is.null(ranges$y)) nrows else ranges$y[2],
        min = 1, max = nrows
      )
    })

    # Reactive clustering calculations
    clustering_data <- reactive({
      # Convert to matrix if not already
      if (!is.matrix(df)) {
        df_matrix <- as.matrix(df)
      } else {
        df_matrix <- df
      }

      # Calculate distance matrices and hierarchical clustering
      row_dist <- dist(df_matrix, method = input$distance_method)
      col_dist <- dist(t(df_matrix), method = input$distance_method)
      row_hclust <- hclust(row_dist, method = input$cluster_method)
      col_hclust <- hclust(col_dist, method = input$cluster_method)

      list(
        df = df_matrix,
        row_hclust = row_hclust,
        col_hclust = col_hclust
      )
    })

    # Reactive plot data
    plot_data <- reactive({
      cluster_data <- clustering_data()
      df_matrix <- cluster_data$df

      # Get the row and column orders from clustering
      row_order <- cluster_data$row_hclust$order
      col_order <- cluster_data$col_hclust$order

      # Create dendrogram data
      df_dendrogram <- as.dendrogram(cluster_data$row_hclust)
      dend_data <- dendro_data(df_dendrogram)

      # Create long format data
      df_long <- df_matrix |>
        as.data.frame() |>
        rownames_to_column("gene") |>
        pivot_longer(2:(ncol(df_matrix) + 1))

      df_long$gene <- factor(df_long$gene, levels = rownames(df_matrix)[row_order])
      df_long$name <- factor(df_long$name, levels = colnames(df_matrix)[col_order])

      # Join with dendrogram positions - ensure all genes are matched
      df_long2 <- df_long |> rename(label = gene)
      df_long2 <- left_join(df_long2, dend_data$labels, by = "label")

      # Check for missing matches and handle them
      if (any(is.na(df_long2$x))) {
        warning("Some genes not found in dendrogram labels")
        # Remove rows with missing dendrogram positions
        df_long2 <- df_long2 |> filter(!is.na(x))
      }

      # Check for missing values in the data
      na_count <- sum(is.na(df_long2$value))
      if (na_count > 0) {
        message(sprintf("Found %d missing values in heatmap data", na_count))
      }

      # Invert y-coordinates so that row 1 is at the top (as users expect)
      df_long2 <- df_long2 |>
        mutate(row_order = as.integer(nrow(df_matrix) + 1 - x)) |>
        select(-x)

      list(
        df_long = df_long2,
        dend_data = dend_data,
        nrows = nrow(df_matrix)
      )
    })

    output$heatmap_plot <- renderPlot({
      plot_info <- plot_data()
      df_long2 <- plot_info$df_long
      dend_data <- plot_info$dend_data
      nrows <- plot_info$nrows

      # Create dendrogram plot with inverted y-coordinates and cluster coloring
      dend_data_inverted <- dend_data
      # Force dendrogram coordinates to integers to match heatmap
      dend_data_inverted$segments$x <- as.integer(nrows + 1 - dend_data_inverted$segments$x)
      dend_data_inverted$segments$xend <- as.integer(nrows + 1 - dend_data_inverted$segments$xend)
      dend_data_inverted$labels$x <- as.integer(nrows + 1 - dend_data_inverted$labels$x)

      # Get cluster assignments at specified height
      clusters <- cutree(clustering_data()$row_hclust, h = input$cut_height)

      # Add cluster colors to labels
      dendro_labels <- dend_data_inverted$labels
      dendro_labels$cluster <- clusters[dendro_labels$label]
      dendro_labels$cluster <- as.factor(dendro_labels$cluster)

      # Find all segments that cross the cut height (create cluster boundaries)
      cluster_segments <- dend_data_inverted$segments |>
        filter((y <= input$cut_height & yend > input$cut_height) |
          (y > input$cut_height & yend <= input$cut_height)) |>
        mutate(cut_x = ifelse(y <= input$cut_height, x, xend)) |>
        arrange(cut_x)

      # Map intersection points to actual cluster numbers
      if (nrow(cluster_segments) > 0) {
        # For each intersection point, find the cluster it represents
        cluster_segments$cluster_num <- sapply(cluster_segments$cut_x, function(x_pos) {
          # Find the closest leaf to this x position
          closest_leaf <- dendro_labels[which.min(abs(dendro_labels$x - x_pos)), ]
          return(as.numeric(as.character(closest_leaf$cluster)))
        })
      }

      p_dendro <- ggplot() +
        geom_segment(
          data = dend_data_inverted$segments,
          aes(x = y, y = x, xend = yend, yend = xend)
        ) +
        geom_text(
          data = dendro_labels,
          aes(x = y, y = x, label = label, color = cluster),
          hjust = 1.05, vjust = 0.5, angle = 0, size = 3, fontface = "bold", family = "sans"
        ) +
        geom_vline(
          xintercept = input$cut_height,
          linetype = "dashed",
          color = "red",
          alpha = 0.7
        ) +
        theme_dendro() +
        theme(
          axis.text.y = element_blank(),
          axis.text.x = element_text(color = "black", size = 10),
          axis.title.x = element_text(color = "black", size = 12),
          plot.margin = unit(c(0, 0, 15, 0), "mm"),
          legend.position = "none"
        ) +
        scale_y_continuous(limits = c(0.5, nrows + 0.5), expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_color_manual(values = random_colors) +
        labs(x = "Cut Height") +
        coord_cartesian(clip = "off")

      # Add cluster labels at intersection points
      if (nrow(cluster_segments) > 0) {
        p_dendro <- p_dendro + geom_text(
          data = cluster_segments,
          aes(y = cut_x, x = input$cut_height, label = cluster_num),
          size = 4, fontface = "bold", color = "red",
          hjust = 0.5, vjust = -0.5
        )
      }

      # Create heatmap plot with conditional y-axis guide numbers
      p_heatmap <- ggplot(df_long2, aes(x = name, y = row_order, fill = value)) +
        geom_tile(color = "white", size = 0.1) + # Add borders to see missing tiles
        scale_y_continuous(
          limits = c(0.5, nrows + 0.5),
          expand = c(0, 0),
          breaks = if (!is.null(ranges$y)) {
            NULL
          } else {
            function(x) {
              # Only show breaks in full view, not when zoomed
              range_size <- diff(x)
              if (range_size <= 50) {
                # For small ranges, show every 5th or 10th number
                step <- if (range_size <= 20) 2 else 5
                seq(from = ceiling(x[1]), to = floor(x[2]), by = step)
              } else if (range_size <= 200) {
                # For medium ranges, show every 10th or 20th
                step <- if (range_size <= 100) 10 else 20
                seq(from = ceiling(x[1] / step) * step, to = floor(x[2] / step) * step, by = step)
              } else {
                # For large ranges, show every 50th or 100th
                step <- if (range_size <= 1000) 50 else 100
                seq(from = ceiling(x[1] / step) * step, to = floor(x[2] / step) * step, by = step)
              }
            }
          },
          labels = if (!is.null(ranges$y)) {
            NULL
          } else {
            function(breaks) {
              # Only show labels in full view, not when zoomed
              user_coords <- nrows + 1 - breaks
              as.character(user_coords)
            }
          }
        ) +
        coord_cartesian(clip = "off") +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey90") +
        theme_void() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = if (!is.null(ranges$y)) element_blank() else element_text(size = 8, hjust = 1, color = "black"),
          axis.title.y = if (!is.null(ranges$y)) element_blank() else element_text(angle = 90, vjust = 0.5),
          legend.position = "left",
          plot.margin = margin(
            t = 0,
            r = 30,
            b = 0,
            l = if (!is.null(ranges$y)) 0 else 20, unit = "mm"
          ) # Adjust margin based on zoom
        ) +
        labs(fill = "Value", y = if (!is.null(ranges$y)) NULL else "Row")

      # Apply zoom if ranges are set
      if (!is.null(ranges$x)) {
        p_heatmap <- p_heatmap + xlim(ranges$x[1], ranges$x[2])
        p_dendro <- p_dendro + xlim(ranges$x[1], ranges$x[2])
      }
      if (!is.null(ranges$y)) {
        # Convert user coordinate ranges to internal coordinate ranges
        # User coordinates: 1 = top row, nrows = bottom row
        # Internal coordinates: 1 = bottom row, nrows = top row
        # So user coordinate Y maps to internal coordinate (nrows + 1 - Y)
        y_min_internal <- nrows + 1 - ranges$y[2] # User's max becomes internal min
        y_max_internal <- nrows + 1 - ranges$y[1] # User's min becomes internal max

        # Debug: count how many genes are in this range
        genes_in_range <- df_long2 |>
          filter(row_order >= ranges$y[1] & row_order <= ranges$y[2]) |>
          select(label, row_order, value) |>
          distinct(label, row_order)

        missing_values_in_range <- df_long2 |>
          filter(row_order >= ranges$y[1] & row_order <= ranges$y[2] & is.na(value)) |>
          nrow()

        message(sprintf(
          "Zoom range %d-%d: %d unique genes, %d missing values",
          ranges$y[1], ranges$y[2],
          nrow(genes_in_range), missing_values_in_range
        ))

        p_heatmap <- p_heatmap + ylim(y_min_internal - 0.5, y_max_internal + 0.5)
        p_dendro <- p_dendro + ylim(y_min_internal - 0.5, y_max_internal + 0.5)
      }

      # Combine plots - adjusted widths to account for y-axis labels
      wrap_plots(p_heatmap, p_dendro, nrow = 1, widths = c(1, 2))
    })

    # Apply Y-axis zoom based on numeric inputs
    observeEvent(input$apply_zoom, {
      if (!is.null(input$y_min) && !is.null(input$y_max)) {
        plot_info <- plot_data()
        nrows <- plot_info$nrows

        # Validate ranges
        y_min <- max(1, min(input$y_min, nrows))
        y_max <- min(nrows, max(input$y_max, 1))

        if (y_min < y_max) {
          ranges$y <- c(y_min, y_max)

          # Debug information
          message(sprintf(
            "Zooming to rows %d-%d (total: %d rows)",
            y_min, y_max, y_max - y_min + 1
          ))
        } else {
          showNotification("Y-axis minimum must be less than maximum!",
            type = "warning", duration = 3000
          )
        }
      }
    })

    # Function to zoom to a specific label
    observeEvent(input$zoom_to_label, {
      if (input$search_label != "") {
        plot_info <- plot_data()
        df_long2 <- plot_info$df_long

        # Find the gene in our processed data (which has the correct inverted coordinates)
        gene_match <- df_long2[df_long2$label == input$search_label, ]

        if (nrow(gene_match) > 0) {
          # Get the inverted y position directly from our processed data
          gene_y_pos <- gene_match$row_order[1] # This is already in our inverted coordinate system
          nrows <- plot_info$nrows

          # Create zoom window around the gene
          zoom_width_y <- 50 # Show 50 genes around the target gene

          y_min <- max(1, gene_y_pos - zoom_width_y / 2)
          y_max <- min(nrows, gene_y_pos + zoom_width_y / 2)

          ranges$y <- c(y_min, y_max)

          # Update the numeric inputs to reflect the zoom
          updateNumericInput(session, "y_min", value = y_min)
          updateNumericInput(session, "y_max", value = y_max)
        } else {
          showNotification("Gene not found!", type = "warning", duration = 3000)
        }
      }
    })

    # Reset zoom button
    observeEvent(input$reset_zoom, {
      ranges$x <- NULL
      ranges$y <- NULL

      # Reset numeric inputs
      plot_info <- plot_data()
      nrows <- plot_info$nrows
      updateNumericInput(session, "y_min", value = 1)
      updateNumericInput(session, "y_max", value = nrows)
    })

    # Download handler for cluster data
    output$download_clusters <- downloadHandler(
      filename = function() {
        paste0("cluster_assignments_", input$distance_method, "_", input$cluster_method, "_h", input$cut_height, ".csv")
      },
      content = function(file) {
        # Get cluster assignments at specified height
        clusters <- cutree(clustering_data()$row_hclust, h = input$cut_height)

        # Create data frame with row labels and cluster assignments
        cluster_data <- data.frame(
          gene_name = names(clusters),
          cluster = clusters,
          stringsAsFactors = FALSE
        )

        write.csv(cluster_data, file, row.names = FALSE)
      }
    )

    # Download handler for image export
    output$download_image <- downloadHandler(
      filename = function() {
        zoom_suffix <- if (!is.null(ranges$y)) paste0("_zoom_", ranges$y[1], "-", ranges$y[2]) else "_full"
        paste0("heatmap_", input$distance_method, "_", input$cluster_method, "_h", input$cut_height, zoom_suffix, ".png")
      },
      content = function(file) {
        # Generate the same plot as displayed
        plot_info <- plot_data()
        df_long2 <- plot_info$df_long
        dend_data <- plot_info$dend_data
        nrows <- plot_info$nrows

        # Create dendrogram plot with inverted y-coordinates and cluster coloring
        dend_data_inverted <- dend_data
        # Force dendrogram coordinates to integers to match heatmap
        dend_data_inverted$segments$x <- as.integer(nrows + 1 - dend_data_inverted$segments$x)
        dend_data_inverted$segments$xend <- as.integer(nrows + 1 - dend_data_inverted$segments$xend)
        dend_data_inverted$labels$x <- as.integer(nrows + 1 - dend_data_inverted$labels$x)

        # Get cluster assignments at specified height
        clusters <- cutree(clustering_data()$row_hclust, h = input$cut_height)

        # Add cluster colors to labels
        dendro_labels <- dend_data_inverted$labels
        dendro_labels$cluster <- clusters[dendro_labels$label]
        dendro_labels$cluster <- as.factor(dendro_labels$cluster)

        # Find all segments that cross the cut height (create cluster boundaries)
        cluster_segments <- dend_data_inverted$segments |>
          filter((y <= input$cut_height & yend > input$cut_height) |
            (y > input$cut_height & yend <= input$cut_height)) |>
          mutate(cut_x = ifelse(y <= input$cut_height, x, xend)) |>
          arrange(cut_x)

        # Map intersection points to actual cluster numbers
        if (nrow(cluster_segments) > 0) {
          # For each intersection point, find the cluster it represents
          cluster_segments$cluster_num <- sapply(cluster_segments$cut_x, function(x_pos) {
            # Find the closest leaf to this x position
            closest_leaf <- dendro_labels[which.min(abs(dendro_labels$x - x_pos)), ]
            return(as.numeric(as.character(closest_leaf$cluster)))
          })
        }

        p_dendro <- ggplot() +
          geom_segment(
            data = dend_data_inverted$segments,
            aes(x = y, y = x, xend = yend, yend = xend)
          ) +
          geom_text(
            data = dendro_labels,
            aes(x = y, y = x, label = label, color = cluster),
            hjust = 1.05, vjust = 0.5, angle = 0, size = 3, fontface = "bold", family = "sans"
          ) +
          geom_vline(
            xintercept = input$cut_height,
            linetype = "dashed",
            color = "red",
            alpha = 0.7
          ) +
          theme_dendro() +
          theme(
            axis.text.y = element_blank(),
            axis.text.x = element_text(color = "black", size = 10),
            axis.title.x = element_text(color = "black", size = 12),
            plot.margin = unit(c(0, 0, 15, 0), "mm"),
            legend.position = "none"
          ) +
          scale_y_continuous(limits = c(0.5, nrows + 0.5), expand = c(0, 0)) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_color_manual(values = random_colors) +
          labs(x = "Cut Height") +
          coord_cartesian(clip = "off")

        # Add cluster labels at intersection points
        if (nrow(cluster_segments) > 0) {
          p_dendro <- p_dendro + geom_text(
            data = cluster_segments,
            aes(y = cut_x, x = input$cut_height, label = cluster_num),
            size = 4, fontface = "bold", color = "red",
            hjust = 0.5, vjust = -0.5
          )
        }

        # Create heatmap plot with conditional y-axis guide numbers
        p_heatmap <- ggplot(df_long2, aes(x = name, y = row_order, fill = value)) +
          geom_tile(color = "white", size = 0.1) + # Add borders to see missing tiles
          scale_y_continuous(
            limits = c(0.5, nrows + 0.5),
            expand = c(0, 0),
            breaks = if (!is.null(ranges$y)) {
              NULL
            } else {
              function(x) {
                # Only show breaks in full view, not when zoomed
                range_size <- diff(x)
                if (range_size <= 50) {
                  # For small ranges, show every 5th or 10th number
                  step <- if (range_size <= 20) 2 else 5
                  seq(from = ceiling(x[1]), to = floor(x[2]), by = step)
                } else if (range_size <= 200) {
                  # For medium ranges, show every 10th or 20th
                  step <- if (range_size <= 100) 10 else 20
                  seq(from = ceiling(x[1] / step) * step, to = floor(x[2] / step) * step, by = step)
                } else {
                  # For large ranges, show every 50th or 100th
                  step <- if (range_size <= 1000) 50 else 100
                  seq(from = ceiling(x[1] / step) * step, to = floor(x[2] / step) * step, by = step)
                }
              }
            },
            labels = if (!is.null(ranges$y)) {
              NULL
            } else {
              function(breaks) {
                # Only show labels in full view, not when zoomed
                user_coords <- nrows + 1 - breaks
                as.character(user_coords)
              }
            }
          ) +
          coord_cartesian(clip = "off") +
          scale_fill_gradient2(
            low = "blue",
            mid = "white",
            high = "red",
            midpoint = 0,
            na.value = "grey90"
          ) +
          theme_void() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.text.y = if (!is.null(ranges$y)) {
              element_blank()
            } else {
              element_text(
                size = 8,
                hjust = 1,
                color = "black"
              )
            },
            axis.title.y = if (!is.null(ranges$y)) {
              element_blank()
            } else {
              element_text(
                angle = 90,
                vjust = 0.5
              )
            },
            legend.position = "left",
            plot.margin = margin(
              t = 0, r = 30, b = 0, l = if (!is.null(ranges$y)) 0 else 20,
              unit = "mm"
            ) # Adjust margin based on zoom
          ) +
          labs(fill = "Value", y = if (!is.null(ranges$y)) NULL else "Row")

        # Apply zoom if ranges are set
        if (!is.null(ranges$x)) {
          p_heatmap <- p_heatmap + xlim(ranges$x[1], ranges$x[2])
          p_dendro <- p_dendro + xlim(ranges$x[1], ranges$x[2])
        }
        if (!is.null(ranges$y)) {
          # Convert user coordinate ranges to internal coordinate ranges
          y_min_internal <- nrows + 1 - ranges$y[2] # User's max becomes internal min
          y_max_internal <- nrows + 1 - ranges$y[1] # User's min becomes internal max

          p_heatmap <- p_heatmap + ylim(y_min_internal - 0.5, y_max_internal + 0.5)
          p_dendro <- p_dendro + ylim(y_min_internal - 0.5, y_max_internal + 0.5)
        }

        # Combine plots - adjusted widths to account for y-axis labels
        combined_plot <- wrap_plots(p_heatmap, p_dendro, nrow = 1, widths = c(1, 2))

        # Save the plot as PNG
        ggsave(file, combined_plot, width = 16, height = 10, dpi = 300, units = "in")
      }
    )
  }

  shinyApp(ui, server)
}

# Run the app
create_zoom_heatmap()
