####################### SET UP ENV #######################   

library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(data.table)
library(Biostrings)
library(readr)
library(stringr)

setwd("OneDrive - KU Leuven/FWO/organellarMetagenomics/4_QualityAssembly/")

######################### FUNCTIONS #########################


################# READ IN DEPTH

read_depth_percent_id <- function(sample_dir, ID = TRUE) {
  if (ID) {
    # Locate depth files that match the pattern "mapped_pp_depth.txt"
    depth_files <- list.files(
      path = sample_dir,
      pattern = "_depth.txt$",
      full.names = TRUE
    )
    
    # Extract percent identities from file names
    percent_id <- sub(".*_id(\\d+)_depth.txt", "id\\1", basename(depth_files))
    
    # Ensure files are sorted by percent identity
    depth_files <- depth_files[order(as.numeric(percent_id))]
    percent_id <- percent_id[order(as.numeric(percent_id))]
    
    # Initialize an empty list to store depth data
    depth_list <- list()
    
    # Read and combine all depth data into a single dataframe
    for (i in seq_along(depth_files)) {
      file <- depth_files[i]
      depth_data <- tryCatch({
        # Read data and coerce to numeric while handling potential issues
        fread(file, select = c("V2", "V3"), col.names = c("Position", percent_id[i]))
      }, error = function(e) {
        warning(paste("Failed to read file:", file, "-", e$message))
        return(NULL)
      })
      
      # Handle non-numeric or missing values
      if (!is.null(depth_data)) {
        depth_data <- depth_data %>%
          mutate(
            Position = as.numeric(Position),
            across(everything(), ~ as.numeric(.))
          )
        depth_data <- depth_data[complete.cases(depth_data), ] # Remove rows with NA
        depth_list[[i]] <- depth_data
      }
    }
    
    # Combine all depth data
    combined_depth_data <- Reduce(function(x, y) merge(x, y, by = "Position", all = TRUE), depth_list)
    
  } else {
    # Read the single file without ID in filename
    depth_file <- list.files(
      path = sample_dir,
      pattern = "_depth.txt$",
      full.names = TRUE
    )
    
    if (length(depth_file) == 1) {
      combined_depth_data <- tryCatch({
        fread(depth_file, select = c("V2", "V3"), col.names = c("Position", "Depth"))
      }, error = function(e) {
        stop(paste("Failed to read file:", depth_file, "-", e$message))
      })
      
      # Ensure data is numeric and clean
      combined_depth_data <- combined_depth_data %>%
        mutate(
          Position = as.numeric(Position),
          Depth = as.numeric(Depth)
        ) %>%
        filter(complete.cases(.))
      
    } else {
      stop("Multiple or no files found for single-depth mode.")
    }
  }
  
  combined_depth_data <- combined_depth_data %>%
    rename_with(
      ~ gsub("denovo_(id\\d+)_depth\\.txt", "\\1", .x), 
      .cols = 2:ncol(combined_depth_data)
    )
  
  # Return the combined dataframe
  return(combined_depth_data)
}


################# READ IN VCF


# Function to process VCF files
read_vcf_data <- function(sample_dir, ID = TRUE) {
  if (ID) {
    # Locate VCF files with percent IDs in their filenames
    vcf_files <- list.files(
      path = sample_dir,
      pattern = "filtered.vcf.gz",
      full.names = TRUE
    )
    
    # Extract percent identities from file names
    percent_id <- sub(".*_id(\\d+)_filtered.vcf.gz", "id\\1", basename(vcf_files))
    
    # Ensure files are sorted by percent identity
    vcf_files <- vcf_files[order(as.numeric(percent_id))]
    percent_id <- percent_id[order(as.numeric(percent_id))]
    
    # Function to process individual VCF file and count variants by position
    process_vcf <- function(file) {
      if (!file.exists(file) || file.info(file)$size == 0) {
        return(NULL)
      }
      tryCatch({
        vcf_data <- fread(
          cmd = paste("gzip -cd", file),
          skip = "#CHROM",
          select = c(2),  # Select only the Position column
          col.names = c("Position")
        )
        # Count variants per position
        vcf_counts <- vcf_data %>%
          group_by(Position) %>%
          summarise(Count = n(), .groups = "drop") %>%
          arrange(Position)
        return(vcf_counts)
      }, error = function(e) {
        warning(paste("Failed to process file:", file, "-", e$message))
        return(NULL)
      })
    }
    
    # Process all VCF files and assign proper column names
    all_data <- lapply(seq_along(vcf_files), function(i) {
      vcf_counts <- process_vcf(vcf_files[i])
      if (!is.null(vcf_counts)) {
        colnames(vcf_counts)[2] <- percent_id[i]  # Rename Count column with percent ID
      }
      return(vcf_counts)
    })
    
    # Merge all results by Position
    final_data <- Reduce(function(x, y) full_join(x, y, by = "Position"), all_data)
    
  } else {
    # Handle a single VCF file without percent IDs
    vcf_file <- list.files(
      path = sample_dir,
      pattern = "_mapped.vcf.gz",
      full.names = TRUE
    )
    
    if (length(vcf_file) == 1) {
      final_data <- tryCatch({
        vcf_data <- fread(
          cmd = paste("gzip -cd", vcf_file),
          skip = "#CHROM",
          select = c(2),  # Select only the Position column
          col.names = c("Position")
        )
        # Count variants per position
        vcf_data %>%
          group_by(Position) %>%
          summarise(Count = n(), .groups = "drop") %>%
          arrange(Position)
      }, error = function(e) {
        stop(paste("Failed to process file:", vcf_file, "-", e$message))
      })
    } else {
      stop("Multiple or no files found for single-VCF mode.")
    }
  }
  
  
  final_data <- final_data %>%
    rename_with(
      ~ gsub("denovo_(id\\d+).vcf\\.gz", "\\1", .x), 
      .cols = 2:ncol(final_data)
    )
  
  # Return the final data
  return(final_data)
}


################# READ IN FASTA FILE


read_fasta <- function(fasta_file) {
  # Read the FASTA file
  sequences <- readDNAStringSet(fasta_file)
  
  # Convert to a named vector
  seq_data <- setNames(as.character(sequences), names(sequences))
  
  return(seq_data)
}


################# CALCULATE GC CONTENT


calculate_gc_content <- function(seq_data, window_size = 150) {
  gc_data <- list()
  
  for (seq_name in names(seq_data)) {
    sequence <- seq_data[[seq_name]]
    seq_length <- nchar(sequence)
    
    # Generate windows
    positions <- seq(1, seq_length, by = window_size)
    gc_values <- sapply(positions, function(pos) {
      window_seq <- substr(sequence, pos, min(pos + window_size - 1, seq_length))
      gc_count <- sum(str_count(window_seq, "[GCgc]"))
      return(gc_count / nchar(window_seq))  # GC content fraction
    })
    
    # Store as a dataframe
    gc_df <- data.frame(Position = positions, GC_Content = gc_values, Sequence = seq_name)
    gc_data[[seq_name]] <- gc_df
  }
  
  return(bind_rows(gc_data))  # Combine into one dataframe
}


################# PLOT CIRCOS

# 1. Function to add annotation tracks and circos axis
plot_annotations <- function(annotation_data) {
  
  # Define colors for annotation types
  dark2_colors <- brewer.pal(8, "Dark2")
  feature_colors <- c(
    "gene" = dark2_colors[8],
    "rRNA" = dark2_colors[6],
    "misc_feature" = dark2_colors[1],
    "repeat_region" = dark2_colors[2]
  )
  
  circos.track(
    ylim = c(0, 1), 
    bg.border = NA,
    track.height = 0.1, 
    panel.fun = function(x, y) {
      
      for (i in seq_len(nrow(annotation_data))) {
        start <- annotation_data$Minimum[i]
        end <- annotation_data$Maximum[i]
        strand <- annotation_data$Direction[i]
        feature <- annotation_data$Type[i]
        feature_name <- annotation_data$Name[i]
        
        arrow_position <- ifelse(strand == "reverse", "start", "end")
        
        # Draw arrow
        circos.arrow(
          x1 = start, x2 = end,
          arrow.head.width = 0.5,
          arrow.head.length = cm_x(0.1),
          col = adjustcolor(feature_colors[feature], alpha.f = 0.8),
          arrow.position = arrow_position,
          border = "gray20",
          lwd = 0.5
        )
        
        # Add text inside the arrow
        circos.text(
          x = (start + end) / 2,  # Center text in arrow
          y = 0,  # Place text inside the arrow
          labels = feature_name,
          facing = "bending.inside",  # Align text inside the arrow
          niceFacing = TRUE,
          cex = 1,  # Adjust text size
          font = 2,  # Make text bold
          adj = c(0.5, 0.5)  # Center text
        )
      }
    }
  )
}

# 2. Function to add depth and SNPs
plot_depth_and_snps <- function(depth_data, vcf_data = NULL, total_length, adapt_max_depth = FALSE) {
  
  if (ncol(depth_data) < 3) {
    depth_colors <- "blue"
  } else {
    depth_colors <- brewer.pal(9, "Blues")
    depth_colors <- tail(depth_colors, ncol(depth_data) - 2)
  }
  
  if (ncol(vcf_data) < 3) {
    vcf_colors <- "#DE2D26"
  } else {
    vcf_colors <- brewer.pal(ncol(vcf_data) - 1, "Reds")
  } 
  
  if (adapt_max_depth == TRUE) {
    max_depth <- median(as.matrix(depth_data[, 3])) + sd(as.matrix(depth_data[, 3]))
  } else {
    max_depth <- max(as.matrix(depth_data[, -c(1:2)]))
  }
  
  max_depth <- median(as.matrix(depth_data[, 9])) + median(as.matrix(depth_data[, 9]))
  
  circos.trackPlotRegion(
    factors = "genome",
    ylim = c(min(as.matrix(depth_data[, -c(1:2)])), max_depth),
    bg.border = NA,
    bg.col = 'gray90',
    track.height = 0.3,
    panel.fun = function(x, y) {
      # Add grid lines
      grid_intervals <- pretty(c(min(as.matrix(depth_data[, -c(1:2)])), max_depth), n = 5)
      for (grid_value in grid_intervals) {
        circos.segments(
          x0 = 0, x1 = total_length,
          y0 = grid_value, y1 = grid_value,
          col = "lightgray", lwd = 0.5
        )
      }
      
      # Plot depth data
      for (i in seq_along(depth_colors)) {
        circos.lines(
          x = depth_data$Position,
          y = depth_data[[i + 2]],
          col = depth_colors[i],
          lwd = 1
        )
      }
      
      # Find the closest window for each SNP position
      vcf_data$matched_depth_pos <- depth_data$Position[findInterval(vcf_data$Position, depth_data$Position)]
      
      # Plot SNPs if vcf_data is provided
      if (!is.null(vcf_data)) {
        for (vcf_col in names(vcf_data)[-1]) {  # Skip Position column
          matching_col <- names(depth_data)[names(depth_data) == vcf_col]  # Find matching column name
          
          if (length(matching_col) > 0) {  # If a match is found
            circos.points(
              x = vcf_data$matched_depth_pos,  # Use nearest 10-bp window
              y = depth_data[[matching_col]][match(vcf_data$matched_depth_pos, depth_data$Position)],  
              pch = 21,
              col = "black",  # Border color
              bg = vcf_colors[which(names(vcf_data) == vcf_col) - 1],  # Match colors correctly
              cex = 1.5
            )
          }
        }
      }
    }
  )
  
  # Add y-axis
  circos.yaxis(
    side = "left",
    at = pretty(c(min(as.matrix(depth_data[, -c(1:2)])), max_depth), n = 5),
    labels.cex = 0.6
  )
}

# 3. Function to add metadata text with conditional formatting
plot_metadata <- function(origin, total_length, depth_data) {
  
  mean_depth <- mean(depth_data)
  median_depth <- median(depth_data)
  sd_depth <- sd(depth_data)
  
  text(
    x = 0, y = 0.1,
    labels = bquote(bold(.(origin))),
    cex = 1.5,
    col = "black"
  )
  
  # Plot the total length without bold formatting
  text(
    x = 0, y = 0,
    labels = paste0(total_length, " bp"),
    cex = 1,
    col = "black"
  )
  
  # Plot mean, median, and standard deviation
  text(
    x = 0, y = -0.05,
    labels = paste0("Mean: ", round(mean_depth,2)),
    cex = 0.8,
    col = "black"
  )
  text(
    x = 0, y = -0.1,
    labels = paste0("Median: ", round(median_depth,2)),
    cex = 0.8,
    col = "black"
  )
  text(
    x = 0, y = -0.15,
    labels = paste0("SD: ", round(sd_depth,2)),
    cex = 0.8,
    col = "black"
  )
}

# 4. Function to add GC content
plot_gc_content <- function(gc_data, total_length) {
  circos.trackPlotRegion(
    factors = "genome",
    ylim = c(min(gc_data$GC_Content), max(gc_data$GC_Content, na.rm = TRUE)),  # Adjusted GC content range
    bg.border = NA,
    track.height = 0.3,
    panel.fun = function(x, y) {
      
      # Add grid lines
      # grid_intervals <- pretty(c(min(gc_data$GC_Content), max(gc_data$GC_Content)), n = 5)
      # for (grid_value in grid_intervals) {
      #   circos.segments(
      #     x0 = 0, x1 = total_length,
      #     y0 = grid_value, y1 = grid_value,
      #     col = "lightgray", lwd = 0.5
      #   )
      # }
      
      # Plot GC content
      circos.lines(
        x = gc_data$Position,
        y = gc_data$GC_Content,
        col = "darkgreen",
        lwd = 0.5
      )
    }
  )
  
  # Add y-axis labels
  circos.yaxis(side = "left", at = seq(0, 1, by = 0.2), labels.cex = 0.6)
  
}


######################### EXECUTE SCRIPT #########################



pdf("circosPlots.pdf", width = 12, height = 17)
par(mfrow = c(3, 2), mar = c(2, 2, 2, 2))



################## Ave River

depth_data <- read_depth_percent_id("AveRiv", ID = TRUE)

# Create a window column by dividing Position by 10 and taking the floor
depth_data_w10 <- depth_data %>%
  mutate(Window = floor(Position / 10)) %>%  # Group positions into windows
  group_by(Window) %>%
  summarise(across(where(is.numeric), ~median(.x, na.rm = TRUE)), .groups = "drop") 

total_length <- nrow(depth_data)

vcf_data <- read_vcf_data("AveRiv", ID = TRUE)

annotations_data <- read.csv("../5_Annotations/AveRiv_plastome_features.csv")[,1:7]

circos.clear()
circos.par(start.degree = 90)
circos.initialize(factors = "genome", xlim = c(0, total_length))

plot_annotations(annotations_data)

# Add genome axis
circos.track(ylim = c(-0.1, 0), bg.border = NA, track.height = 0.05, panel.fun = function(x, y) {
  circos.axis(
    h = "bottom",
    major.at = seq(0, total_length, by = 5000),
    labels = seq(0, total_length, by = 5000),
    labels.cex = 0.7,
    col = "black",
    lwd = 0.6,
    labels.col = "black"
  )
})

plot_depth_and_snps(depth_data_w10, na.omit(vcf_data[,c(1,9)]), total_length, adapt_max_depth = TRUE)

# Define the positions for the vertical lines
line_positions <- c(16800, 49500, 73500, 
                    90150, 125000, 133700)

# Loop through positions and add dashed pink lines
for (pos in line_positions) {
  circos.segments(
    x0 = pos, x1 = pos,  
    y0 = 0, y1 = 250,   # Adjust Y range to match depth track height
    sector.index = "genome",  
    col = "#FF69B4",        
    lty = 2,  # Dashed line
    lwd = 2   # Line thickness
  )
}

# Define the positions for the vertical lines
line_positions <- c(16800, 73500, 125000)

# Loop through positions and add dashed pink lines
for (pos in line_positions) {
  circos.segments(
    x0 = pos, x1 = pos,  
    y0 = 0, y1 = 250,   # Adjust Y range to match depth track height
    sector.index = "genome",  
    col = "#FF69B4",        
    lty = 1,  # full line
    lwd = 2   # Line thickness
  )
}



fasta_sequences <- read_fasta("AveRiv/AveRiv_plastome.fasta")
gc_data <- calculate_gc_content(fasta_sequences)
plot_gc_content(gc_data, total_length = max(gc_data$Position))
plot_metadata("AveRiv", total_length, depth_data$id99)


################## LakeMen


depth_data <- read_depth_percent_id("LakeMen", ID = TRUE)

# Create a window column by dividing Position by 10 and taking the floor
depth_data_w10 <- depth_data %>%
  mutate(Window = floor(Position / 10)) %>%  # Group positions into windows
  group_by(Window) %>%
  summarise(across(where(is.numeric), ~median(.x, na.rm = TRUE)), .groups = "drop") 

total_length <- nrow(depth_data)

annotations_data <- read.csv("../5_Annotations/LakeMen_plastome_features.csv")[,1:7]


vcf_data <- read_vcf_data("LakeMen", ID = TRUE)

circos.clear()
circos.par(start.degree = 90)
circos.initialize(factors = "genome", xlim = c(0, total_length))

plot_annotations(annotations_data)


# Add genome axis
circos.track(ylim = c(-0.1, 0), bg.border = NA, track.height = 0.05, panel.fun = function(x, y) {
  circos.axis(
    h = "bottom",
    major.at = seq(0, total_length, by = 5000),
    labels = seq(0, total_length, by = 5000),
    labels.cex = 0.7,
    col = "black",
    lwd = 0.6,
    labels.col = "black"
  )
})

plot_depth_and_snps(depth_data_w10, vcf_data = na.omit(vcf_data[,c(1,9)]), total_length, adapt_max_depth = FALSE)

fasta_sequences <- read_fasta("LakeMen/LakeMen_plastome.fasta")
gc_data <- calculate_gc_content(fasta_sequences)
plot_gc_content(gc_data, total_length = max(gc_data$Position))
plot_metadata("LakeMen", total_length, depth_data$id99)


################## ResCRep


depth_data <- read_depth_percent_id("ResCRep", ID = TRUE)
# Create a window column by dividing Position by 10 and taking the floor
depth_data_w10 <- depth_data %>%
  mutate(Window = floor(Position / 10)) %>%  # Group positions into windows
  group_by(Window) %>%
  summarise(across(where(is.numeric), ~median(.x, na.rm = TRUE)), .groups = "drop") 

total_length <- nrow(depth_data)

vcf_data <- read_vcf_data("ResCRep", ID = TRUE)

annotations_data <- read.csv("../5_Annotations/ResCRep_plastome_features.csv")[,1:7]


circos.clear()
circos.par(start.degree = 90)
circos.initialize(factors = "genome", xlim = c(0, total_length))

plot_annotations(annotations_data)

# Add genome axis
circos.track(ylim = c(-0.1, 0), bg.border = NA, track.height = 0.05, panel.fun = function(x, y) {
  circos.axis(
    h = "bottom",
    major.at = seq(0, total_length, by = 5000),
    labels = seq(0, total_length, by = 5000),
    labels.cex = 0.7,
    col = "black",
    lwd = 0.6,
    labels.col = "black"
  )
})

plot_depth_and_snps(depth_data_w10, vcf_data = na.omit(vcf_data[,c(1,9)]), total_length, adapt_max_depth = TRUE)

fasta_sequences <- read_fasta("ResCRep/ResCRep_plastome.fasta")
gc_data <- calculate_gc_content(fasta_sequences)
plot_gc_content(gc_data, total_length = max(gc_data$Position))
plot_metadata("ResCRep", total_length, depth_data$id99)







dev.off()

