# Interactive Heatmap and Dendrogram Explorer

An interactive R Shiny application for exploring hierarchical clustering patterns and associated data shown in heatmaps using various distance and clustering methods through zoomable plots.

## Background

This application combines hierarchical clustering dendrograms with interactive heatmaps to provide a comprehensive visualization of high-dimensional data patterns. Users can simultaneously explore both the clustering structure (via dendrogram) and the underlying data values (via heatmap) with synchronized zooming, cluster visualization, and multiple export options. The application is designed to handle any numerical dataset suitable for hierarchical clustering analysis.

The example data provided for application exploration was sourced from this publication: https://doi.org/10.1126/sciadv.abf5733. Lopes, et. al. 2021. Systematic dissection of transcriptional regulatory networks by genome-scale and single-cell CRISPR screens. Science Advances.

## Features

### Interactive Clustering
- **Multiple Distance Methods**: Choose from Manhattan, Euclidean, Maximum, or Binary distance calculations
- **Multiple Clustering Methods**: Select from Complete, Single, Average, or Ward.D2 linkage methods
- **Dynamic Cut Height**: Adjust clustering granularity with real-time visualization
- **Color-coded Clusters**: Each cluster gets a unique, high-contrast color applied to both dendrogram labels and gene names
- **Cluster Intersection Labels**: See cluster numbers displayed at cut height intersections

### Dual Visualization
- **Side-by-side Layout**: Heatmap and dendrogram displayed together with synchronized coordinates
- **Coordinated Views**: Both plots maintain alignment for easy pattern identification
- **Conditional Y-axis Labels**: Row numbers shown in full view, hidden when zoomed
- **Value-based Coloring**: Heatmap uses blue-white-red gradient for intuitive data interpretation

### Advanced Navigation & Zoom
- **Numeric Input Zoom**: Precise Y-axis zoom control using minimum and maximum row number inputs
- **Gene Search & Zoom**: Find and zoom to specific genes by name
- **Reset Functionality**: Easily return to full view with reset

### Comprehensive Export Options
- **Cluster Data Export**: Download cluster assignments as CSV files with descriptive filenames
- **High-Quality Image Export**: Export current view (full or zoomed) as publication-ready 300 DPI PNG files
- **Dynamic Filenames**: All exports automatically include parameters and zoom status for easy organization

### Visual Enhancement Features
- **Cut Height Visualization**: Red dashed line shows current clustering threshold
- **Missing Data Handling**: Graceful handling of NA values with heatmap cells shown in grey
- **Responsive Design**: Clean, professional layout optimized for data exploration
- **Real-time Updates**: All visualizations update instantly as parameters change

## Getting Started

### Prerequisites
```r
# Required R packages
install.packages(c("shiny", "dplyr", "tibble", "tidyr", "ggplot2", "ggdendro", "patchwork"))
```

### Running the Application
1. Clone this repository:
   ```bash
   git clone https://github.com/wendysphillips/interactive_dendrogram.git
   cd interactive_dendrogram
   ```

2. Use the provided example data OR prepare your data file:
   - Replace `"tre_data.tsv"` in the script with your actual data file path
   - Ensure data has items to be clustered as rows and conditions/variables as columns
   - First column should contain labels of items to be clustered

3. Launch the Shiny app:
   ```r
   # In R console
   source("shiny_heatmap_dendrogram.r")
   create_zoom_heatmap()
   ```
   The application will open in your default web browser.  
   
   OR 
   
   ```r
   # In terminal
   Rscript shiny_heatmap_dendrogram.r
   ```
   Then paste provided link into web browser.


## Project Structure

```
interactive_dendrogram/
├── shiny_heatmap_dendrogram.r  # Main heatmap Shiny application
├── generate_colors.R            # Color palette generator for clusters
├── tre_data.tsv          # Or replace with your dataset
├── README.md                  # This documentation
└── LICENSE                     # MIT License
```

## How to Use

### Basic Operation
1. **Select Methods**: Choose your preferred distance calculation and clustering method from the dropdowns
2. **Adjust Cut Height**: Use the numeric input to set clustering granularity and see cluster boundaries
3. **Explore Patterns**: Observe how clustering patterns relate to heatmap data values

### Zoom and Navigation
4. **Numeric Zoom**: Enter specific row numbers in Y-axis minimum and maximum fields, then click "Apply Y-axis Zoom"
5. **Gene Search**: Enter a gene name and click "Zoom to Gene" for automatic centering and zoom
6. **Reset View**: Click "Reset Zoom" to return to the full dataset view

### Data Export
7. **Export Clusters**: Click "Export Cluster Data" to download CSV with cluster assignments
8. **Export Images**: Click "Export Image (PNG)" to save current view as high-resolution image
9. **File Organization**: All exports include parameter settings and zoom status in filenames


### Using Your Own Data

This application works with any numerical dataset suitable for hierarchical clustering:

**Data Format Requirements:**
- Tab-separated values (.tsv)
- Row names in first column (genes, samples, conditions, etc.)
- Numerical data in subsequent columns
- Data should be pre-processed (scaled, normalized) as appropriate for your analysis

**Example Data Structure:**
```
Gene_ID    Condition_1    Condition_2    Condition_3    ...
Gene_A     2.5           -1.2           0.8            ...
Gene_B     -0.3          3.1            -2.1           ...
Gene_C     1.7           0.4            1.9            ...
...        ...           ...            ...            ...
```

**Setup Steps:**
1. Place your data file in the same directory as the script
2. Update the data loading line: `df <- read.table("your_data_filename.tsv", header = TRUE, row.names = 1)`
3. Ensure proper data preprocessing (scaling, normalization) before clustering

## Technical Details

- **Distance Methods**: Manhattan, Euclidean, Maximum, Binary
- **Clustering Methods**: Complete linkage, Single linkage, Average linkage, Ward.D2 (minimum variance)
- **Visualization Framework**: ggplot2 with ggdendro for dendrogram extraction and patchwork for plot combination
- **Coordinate System**: Custom inversion logic to align ggplot2 coordinates with intuitive row numbering
- **Reactive Programming**: Shiny reactive framework ensures efficient real-time updates
- **Export Formats**: CSV (cluster data) and PNG (images) with customizable parameters
- **Image Quality**: 300 DPI, 16"×10" dimensions optimized for publications

### Performance Considerations
- Application scales well with datasets up to several thousand rows
- Large datasets may experience slower rendering during zoom operations
- Memory usage scales with dataset dimensions and number of clusters

## Contributing

Contributions are welcome! Areas for improvement include:
- Additional distance/clustering methods
- Enhanced zoom controls (X-axis zoom)
- Alternative color schemes
- Performance optimizations for very large datasets
- Additional export formats

Please feel free to submit issues, feature requests, or pull requests.

## Authors

This application was developed jointly by **Wendy Phillips** and **GitHub Copilot** through collaborative AI-assisted programming.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Links

- [Example data source publication](https://doi.org/10.1126/sciadv.abf5733)
- [R Shiny Documentation](https://shiny.rstudio.com/)
- [ggplot2 Documentation](https://ggplot2.tidyverse.org/)
- [ggdendro Package](https://cran.r-project.org/package=ggdendro)
- [patchwork Package](https://patchwork.data-imaginist.com/)
- [Hierarchical Clustering in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust)