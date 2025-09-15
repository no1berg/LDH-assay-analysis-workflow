# LDH Cytotoxicity Assay Analysis Shiny App

![Language](https://img.shields.io/badge/Language-R-blue.svg)
![Framework](https://img.shields.io/badge/Framework-Shiny-brightgreen.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

A flexible and robust R Shiny application for the automated analysis and visualization of 96-well plate LDH cytotoxicity assay data. This tool is designed to move from raw plate reader output to a full statistical report with minimal user intervention.



---

## Features

- **Interactive Visualization**: Plots are generated with `ggiraph` to allow hovering over data points for precise values.
- **Dynamic Control Assignment**: Use your own naming conventions for controls; simply assign them to their roles (e.g., Spontaneous Lysis, Blank) within the app.
- **Group-Specific Normalization**: Accurately normalizes data from complex experimental designs where each treatment group has its own dedicated maximum lysis control.
- **Robust Statistical Analysis**: Automatically checks statistical assumptions (normality and homogeneity of variances) and selects the appropriate test (ANOVA or Kruskal-Wallis), or allows for manual override.
- **Dynamic Interpretation**: The "Assumption Checks" and "Interpretation Guide" tabs provide clear, data-driven explanations of your results.
- **Full Reproducibility**:
    - **Downloadable Reports**: Generate a complete HTML report containing the plots, tables, and a summary of the analysis.
    - **Exportable Data**: Download the final, processed cytotoxicity data as a CSV file for use in other software.

---

## How to Use

### Running the App Locally

To run this application on your own machine, you will need R and RStudio installed.

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/no1berg/LDH-assay-analysis-workflow-assay-workflow.git](https://github.com/no1berg/LDH-assay-analysis-workflow)
    cd ldh-assay-analysis-workflow
    ```

2.  **Install Required Packages:**
    Open R or RStudio and run the following command to install all necessary packages:
    ```R
    install.packages(c(
      "shiny", "tidyverse", "car", "broom", "ggpubr",
      "dunn.test", "DT", "ggiraph", "rmarkdown"
    ))
    ```

3.  **Run the App:**
    Open the `app.R` file in RStudio and click the "Run App" button, or run the following command in the R console:
    ```R
    shiny::runApp("app.R")
    ```

---

## Input File Format

The application requires two `.csv` files with specific formatting.

### 1. Plate Layout File (`layout.csv`)

This file maps the contents of each well on the 96-well plate.

- The file **must not** have row names (e.g., A, B, C...).
- The first row is the column header (e.g., 1, 2, 3...).
- The cells should contain your custom group names.

**Important Naming Convention:**
For group-specific normalization, the maximum lysis control must be named using a consistent prefix followed by the exact name of the treatment group it corresponds to. The default prefix is `Tritonx100-`.

**Example `layout.csv`:**
```csv
1,2,3,4,5,6
empty,empty,empty,Vehicle,Vehicle,Vehicle
empty,empty,empty,DrugA,DrugA,DrugA
untreated,untreated,untreated,Tritonx100-Vehicle,Tritonx100-Vehicle,Tritonx100-Vehicle
untreated,untreated,untreated,Tritonx100-DrugA,Tritonx100-DrugA,Tritonx100-DrugA
```

### 2. Raw Absorbance File (`absorbance.csv`)

This file contains the raw optical density (OD) readings from the plate reader.

- The file **must not** have any headers or row names.
- Data must be in an **alternating row format**:
    - Row 1: 490nm readings for plate row A
    - Row 2: 680nm readings for plate row A
    - Row 3: 490nm readings for plate row B
    - Row 4: 680nm readings for plate row B
    - ...and so on.

---

## Methodology

1.  **Data Correction**: Raw absorbance is background-corrected by subtracting the 680nm reading from the 490nm reading for each well.
2.  **Normalization**: Percent cytotoxicity is calculated using the standard formula, with options for both plate-wide and group-specific controls:
![Cytotoxicity Formula](https://latex.codecogs.com/svg.latex?%5Ctext%7BCytotoxicity%20%28%5C%25%29%7D%20%3D%20%5Cfrac%7B%28%5Ctext%7BOD%7D_%7B%5Ctext%7BSample%7D%7D%20-%20%5Ctext%7BOD%7D_%7B%5Ctext%7BSpontaneous%7D%7D%29%7D%7B%28%5Ctext%7BOD%7D_%7B%5Ctext%7BMaximum%7D%7D%20-%20%5Ctext%7BOD%7D_%7B%5Ctext%7BSpontaneous%7D%7D%29%7D%20%5Ctimes%20100)
3.  **Statistical Testing**: The app uses a linear model to test for assumptions of normality (Shapiro-Wilk test) and homogeneity of variances (Levene's test). Based on these results, it performs either a One-Way ANOVA with a Tukey HSD post-hoc test or a Kruskal-Wallis test with a Dunn's post-hoc test.

---

## License

This project is licensed under the MIT License.
