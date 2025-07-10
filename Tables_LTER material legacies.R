#### Tables for LTER material legacy synthesis project

# Load required packages
library(tidyverse)
library(gt)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

#### Table S1: Data descriptions and sources ----
library(tidyverse)
library(flextable)
library(officer)
library(pagedown)

# Read the data
data_sources <- read_csv("Figures/Tables/Table of data sources_LTER material legacies.csv")

# Create the flextable
data_sources.table <- flextable(data_sources) %>%
  compose(
    j = "URL",
    i = 1:nrow(data_sources),
    value = as_paragraph(
      as_chunk("Link", props = fp_text(color = "blue", underline = TRUE),
               hyperlink = data_sources$URL)
    )
  ) %>%
  merge_v(j = c("Ecosystem", "LTER (or other) site")) %>%
  valign(j = c("Ecosystem", "LTER (or other) site"), valign = "top") %>%
  set_header_labels(
    `LTER (or other) site` = "Site",
    `Description of dataset` = "Dataset Description",
    `Timespan of dataset` = "Timespan",
    URL = "Data DOI/URL"
  ) %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  width(j = "Timespan of dataset", width = 1.2) %>%  # ⬅️ Wider Timespan column
  bold(part = "header") %>%
  align(align = "center", part = "header") %>%
  autofit() %>%
  border_remove() %>%
  border_outer() %>%
  border_inner_h()

data_sources.table

# Export as PDF
flextable::save_as_pdf(data_sources.table, path = "data_sources_table.pdf")


#### Table S2: GLMM structures and outputs ----
library(tidyverse)
library(gt)

# Read in the data
model_table <- read_csv("Figures/GLMM structures.csv") %>% 
  gt() %>%
  tab_header(
    title = md("**Model Structures by Ecosystem (Raw Relationships)**")
  ) %>%
  cols_label(
    ecosystem_type = "Ecosystem Type",
    material_legacy_predictor = "Material Legacy Predictor (Fixed Effect)",
    demographic_response = "Demographic Response",
    'distribution_(link)' = "Distribution (Link Function)",
    random_effects = "Random Effects"
  ) %>%
  tab_options(
    table.font.names = "Times New Roman",
    table.width = pct(100),
    column_labels.font.weight = "bold",
    table.font.size = 12,
    heading.title.font.size = 14,
    data_row.padding = px(5)
  ) 

