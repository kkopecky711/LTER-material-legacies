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

#### Table S2: GLMM structures and outputs ----
library(tidyverse)
library(gt)

# Order ecosystem types according to main figures
ecosystem_order <- c(
  "Temperate forest",
  "Mangrove forest",
  "Oyster reef",
  "Temperate rainforest",
  "Boreal forest",
  "Kelp forest",
  "Saltmarsh",
  "Coral reef",
  "Tropical forest",
  "Tallgrass prairie"
)

# Read in the data
model_table <- read_csv("Figures/Tables/GLMM structures.csv") %>% 
  mutate(ecosystem_type = factor(ecosystem_type, levels = ecosystem_order)) %>%
  arrange(ecosystem_type) %>% 
  gt() %>%
  tab_header(
    title = md("**Model Structures by Ecosystem (Raw Relationships)**")
  ) %>%
  cols_label(
    ecosystem_type = "Ecosystem Type",
    material_legacy_predictor = "Material Legacy Predictor (Fixed Effect)",
    demographic_response = "Demographic Response",
    'distribution_(link function)' = "Distribution (Link Function)",
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

model_table

#### Table S3: GLMM model outputs ----
library(tidyverse)
library(flextable)
library(officer)

# Define the folder and matching file names with ecosystem labels
file_paths <- list.files("Datasets/Effect sizes/Raw/", pattern = "*.csv", full.names = TRUE)
ecosystems <- c("Boreal forest", "Mangrove forest", "Saltmarsh", "Temperate forest", "Temperate rainforest", "Tallgrass prairie", "Tropical forest",
                  "Coral reef", "Kelp forest", "Oyster reef")

# Read and combine
model_outputs <- map2_dfr(file_paths, ecosystems, ~ {
  read_csv(.x) %>% mutate(Ecosystem = .y)
})

# Combine all model outputs into one table
model_outputs <- map2_dfr(file_paths, ecosystems, ~ {
  read_csv(.x) %>% mutate(Ecosystem = .y)
})

# Clean and format model output
model_outputs_clean <- model_outputs %>%
  filter(effect == "fixed") %>%
  select(Ecosystem, estimate, conf.low, conf.high, p.value) %>%
  mutate(
    across(
      .cols = c(estimate, conf.low, conf.high),
      .fns = ~ ifelse(abs(.) < 0.001,
                      formatC(., format = "e", digits = 2),
                      formatC(., format = "f", digits = 3))
    ),
    p.value = case_when(
      is.na(p.value) ~ NA_character_,
      p.value < 0.001 ~ "< 0.001",
      TRUE ~ formatC(p.value, format = "f", digits = 3)
    )
  )


# Order ecosystem types according to main figures
ecosystem_order <- c(
  "Temperate forest",
  "Mangrove forest",
  "Oyster reef",
  "Temperate rainforest",
  "Boreal forest",
  "Kelp forest",
  "Saltmarsh",
  "Coral reef",
  "Tropical forest",
  "Tallgrass prairie"
)

model_outputs_clean <- model_outputs_clean %>%
  mutate(Ecosystem = factor(Ecosystem, levels = ecosystem_order)) %>%
  arrange(Ecosystem)

# Create flextable
flextable_outputs <- flextable(model_outputs_clean) %>%
  set_header_labels(
    estimate = "Estimate",
    conf.low = "95% CI (Low)",
    conf.high = "95% CI (High)",
    p.value = "p-value"
  ) %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  bold(part = "header") %>%
  align(part = "header", align = "center") %>%
  autofit() %>%
  border_remove() %>%
  border_outer() %>%
  border_inner_h()

# View the table
flextable_outputs

# Optional: Save to Word
save_as_docx(flextable_outputs, path = "model_outputs_table.docx")

# Optional: Save to PDF (requires Chrome + pagedown)
# install.packages("pagedown")
# pagedown::chrome_print("model_outputs_table.html")


#### Table S4: Outputs from standardized models ----
library(tidyverse)
library(flextable)

# List and load standardized output files
std_paths <- list.files("Datasets/Effect sizes/Standardized/", full.names = TRUE, pattern = "\\.csv$")

# Update this list to match the order of standardized model files:
ecosystems <- c(
  "Boreal forest", "Mangrove forest", "Saltmarsh", "Temperate forest", "Temperate rainforest", "Tallgrass prairie", "Tropical forest",
  "Coral reef", "Kelp forest", "Oyster reef"
)

# Read and combine files with ecosystem labels
std_outputs <- map2_dfr(std_paths, ecosystems, ~ {
  read_csv(.x) |> mutate(Ecosystem = .y)
})

# Filter fixed effects and select relevant columns
std_outputs_clean <- std_outputs %>%
  filter(term == "bs_stg_ba" | term == "standardized_predictor" | term == "dead_holdfast_cover" | effect == "fixed") %>%
  select(Ecosystem, estimate, conf.low, conf.high, p.value)

# Format numbers (use scientific notation for small values)
std_outputs_clean <- std_outputs_clean %>%
  mutate(
    across(
      .cols = c(estimate, conf.low, conf.high),
      .fns = ~ ifelse(abs(.) < 0.001,
                      formatC(., format = "e", digits = 2),
                      formatC(., format = "f", digits = 3))
    ),
    p.value = case_when(
      is.na(p.value) ~ NA_character_,
      p.value < 0.001 ~ "< 0.001",
      TRUE ~ formatC(p.value, format = "f", digits = 3)
    )
  )

# Reorder ecosystems to match the forest plot
ecosystem_order <- c(
  "Temperate forest", "Mangrove forest", "Oyster reef", "Temperate rainforest",
  "Boreal forest", "Kelp forest", "Saltmarsh", "Coral reef",
  "Tropical forest", "Tallgrass prairie"
)

std_outputs_clean <- std_outputs_clean %>%
  mutate(Ecosystem = factor(Ecosystem, levels = ecosystem_order)) %>%
  arrange(Ecosystem)

# Create the flextable
flextable_std <- flextable(std_outputs_clean) %>%
  set_header_labels(
    estimate = "Standardized estimate",
    conf.low = "95% CI (Low)",
    conf.high = "95% CI (High)",
    p.value = "p-value"
  ) %>%
  fontsize(size = 11, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  bold(part = "header") %>%
  align(part = "header", align = "center") %>%
  autofit() %>%
  border_remove() %>%
  border_outer() %>%
  border_inner_h()

flextable_std
  