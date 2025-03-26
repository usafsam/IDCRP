#Shiny app for USNA WW

library(shiny)
library(bslib)
library(plotly)
library(readxl)
library(tidyverse)

#######Load Data
#Upload twist output from WW_Twist_list_from_kraken_10Jan25.rmd
#Add _USNA to control names (of text) before running Twist list script, or only run on samples of interest and remove "filter(grepl("USNA", SampeID))" to ensure controls are included

#Rank includes S, S1, and S2 from twist_list - (check for double counting))

Twist_output <- read_xlsx("/Users/vhogan/Documents/WW/Data/WW_twist_output_29Jan25.xlsx")

#Select USNA samples and fields of interest, simplify SampleID to date
WW_twist <- Twist_output %>%
  filter(grepl("USNA", SampleID)) %>%
  select(SampleID, Name, NCBI_ID, Num_fragments_covered) %>% #select only needed columns
  mutate(NCBI_ID = as.character(NCBI_ID)) %>% #read as character to avoid any calculations
  mutate(SampleID = ifelse(startsWith(SampleID, 'USNA'), str_split_i(SampleID,'_', i =3), str_split_i(SampleID, '_', 1)))

#Upload list of Twist targets in RPP from ARIA
ARIA_viruses <- read_xlsx("/Users/vhogan/Documents/WW/Data/ARIA_resp_viruses.xlsx") %>%
  mutate(NCBI_ID = as.character(NCBI_ID))

#Upload list of GI viruses from Frontier paper (Singh et al 2024)
GI_viruses <- read_xlsx("/Users/vhogan/Documents/WW/Data/GI_viruses.xlsx") %>%
  mutate(NCBI_ID = as.character(NCBI_ID))

#Upload weekly_counts converted week manually in excel to week_starting column (e.g. week 14, 2024 is week starting April 1 2024)
weekly_counts <- read_xlsx("/Users/vhogan/Documents/WW/Data/Weekly_counts_Apr_Nov24.xlsx") 

######Transform data
#pivot and sum reads by column
reads_by_name <- WW_twist %>%
  pivot_wider(
    names_from = SampleID,
    values_from = Num_fragments_covered) 

#create single row df with total reads
totals <- reads_by_name %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
  mutate_all(as.numeric) 

#divide reads df by totals - should direcly compare to Agans' excel workbook
percent <- reads_by_name %>%
  rowwise() %>% 
  mutate(c_across(3:ncol(reads_by_name)) / totals, .keep = "unused") %>% #calculate percent rowwise for each column starting at column 3
  ungroup()

#long version for plotting bar chart
percent_long <- percent %>%
  pivot_longer(
    cols = 3:ncol(percent),
    names_to = 'SampleID',
    values_to = 'Percent'
  ) %>%
  mutate(Percent = Percent *100) %>%
  filter(str_detect(SampleID, "^2")) %>% #remove controls for plotting 
  mutate(date = as.Date(SampleID, "%y%b%d"))

#Filter for viruses in ARIA
ARIA_resp_reads <- percent_long %>%
  filter(NCBI_ID %in% ARIA_viruses$NCBI_ID) 

#Rearrange to match ARIA_resp_reads df
weekly_counts_new <-weekly_counts %>%
  select(week_starting, n_encounters) %>%
  rename(date = week_starting) %>%
  mutate(Name = 'Weekly Count',
         NCBI_ID = NA,
         SampleID = NA, 
         Percent = NA)

#Bind rows with weekly counts to ARIA_resp_reads to plot together
ARIA_resp <- ARIA_resp_reads %>%
  mutate(n_episodes = NA) %>%
  bind_rows(weekly_counts_new)


###########Shiny
ui <- page_navbar(
  theme = bs_theme(version = 4),
  title = "USNA Wastewater",
  fillable = FALSE,
  
  nav_panel(
    title = 'All',
    h5("Major species detected in wastewater that are present in Twist Viral Pathogen Panel"),
    "To calculate percent abundance, number of reads was divided by total number of Twist 
                pathogen reads by sample/date. This chart includes all species with >5% abundance in each sample.
                To adjust the percent cutoff, enter a percent below. E.g. 100 will show all species detected 
                that are present in Twist Viral Pathogen Panel",
    numericInput("obs_percent_all", "Choose Percent:", 5, min = 0, max = 100),
    plotlyOutput("graph_all")
  ),
  
  nav_panel(
    title = 'Respiratory',
    h5("Species detected in wastewater that are present in Twist Viral Pathogen Panel 
                 and tested for in Respiratory Pathogen Panel (RPP) used in ARIA"),
    "To calculate percent abundance, number of reads was divided by total number of Twist 
                pathogen reads by sample/date. This chart includes all species with >0.01% abundance in each sample.
                To adjust the percent cutoff, enter a percent below. E.g. 100 will show all species detected 
                that are present in Twist Viral Pathogen Panel",
    numericInput("obs_percent_resp", "Choose Percent:", 0.01, min = 0, max = 100),
    plotlyOutput("graph_resp")
  ),
  
  nav_panel(
    title = 'GI',
    h5("Species detected in wastewater that are present in Twist Viral Pathogen Panel and common GI viruses of concern"),
    "To calculate percent abundance, number of reads was divided by total number of Twist 
                pathogen reads by sample/date. This chart includes all species with >0.01% abundance in each sample.
                To adjust the percent cutoff, enter a percent below. E.g. 100 will show all species detected 
                that are present in Twist Viral Pathogen Panel",
    numericInput("obs_percent_gi", "Choose Percent:", 0.01, min = 0, max = 100),
    plotlyOutput("graph_gi")
  )
    
)

server <- function(input, output, session){
  
  
  
  #filter for chosen percent (all viruses)
  percent_all <- reactive({percent_long %>%
    filter(Percent > input$obs_percent_all)})
  

  
  output$graph_all <- renderPlotly({
    #Stacked bar chart with Percent > 5
    plot_ly(percent_all(), x = ~date, y = ~Percent, color = ~Name, type = 'bar',
            hovertemplate = paste(
              "%{yaxis.title.text}: %{y:.0}<br>",
              "%{xaxis.title.text}: %{x:%W}<br>")) %>% 
      layout(xaxis = list(title = 'Date'),
             yaxis = list(title = 'Percent Abundance'),
             legend = list(orientation = 'h', y = -0.4), 
             margin = list(r = 50),
             barmode = 'stack')
  })
  
  
  #filter for chosen percent (respiratory viruses)
  percent_resp <- reactive({ARIA_resp %>%
      filter(Percent > input$obs_percent_resp | Name =='Weekly Count')})
  
  
  
  output$graph_resp <- renderPlotly({
    #Stacked bar chart with Percent > 0.01
    plot_ly(percent_resp(), x = ~date, y = ~Percent, color = ~Name, type = 'bar') %>% 
      add_trace(x = ~date, y = ~n_encounters, type = 'scatter', mode = 'lines', yaxis = "y2", marker = list(
        color = 'rgb(17, 157, 255)')) %>%
      
      layout(xaxis = list(title = 'Date'),
             yaxis = list(title = 'Percent Abundance'),
             yaxis2 = list(title = "Weekly Encounters",
                           overlaying = "y", side = "right", 
                           showgrid = FALSE, range = list(0,250)),
             legend = list(orientation = 'h', y = -0.4), 
             margin = list(r = 50),
             barmode = 'stack')
  })
  
  #filter for chosen percent (GI viruses)
  percent_gi <- reactive({ARIA_resp %>%
      filter(Percent > input$obs_percent_gi)})
  
  
  
  output$graph_gi <- renderPlotly({
    #Stacked bar chart with Percent > 0.01
    plot_ly(percent_resp(), x = ~date, y = ~Percent, color = ~Name, type = 'bar') %>%
      
      layout(xaxis = list(title = 'Date'),
             yaxis = list(title = 'Percent Abundance'),
             legend = list(orientation = 'h', y = -0.4), 
             margin = list(r = 50),
             barmode = 'stack')
  })  
}

shinyApp(ui, server)