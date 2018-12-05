# fantasy_scoring
Evaluate how scoring type and roster settings affect fantasy positional value.

Explore an [interactive demo](http://rynecarbone.github.io/ff/scoring/) showcasing the final result.

## Installation
```bash
git clone https://github.com/rynecarbone/fantasy_scoring
cd fantasy_scoring
pip install -r requirements.txt
```

## Quick Run
All of the fantasy data is available in `data/espn_fantasy_data_small.csv`. The bokeh dashboard can be
generated with :
```bash
python scripts/create_bokeh_dash.py
```


## Retreiving data from ESPN
An R script is available to scrape the raw data from ESPN. It includes examples of how to create 
the same plots as the dashboard for quick data exploration. The script is located in `R/scrape_espn_ff_raw.R`.
You should download RStudio, and install the packages: `tidyverse`, `Rcurl`, and `XML`.
