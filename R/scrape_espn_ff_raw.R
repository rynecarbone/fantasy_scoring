library(XML)
library(RCurl)
library(tidyverse)

##
## DATA SCRAPING
##
# Retreive data from ESPN
scrapeDataESPN <- function(urlprefix, urlend, startyr, endyr, tot_positions, data_type){
  master <- data.frame()
  # Loop through specified years
  for (j in startyr:endyr){
    cat('Loading Year', j, '\n')
    # Loop through pagination (in steps of 40 on ESPN)
    for (i in seq(1, tot_positions, 40)){
      #cat(' > Loading player first downs from rank: ', i, 'to ', (i+20), '\n')
      # Grab the table
      URL <- paste(urlprefix, as.character(j), urlend, as.character(i), sep="")
      link_data <- getURL(URL)
      table <- readHTMLTable(link_data, stringsAsFactors=F)[[1]]
      # If table isn't empty, add to master table
      if (!is.null(table)){
        # Scoring table is weird, fix column names
        if(data_type=='Scoring'){
          # Rename header/remove first row
          new_names <- as.character(unlist(table[1,]))
          colnames(table) <- new_names
          table <- table[-1,] 
        }
        # Add column for year
        table$Year <- j
        master <- rbind(table, master)
      }
    }
  }
  # Fix the table data, convert to tibble
  return_data <- fix_data(master, data_type)
  return(return_data)
}

# Fix column names, and data
fix_data <- function(data, data_type){
  # Scoring data is weird
  if(data_type=='Scoring'){ data <- data[which(data$PLAYER != 'TOUCHDOWNS'),]}
  # Remove repeated header in rows
  data <- data[which(data$RK != 'RK'),]
  # Separate Player/Pos column, rename multi-team player team names, change to numeric
  r <- data %>% as.tibble() %>%
    separate(PLAYER, into = c('Player','Pos'), sep=', ') %>%
    mutate(n_slashes=str_count(TEAM, '/'), 
           Team=if_else(n_slashes > 0, paste0((n_slashes+1),'TM'), TEAM))%>%
    select(-RK, -TEAM) %>%
    mutate_at(vars(-Player, -Team, -Pos, -Year),
              funs(as.numeric(gsub(",","",.))))
  # Rename and select columns to return
  if (data_type=='Passing'){
    r <- r %>% 
      select(Player, Pos, Team, Passing.Att=ATT, Passing.Comp=COMP, 
             Passing.Yds=YDS, Passing.TD=TD, Passing.Int=INT, Year)
  } else if (data_type=='Receiving'){
    r <- r %>% select(Player, Pos, Team, Receiving.Tar=TAR, Receiving.Rec=REC, Receiving.Yds=YDS, 
                      Receiving.TD=TD, Receiving.FD=`1DN`, Receiving.Fum=FUM, Year)
  } else if (data_type=='Rushing'){
    r <- r %>% select(Player, Pos, Team, Rushing.Att=ATT, Rushing.Yds=YDS, 
                      Rushing.TD=TD, Rushing.FD=`1DN`, Rushing.Fum=FUM, Year)
  } else if (data_type=='Scoring'){
    r <- r %>% select(Player, Pos, Team, Ret.TD=RET, `2PT`, Year)
  }
  # Sometimes more than one row per player if on multiple teams?
  r <- r %>%  
    group_by(Player, Pos, Team, Year) %>% 
    summarise_all(sum) %>% ungroup()
  return(r)
}

##
## CALCULATE FANTASY POINTS/RANKINGS
##
# Calculate raw fantasy points
calc_points <- function(data, PPR, FD, score){
  data <- data %>% 
    mutate(!!paste0(score):=(
      Passing.Yds/25. +
        Passing.TD*4 +
        Passing.Int*(-2) +
        Rushing.Yds/10. +
        Rushing.TD*6 +
        Rushing.Fum*(-2) +
        Rushing.FD*FD +
        Receiving.Rec*PPR +
        Receiving.Yds/10. +
        Receiving.TD*6 +
        Receiving.Fum*(-2) +
        Receiving.FD*FD +
        Ret.TD*6 +
        `2PT`*(2) 
    ))
  return(data %>% as.tibble())
}

# Calculate position rankings by league scoring
calc_ranks <- function(data, scoring){
  pts_var <- sym(paste0('Pts.', scoring))
  rank_pos_var <- sym(paste0('RankPos.', scoring))
  data <- data %>%  arrange(Year, Pos, desc((!!pts_var)))%>%
    group_by(Pos, Year) %>%
    mutate( (!!rank_pos_var) := rank(-(!!pts_var), ties="first"))
  return(data)
}

# Calculate Value over Replacement
calc_rep_val <- function(data, n_qb, n_rb, n_wr, n_te, n_flex, n_teams){
  # Define how many starters by position 
  n_starting <- c('QB'=n_teams*n_qb,
                  'WR'=n_teams*n_wr,
                  'RB'=n_teams*n_rb,
                  'TE'=n_teams*n_te,
                  'FLEX'=n_teams*n_flex)
  roster <-bind_rows(c('nQB'=n_qb, 'nWR'=n_wr, 'nRB'=n_rb, 
                       'nTE'=n_te, 'nFLEX'=n_flex, 'nTEAMS'=n_teams))
  # Find flex values by scoretype
  n_pos_flex <- data %>% 
    filter((Pos=='RB' & RankPos > n_starting['RB'][[1]]) |
             (Pos=='WR' & RankPos > n_starting['WR'][[1]]) |
             (Pos=='TE' & RankPos > n_starting['TE'][[1]])) %>%  
    arrange(Year, Pos, PPR_type, PFD_type, RankPos, Pts) %>%
    group_by(Year, PPR_type, PFD_type )%>% 
    mutate(flex_rank = rank(-Pts, ties="first")) %>% 
    group_by(PPR_type, PFD_type, Pos, Year) %>% 
    summarise(n = length(flex_rank[flex_rank<=n_starting['FLEX'][[1]]])) %>% 
    group_by(PPR_type, PFD_type, Pos) %>% 
    summarise(nFlex=round(mean(n)))
  # Calculate relative value % 
  data <- data %>% 
    left_join(n_pos_flex, by=c('PPR_type','PFD_type','Pos')) %>%
    group_by(Pos, Year, PPR_type, PFD_type)%>%
    mutate(n_start = n_starting[Pos],
           pos_rep = if_else(!is.na(nFlex), n_start[1]+nFlex, n_start[1]),
           rep_pts = nth(Pts, pos_rep[1]),
           pts_ovr_rep = if_else(Pts-rep_pts>0, Pts-rep_pts, 0))%>% 
    ungroup() %>%
    group_by(Year, PPR_type, PFD_type) %>%
    mutate(rel_val = pts_ovr_rep/sum(pts_ovr_rep)) %>%
    select(-nFlex, -n_start,-pos_rep, -rep_pts,-pts_ovr_rep) 
  # Add roster rules to data
  data <- crossing(data, roster)
  return(data)
}


##
## RUN CODE
##

# ESPN urls
espn_url_prefix_pass <- "http://www.espn.com/nfl/statistics/player/_/stat/passing/sort/passingYards/year/"
espn_url_prefix_rush <- "http://www.espn.com/nfl/statistics/player/_/stat/rushing/sort/rushingYards/year/"
espn_url_prefix_rec <- "http://www.espn.com/nfl/statistics/player/_/stat/receiving/sort/receivingYards/year/"
espn_url_prefix_score <- "http://www.espn.com/nfl/statistics/player/_/stat/scoring/sort/totalPoints/year/"
espn_url_suffix <- "/qualified/false/count/"

# Date params
year_start = 2002
year_end = 2017

# Get data
espn_fantasy_pass <- scrapeDataESPN(espn_url_prefix_pass, espn_url_suffix, year_start, year_end, 300, 'Passing')
espn_fantasy_rush <- scrapeDataESPN(espn_url_prefix_rush, espn_url_suffix, year_start, year_end, 300, 'Rushing')
espn_fantasy_rec <- scrapeDataESPN(espn_url_prefix_rec, espn_url_suffix, year_start, year_end, 300, 'Receiving')
espn_fantasy_score <- scrapeDataESPN(espn_url_prefix_score, espn_url_suffix, year_start, year_end, 500, 'Scoring')

# Join all the data, replace NA values with 0
espn_fantasy <- plyr::join_all(list(espn_fantasy_pass, espn_fantasy_rec, espn_fantasy_rush, espn_fantasy_score),
                               by=c('Player','Pos','Team','Year'),
                               type='full') %>% 
  replace(., is.na(.), 0)

# Add fantasy points
espn_fantasy <- espn_fantasy %>%
  calc_points(PPR=0, FD=0, score='Pts.STD_STD') %>%
  calc_points(PPR=0, FD=0.5, score='Pts.STD_HPFD') %>%
  calc_points(PPR=0, FD=1, score='Pts.STD_PFD') %>%
  calc_points(PPR=0.5, FD=0, score='Pts.HPPR_STD') %>%
  calc_points(PPR=0.5, FD=0.5, score='Pts.HPPR_HPFD') %>%
  calc_points(PPR=0.5, FD=1, score='Pts.HPPR_PFD') %>%
  calc_points(PPR=1, FD=0, score='Pts.PPR_STD') %>%
  calc_points(PPR=1, FD=0.5, score='Pts.PPR_HPFD') %>%
  calc_points(PPR=1, FD=1, score='Pts.PPR_PFD')

# Add positional ranks
espn_fantasy <- espn_fantasy %>%
  calc_ranks('STD_STD') %>%
  calc_ranks('STD_HPFD') %>%
  calc_ranks('STD_PFD') %>%
  calc_ranks('HPPR_STD') %>%
  calc_ranks('HPPR_HPFD') %>%
  calc_ranks('HPPR_PFD') %>%
  calc_ranks('PPR_STD') %>%
  calc_ranks('PPR_HPFD') %>%
  calc_ranks('PPR_PFD') 

# Only keep QB, RB, WR, TE
espn_fantasy <- espn_fantasy %>% filter(Pos %in% c('QB', 'RB', 'WR', 'TE'))

# Tidied
espn_tidy <- espn_fantasy %>% 
  gather(matches('Pts|Rank'), key='Scoring', value='Value') %>% 
  separate(Scoring, into=c('pts_rk', 'PPR_type','PFD_type')) %>% 
  spread(key='pts_rk', value='Value')
# Order scoring type factors
espn_tidy$PPR_type_f = factor(espn_tidy$PPR_type, levels=c('STD','HPPR','PPR'))
espn_tidy$PFD_type_f = factor(espn_tidy$PFD_type, levels=c('STD','HPFD','PFD'))

##
## SAVE TO CSV
##
# At this point you can save the espn_tidy data frame
espn_tidy %>% to_csv('espn_tidy_full.csv')


##
## PLOTTING
##

MIN_YEAR <- 2012
MIN_PTS <- 50
# Plot all years together -- Avg is line, range is 'error bars'
# Plotted on log-y: most positions have log drop off after "top" tier
espn_tidy %>% 
  filter(Year> MIN_YEAR) %>%
  group_by(Pos, RankPos, PPR_type_f, PFD_type_f)%>% 
  summarise(avg_pts_pos=mean(Pts),
            min_pts_pos=min(Pts),
            max_pts_pos=max(Pts))%>%
  filter( max_pts_pos > MIN_PTS)%>%
  ggplot(aes(x=RankPos, group=Pos)) + 
  geom_ribbon(aes(ymin=if_else(min_pts_pos>MIN_PTS, min_pts_pos, MIN_PTS), ymax=max_pts_pos, fill=Pos), alpha=0.3)+
  geom_line(aes(y=avg_pts_pos, color=Pos))+
  facet_grid(PPR_type_f~PFD_type_f)+
  labs(x='Position Rank', y='Fantasy Points')+
  theme_bw() + scale_y_continuous(trans='log1p', limits=c(MIN_PTS,500))


# Calculate relative value
# Need to specify the roster settings
espn_tidy_rel <- calc_rep_val(espn_tidy, n_qb=1, n_rb=2, n_wr=2, 
                              n_te=1, n_flex=2, n_teams=10) 
# Plot relative value
espn_tidy_rel %>%   
  filter(rel_val > 0, Year > MIN_YEAR)%>%
  group_by(Pos, RankPos, PPR_type_f, PFD_type_f)%>%
  summarise(avg_val_pos=mean(rel_val), 
            min_val_pos=min(rel_val), 
            max_val_pos=max(rel_val))%>%
  ggplot(aes(x=RankPos)) + 
  geom_line(aes(y=100*avg_val_pos, group=Pos, color=Pos))+
  geom_ribbon(aes(ymin=100*min_val_pos, 
                  ymax=100*max_val_pos, group=Pos, fill=Pos), alpha=0.2)+
  facet_grid(PPR_type_f~PFD_type_f)+
  labs(x='Position Rank', y='Relative Value (%)')+
  theme_bw() 

# Rel value with RB as baseline
espn_tidy_rel %>%   
  filter(rel_val > 0, Year > MIN_YEAR)%>%
  group_by(Year, RankPos, PPR_type_f, PFD_type_f)%>%
  mutate(rel_val_baseline = if_else(length(rel_val[Pos=='RB'])>0, rel_val[Pos=='RB'][1], 0),
         rel_val = if_else(rel_val_baseline>0, rel_val / rel_val_baseline, 0)) %>%
  filter(rel_val_baseline > 0.001) %>%
  group_by(Pos, RankPos, PPR_type_f, PFD_type_f)%>%
  summarise(avg_val_pos=mean(rel_val), 
            min_val_pos=min(rel_val), 
            max_val_pos=max(rel_val))%>%
  ggplot(aes(x=RankPos)) + 
  geom_line(aes(y=avg_val_pos, group=Pos, color=Pos))+
  geom_ribbon(aes(ymin=min_val_pos, 
                  ymax=max_val_pos, group=Pos, fill=Pos), alpha=0.2)+
  scale_y_continuous(trans='log10', breaks=c(0.1,.5, 1,5), minor_breaks=c(.2, .3, .4, .6, .7,.8,.9,2,3,4))+
  facet_grid(PPR_type_f~PFD_type_f, scales="free")+
  labs(x='Position Rank', y='Relative Value (RB baseline)')+
  theme_bw() 
