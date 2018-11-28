#!/usr/bin/env python

"""Create a bokeh dashboard with customJS callbacks for interactivity"""

import pandas as pd
from bokeh.layouts import row, column, layout, gridplot, widgetbox 
from bokeh.models import (CustomJS, RangeSlider, HoverTool, Range1d,
                          CheckboxGroup, CheckboxButtonGroup, Select, RadioButtonGroup,
                          ColumnDataSource, CDSView, BooleanFilter, GroupFilter, LinearAxis, LogAxis)
from bokeh.plotting import figure, output_file,show
from bokeh.palettes import Category10  
from bokeh.models.widgets import Tabs, Panel
from bokeh.models.tickers import FixedTicker

__author__ = 'Ryne Carbone'

# Define some global variables
default_roster = dict(nQB=1, nRB=2, nWR=2, nTE=1, nFLEX=2, nTEAMS=10)
POSITIONS = ['QB','RB','TE','WR']
PPR_TYPES = ['STD','HPPR','PPR']
PFD_TYPES = ['STD','HPFD','PFD']
sorts = ['PPR_type','PFD_type','Pos','RankPos']
groups = ['PPR_type','PFD_type','Pos','RankPos']
# Read input data
df = (pd.read_csv('data/espn_fantasy_data_small.csv')
        .sort_values(by=sorts)
        .reset_index(drop=True))
YEARS = sorted(list(df.Year.unique()))
# Split data for determining Flex info: [year_ind][pfd_ind][ppr_ind]
df = df.sort_values(by=['Year','PPR_type','PFD_type','Pts'],ascending=False).reset_index(drop=True)
# pscript behaves weirdly with Float64Array: https://github.com/flexxui/pscript/issues/5
df['Pts']=df['Pts'].astype(str)
l_df_flex = [[[ColumnDataSource(df[(df.Year==y)&(df.PFD_type==pf)&(df.PPR_type==pp)][['Pos','RankPos','Pts']]) 
               for pp in ['STD','HPPR','PPR']] 
              for pf in ['STD','HPFD','PFD']] 
             for y in sorted(list(set(df.Year.tolist())))]

# Define source data
raw_source = ColumnDataSource({'Pos': [], 'RankPos': [], 'PPR_type': [], 'PFD_type': [], 'AVG': []})
raw_source_bands = ColumnDataSource({'Pos': [], 'RankPos': [], 'PPR_type': [], 'PFD_type': [], 'BANDS': []})
rel_source = ColumnDataSource({'Pos': [], 'RankPos': [], 'PPR_type': [], 'PFD_type': [], 'REL': []})
rel_source_bands = ColumnDataSource({'Pos': [], 'RankPos': [], 'PPR_type': [], 'PFD_type': [], 'BANDS': []})
rb_source = ColumnDataSource({'Pos': [], 'RankPos': [], 'PPR_type': [], 'PFD_type': [], 'RBB': []})
rb_source_bands = ColumnDataSource({'Pos': [], 'RankPos': [], 'PPR_type': [], 'PFD_type': [], 'BANDS': []})

# Define the tools
tools = 'box_zoom,wheel_zoom,pan,reset,save'
rel_hover = HoverTool(tooltips=[ 
  ('Pos. Rank','@Pos-@RankPos'),
  ('Rel. Val.', '@REL{0.000 a}')])
rb_hover = HoverTool(tooltips=[ 
  ('Pos. Rank','@Pos-@RankPos'),
  ('Rel. Val. (RB baseline)', '@RBB{0.000 a}')])
hover = HoverTool(tooltips=[
  ('Pos. Rank','@Pos-@RankPos'),
  ('Avg Pts.', '@AVG{0.0 a}')])

# Create filters for using views 
qb = GroupFilter(column_name='Pos', group='QB')
rb = GroupFilter(column_name='Pos', group='RB')
te = GroupFilter(column_name='Pos', group='TE')
wr = GroupFilter(column_name='Pos', group='WR')
ppr = GroupFilter(column_name='PPR_type', group='PPR')
hppr = GroupFilter(column_name='PPR_type', group='HPPR')
sppr = GroupFilter(column_name='PPR_type', group='STD')
pfd = GroupFilter(column_name='PFD_type', group='PFD')
hpfd = GroupFilter(column_name='PFD_type', group='HPFD')
spfd = GroupFilter(column_name='PFD_type', group='STD')
ppr_filters = [sppr, hppr, ppr]
pfd_filters = [spfd, hpfd, pfd]
pos_filters=  [qb, rb, te, wr]


def callback(l_df_flex=l_df_flex, default_roster=None, yr_lo=None, yr_hi=None, init_run=False,
              rel_source=rel_source, rel_source_bands=rel_source_bands,
              raw_source=raw_source, raw_source_bands=raw_source_bands,
              rb_source=rb_source, rb_source_bands=rb_source_bands,
              YEARS=YEARS, PPR_TYPES=PPR_TYPES, PFD_TYPES=PFD_TYPES, window=None):
  """Massive callback written with pscript, converted to javascript
  
  Note: this is a terrible way to calculate averages/min/max by groups,
  but pscript doesn't allow the use of python packages, and bokeh only allows
  python through pscript with standalone html pages :(
  :param l_df_flex: list of data frames in 3D array (by year, ppr type, pfd type)
  :param default_roster: roster settings for initial plot creation
  :param yr_lo: min year for initial plot creation
  :param yr_hi: max year for initial plot creation
  :param init_run: flag to indicate if initial plot creation
  :param rel_source: ColumnDataSource holding data for relative value plot
  :param rel_source_bands: ColumnDataSource holidng data for relative value plot bands
  :param raw_source: ColumnDataSource holding data for raw value plot
  :param raw_source_bands: ColumnDataSource holidng data for raw value plot bands
  :param rb_source: ColumnDataSource holding data for RB baseline plot
  :param rb_source_bands: ColumnDataSource holidng data for RB baseline plot bands
  :param YEARS: list of possible years to select
  :param PPR_TYPES: list of possible ppr setings
  :param PFD_TYPES: list of possible pfd types
  :param window: allows access to javascript functions/variables
  :return: Updated ColumnDataSources on initial run, otherwise bokeh handles JS updates
  """
  # Read in the roster settings for calculating flex replacement values
  if default_roster:
    roster = default_roster
  else:
    roster = dict(nQB=nQB.value, nRB=nRB.value, nWR=nWR.value, 
                nTE=nTE.value, nFLEX=nFLEX.value, nTEAMS=nTEAMS.value)
  # Set the correct year range for selecting data
  if not yr_lo:
    yr_lo = year.value[0]
  if not yr_hi:
    yr_hi = year.value[1]
  yr_range = range(YEARS.index(yr_lo), YEARS.index(yr_hi)+1)
  # Connect to plotting data sources
  data = rel_source.data
  data_bands = rel_source_bands.data
  data_raw = raw_source.data
  data_bands_raw = raw_source_bands.data
  data_rb = rb_source.data
  data_bands_rb = rb_source_bands.data
  bmax=[]; bmin=[]; bmax_raw=[]; bmin_raw=[]; bmax_rb=[]; bmin_rb=[]
  # Reset return data
  return_keys = ['Pos','RankPos','PPR_type','PFD_type','REL']
  return_keys_raw = ['Pos','RankPos','PPR_type','PFD_type','AVG']
  return_keys_rb = ['Pos','RankPos','PPR_type','PFD_type','RBB']
  return_bands_keys = ['Pos','RankPos','PPR_type','PFD_type','BANDS']
  for rk,rkr,rkrb,rbk in zip(return_keys, return_keys_raw, return_keys_rb, return_bands_keys):
    data[rk]=[]; data_bands[rbk]=[];
    data_raw[rkr]=[]; data_bands_raw[rbk]=[]
    data_rb[rkrb]=[]; data_bands_rb[rbk]=[]
  # Start calculating 
  for i_pfd in range(3):
    for i_ppr in range(3):
      raw_val_dict = {}
      rel_val_dict = {}
      for i_y in yr_range:
        # Get data for this year and scoring type
        i_df = l_df_flex[i_y][i_pfd][i_ppr].data
        # For each year, calculate replacement value by position
        rep_val = {'QB':0,'RB':0,'WR':0,'TE':0}
        nflex_left = int(roster['nFLEX'])*int(roster['nTEAMS'])
        # Update values until filled up all flex spots 
        for pos, rk, pts in zip(i_df['Pos'],i_df['RankPos'],i_df['Pts']):
          if pos=='QB' and rk==roster['nQB']*roster['nTEAMS']:
            rep_val['QB'] = pts 
          elif pos!='QB' and int(rk) > int(roster[f'n{pos}'])*int(roster['nTEAMS']) and nflex_left>0:
            rep_val[pos] = pts
            nflex_left -= 1
          elif pos!='QB' and int(rk)==int(roster[f'n{pos}'])*int(roster['nTEAMS']):
            rep_val[pos] = pts
        # Calculate pts over rep
        pts_over_rep = []
        for pos, pts in zip(i_df['Pos'],i_df['Pts']):
          tmp_pts = (float(pts) - float(rep_val[pos])) if (float(pts)-float(rep_val[pos])) > 0 else 0
          pts_over_rep.append(tmp_pts)
        # Convert to relative val
        tot_rep_pts = sum(pts_over_rep)
        i_rel_val = [i_por/tot_rep_pts for i_por in pts_over_rep]
        # Make dict for rel val/raw val by pos and rank
        for pos, rk, pts, rv in zip(i_df['Pos'], i_df['RankPos'], i_df['Pts'], i_rel_val):
          key = f'{pos}_{rk}'
          # Create entry for raw val 
          if not raw_val_dict.get(key):
            raw_val_dict[key] = []
          raw_val_dict[key].append(float(pts))
          # Only create entry for rel val if val not 0
          if rv == 0: continue
          if not rel_val_dict.get(key):
            rel_val_dict[key] = []
          rel_val_dict[key].append(100*rv)
      # Get the average per year, and bands for rel val
      l_pos = []; l_rankpos = []; l_rel = []; l_bmin = []; l_bmax = []
      l_pos_rb = []; l_rankpos_rb = []; l_rel_rb = []; l_bmin_rb = []; l_bmax_rb = []
      for kk, vv in rel_val_dict.items():
        k_pos, k_rk = kk.split('_')
        l_pos.append(k_pos)
        l_rankpos.append(int(k_rk))
        l_rel.append(sum(vv)/len(vv))
        l_bmin.append(min(vv))
        l_bmax.append(max(vv))
        if rel_val_dict.get(f'RB_{k_rk}'):
          vv_rb = rel_val_dict.get(f'RB_{k_rk}')
          l_pos_rb.append(k_pos)
          l_rankpos_rb.append(int(k_rk))
          l_rel_rb.append((sum(vv)/len(vv))/(sum(vv_rb)/len(vv_rb)))
          l_bmin_rb.append(min(vv)/(sum(vv_rb)/len(vv_rb)))
          l_bmax_rb.append(max(vv)/(sum(vv_rb)/len(vv_rb)))
      # Get the average per year, and bands or raw val
      l_pos_raw = []; l_rankpos_raw = []; l_raw = []; l_bmin_raw = []; l_bmax_raw = []
      for kk, vv in raw_val_dict.items():
        k_pos, k_rk = kk.split('_')
        l_pos_raw.append(k_pos)
        l_rankpos_raw.append(int(k_rk))
        l_raw.append(sum(vv)/len(vv))
        l_bmin_raw.append(min(vv))
        l_bmax_raw.append(max(vv))
      # Sort by rankpos?
      sorted_zip = sorted(zip(l_pos, l_rankpos, l_rel, l_bmin, l_bmax), 
                          key=lambda tup: (tup[0], float(tup[1])/1000.))
      sorted_zip_raw = sorted(zip(l_pos_raw, l_rankpos_raw, l_raw, l_bmin_raw, l_bmax_raw), 
                              key=lambda tup: (tup[0], float(tup[1])/1000.))
      sorted_zip_rb = sorted(zip(l_pos_rb, l_rankpos_rb, l_rel_rb, l_bmin_rb, l_bmax_rb), 
                          key=lambda tup: (tup[0], float(tup[1])/1000.))
      for tup in sorted_zip:
        data['Pos'].append(tup[0])
        data['RankPos'].append(tup[1])
        data['PFD_type'].append(PFD_TYPES[i_pfd])
        data['PPR_type'].append(PPR_TYPES[i_ppr])
        data['REL'].append(float(tup[2]))
        bmin.append(float(tup[3]))
        bmax.append(float(tup[4]))
      for tup in sorted_zip_raw:
        data_raw['Pos'].append(tup[0])
        data_raw['RankPos'].append(tup[1])
        data_raw['PFD_type'].append(PFD_TYPES[i_pfd])
        data_raw['PPR_type'].append(PPR_TYPES[i_ppr])
        data_raw['AVG'].append(float(tup[2]))
        bmin_raw.append(float(tup[3]))
        bmax_raw.append(float(tup[4]))
      for tup in sorted_zip_rb:
        data_rb['Pos'].append(tup[0])
        data_rb['RankPos'].append(tup[1])
        data_rb['PFD_type'].append(PFD_TYPES[i_pfd])
        data_rb['PPR_type'].append(PPR_TYPES[i_ppr])
        data_rb['RBB'].append(float(tup[2]))
        bmin_rb.append(float(tup[3]))
        bmax_rb.append(float(tup[4]))
  # Create bands for relative data
  data_bands['Pos'] = list(data['Pos']) + list(reversed(data['Pos']))
  data_bands['RankPos'] = list(data['RankPos']) + list(reversed(data['RankPos']))
  data_bands['PFD_type'] = list(data['PFD_type']) + list(reversed(data['PFD_type']))
  data_bands['PPR_type'] = list(data['PPR_type']) + list(reversed(data['PPR_type']))
  data_bands['BANDS'] = list(bmin) + list(reversed(bmax))
  # Create bands for raw data
  data_bands_raw['Pos'] = list(data_raw['Pos']) + list(reversed(data_raw['Pos']))
  data_bands_raw['RankPos'] = list(data_raw['RankPos']) + list(reversed(data_raw['RankPos']))
  data_bands_raw['PFD_type'] = list(data_raw['PFD_type']) + list(reversed(data_raw['PFD_type']))
  data_bands_raw['PPR_type'] = list(data_raw['PPR_type']) + list(reversed(data_raw['PPR_type']))
  data_bands_raw['BANDS'] = list(bmin_raw) + list(reversed(bmax_raw))
  # Create bands for rb baseline data
  data_bands_rb['Pos'] = list(data_rb['Pos']) + list(reversed(data_rb['Pos']))
  data_bands_rb['RankPos'] = list(data_rb['RankPos']) + list(reversed(data_rb['RankPos']))
  data_bands_rb['PFD_type'] = list(data_rb['PFD_type']) + list(reversed(data_rb['PFD_type']))
  data_bands_rb['PPR_type'] = list(data_rb['PPR_type']) + list(reversed(data_rb['PPR_type']))
  data_bands_rb['BANDS'] = list(bmin_rb) + list(reversed(bmax_rb))
  if not init_run:
    rel_source.change.emit()
    rel_source_bands.change.emit()
    raw_source.change.emit()
    raw_source_bands.change.emit()
    rb_source.change.emit()
    rb_source_bands.change.emit()
  else:
    return raw_source, raw_source_bands,rel_source, rel_source_bands, rb_source, rb_source_bands

# Run initial relative value code
raw_source, raw_source_bands, rel_source, rel_source_bands, rb_source, rb_source_bands  = callback(default_roster=default_roster, yr_lo=2002, yr_hi=2017, init_run=True)

# Define relative value y-axis ticks
rel_ticker =  FixedTicker(ticks=[0.01,0.02,0.05,0.1, 0.2, 0.5, 1,2,5,10,20,50,100], 
                          minor_ticks=[0.03,0.04,0.06,0.07,0.08,0.09,0.3,0.4,0.6,0.7,0.8,0.9,3,4,6,7,8,9,30,40,60,70,80,90])
raw_ticker = FixedTicker(ticks=[1,100,200,300,400, 500,600,700,800, 900,1000],
                         minor_ticks=[20,40,60,80,120,140,160,180,220,240,260,280,320,340,360,380,420,440,460,480,520,
                                      540,560,580,620,640,660,680,720,740,760,780,820,840,860,880,920,940,960,980])

# Keep list of figures and glyphs
plots=[]; rel_plots=[]; rb_plots=[]
lines_list=[]; rel_lines_list=[]; rb_lines_list=[]
patches_list=[]; rel_patches_list=[]; rb_patches_list=[]
# Create Raw Pts plots separately for each score type 
for i_pfd, (FD, pfd_type) in enumerate(zip(pfd_filters, PFD_TYPES)):
  for i_ppr, (PR, ppr_type) in enumerate(zip(ppr_filters, PPR_TYPES)):
    # Keep track of glyphs in each figure
    lines=[]; rel_lines=[];rb_lines=[]
    patches=[]; rel_patches=[]; rb_patches=[]
    # Create figure 
    temp_plot = figure(tools=[hover, tools], x_axis_label='Position Rank', x_range=Range1d(start=-1, end=122, bounds="auto"),
                       y_axis_label='Points', y_axis_type='log', y_range=Range1d(start=20,end=520, bounds=(20,1000)))
    temp_plot2 = figure(tools=[rel_hover, tools], x_axis_label='Position Rank',
                      y_axis_label='Relative Value (%)', y_axis_type='log', y_range=Range1d(start=0.01, end=11, bounds=(0.01,100)))
    temp_plot3 = figure(tools=[rb_hover, tools], x_axis_label='Position Rank', 
                      y_axis_label='Relative Value (RB==1)', y_axis_type='log', y_range=Range1d(start=0.04, end=25, bounds=(0.01,100)))
    temp_plot.yaxis.ticker = raw_ticker
    temp_plot2.yaxis.ticker = rel_ticker 
    temp_plot3.yaxis.ticker = rel_ticker 
    # Add line and patch for each position
    for i, (P, p, c) in enumerate(zip(pos_filters, POSITIONS, reversed(Category10[4]))):
      # Only creae legend in upper right plot
      leg = p if (i_pfd==0 and i_ppr==2) else None
      # Create a line for the avg pts by positional rank 
      l1 = temp_plot.line('RankPos', 'AVG', source=raw_source, legend=leg, color=c, line_width=3,
                          view=CDSView(source=raw_source, filters=[P, FD, PR]))
      l2 = temp_plot2.line('RankPos','REL', source=rel_source, legend=leg, color=c, line_width=3,
              view=CDSView(source=rel_source, filters=[P, FD, PR]))
      l3 = temp_plot3.line('RankPos','RBB', source=rb_source, legend=leg, color=c, line_width=3,
              view=CDSView(source=rb_source, filters=[P, FD, PR]))
      # Create a patched area for range of values 
      p1 = temp_plot.patch('RankPos', 'BANDS', source=raw_source_bands, legend=leg, color=c, alpha=0.1, 
                            view=CDSView(source=raw_source_bands, filters=[P, FD, PR]))
      p2 = temp_plot2.patch('RankPos','BANDS', source=rel_source_bands, legend=leg, color=c, alpha=0.1,
               view=CDSView(source=rel_source_bands, filters=[P, FD, PR]))
      p3 = temp_plot3.patch('RankPos','BANDS', source=rb_source_bands, legend=leg, color=c, alpha=0.1,
               view=CDSView(source=rb_source_bands, filters=[P, FD, PR]))
      # Update lines/patches lists 
      lines.append(l1)
      rel_lines.append(l2)
      rb_lines.append(l3)
      patches.append(p1)
      rel_patches.append(p2)
      rb_patches.append(p3)
    # Add lists to master list of lists for all figures
    lines_list.append(lines)
    rel_lines_list.append(rel_lines)
    rb_lines_list.append(rb_lines)
    patches_list.append(patches)
    rel_patches_list.append(rel_patches)
    rb_patches_list.append(rb_patches)
    # Hide yaxis
    if i_ppr != 0:
      temp_plot.yaxis.visible=False
      temp_plot2.yaxis.visible=False
      temp_plot3.yaxis.visible=False
    # Hide xaxis
    if i_pfd != 2:
      temp_plot.xaxis.visible=False
      temp_plot2.xaxis.visible=False
      temp_plot3.xaxis.visible=False
    # Add Row labels on right
    if i_ppr == 2:
      temp_plot.add_layout(LogAxis(axis_label=pfd_type,axis_label_text_font_style='bold', ticker=raw_ticker), 'right')
      temp_plot2.add_layout(LogAxis(axis_label=pfd_type,axis_label_text_font_style='bold', ticker=rel_ticker), 'right')
      temp_plot3.add_layout(LogAxis(axis_label=pfd_type,axis_label_text_font_style='bold', ticker=rel_ticker), 'right')
    # Add Col labels on top
    if i_pfd == 0:
      temp_plot.title.text = ppr_type
      temp_plot.title.align = 'center'
      temp_plot2.title.text = ppr_type
      temp_plot2.title.align = 'center'
      temp_plot3.title.text = ppr_type
      temp_plot3.title.align = 'center'
    plots.append(temp_plot)
    rel_plots.append(temp_plot2)
    rb_plots.append(temp_plot3)
# Synchronize all axis ranges for linked panning/zooming
for p, p2, p3 in zip(plots, rel_plots, rb_plots):
  p.x_range = plots[8].x_range  
  p.y_range = plots[8].y_range
  p2.x_range = rel_plots[8].x_range
  p2.y_range = rel_plots[8].y_range
  p3.x_range = rb_plots[8].x_range
  p3.y_range = rb_plots[8].y_range
# Arrange into a grid 
grid = gridplot(plots, ncols=3, toolbar_location='right', 
                plot_height=250, plot_width=300)
grid_rel = gridplot(rel_plots, ncols=3, toolbar_location='right', 
                plot_height=250, plot_width=300)
grid_rb = gridplot(rb_plots, ncols=3, toolbar_location='right',
                  plot_height=250, plot_width=300)
# Define the layout
row_plot = row(children=[grid]) 

# This changes the alpha of the lines and fills
def checkbox_callback(checkboxes=None, lines_list=None, patches_list=None, 
                      rel_lines_list=None, rel_patches_list=None,
                      rb_lines_list=None, rb_patches_list=None):
  for l in lines_list:
    l[0].glyph.line_alpha = 1 if 0 in checkboxes.active else 0.15 
    l[1].glyph.line_alpha = 1 if 1 in checkboxes.active else 0.15
    l[2].glyph.line_alpha = 1 if 2 in checkboxes.active else 0.15
    l[3].glyph.line_alpha = 1 if 3 in checkboxes.active else 0.15 
  for l in rel_lines_list:
    l[0].glyph.line_alpha = 1 if 0 in checkboxes.active else 0.15 
    l[1].glyph.line_alpha = 1 if 1 in checkboxes.active else 0.15
    l[2].glyph.line_alpha = 1 if 2 in checkboxes.active else 0.15
    l[3].glyph.line_alpha = 1 if 3 in checkboxes.active else 0.15 
  for l in rb_lines_list:
    l[0].glyph.line_alpha = 1 if 0 in checkboxes.active else 0.15 
    l[1].glyph.line_alpha = 1 if 1 in checkboxes.active else 0.15
    l[2].glyph.line_alpha = 1 if 2 in checkboxes.active else 0.15
    l[3].glyph.line_alpha = 1 if 3 in checkboxes.active else 0.15 
  for p in patches_list:
    p[0].glyph.fill_alpha = .1 if 0 in checkboxes.active else 0.01 
    p[1].glyph.fill_alpha = .1 if 1 in checkboxes.active else 0.01 
    p[2].glyph.fill_alpha = .1 if 2 in checkboxes.active else 0.01 
    p[3].glyph.fill_alpha = .1 if 3 in checkboxes.active else 0.01 
  for p in rel_patches_list:
    p[0].glyph.fill_alpha = .1 if 0 in checkboxes.active else 0.01 
    p[1].glyph.fill_alpha = .1 if 1 in checkboxes.active else 0.01 
    p[2].glyph.fill_alpha = .1 if 2 in checkboxes.active else 0.01 
    p[3].glyph.fill_alpha = .1 if 3 in checkboxes.active else 0.01 
  for p in rb_patches_list:
    p[0].glyph.fill_alpha = .1 if 0 in checkboxes.active else 0.01 
    p[1].glyph.fill_alpha = .1 if 1 in checkboxes.active else 0.01 
    p[2].glyph.fill_alpha = .1 if 2 in checkboxes.active else 0.01 
    p[3].glyph.fill_alpha = .1 if 3 in checkboxes.active else 0.01 

# This selects y-axis
def pts_callback(grid_raw=grid, grid_rel=grid_rel, grid_rb=grid_rb):
  if points.value == 'Raw Points':
    row_plot.children=[grid_raw]
  elif points.value == 'Relative Value':
    row_plot.children=[grid_rel]
  elif points.value == 'RB Baseline':
    row_plot.children = [grid_rb]

# Slider (select year range)
rel_callback = CustomJS.from_py_func(callback)
year_slider = RangeSlider(start=2002, end=2017, step=1, value=(2002, 2017), title='Year Range', 
                          callback=rel_callback, callback_policy='mouseup')
rel_callback.args['year'] = year_slider

# Add checkbox for toggling lines for positions
checkboxes = CheckboxButtonGroup(labels=POSITIONS, active=[0,1,2,3],  
                           callback=CustomJS.from_py_func(checkbox_callback))
checkboxes.callback.args = dict(checkboxes=checkboxes, lines_list=lines_list, patches_list=patches_list,
                                rel_lines_list=rel_lines_list, rel_patches_list=rel_patches_list,
                                rb_lines_list=rb_lines_list, rb_patches_list=rb_patches_list)

# Add select for defining roster
select_qb = Select(options=["1","2"], value="1", title='nQB', callback=rel_callback)
rel_callback.args['nQB'] = select_qb
select_rb = Select(options=["1","2","3","4"], value="2", title='nRB', callback=rel_callback)
rel_callback.args['nRB'] = select_rb
select_wr = Select(options=["1","2","3","4"], value="2", title='nWR',callback=rel_callback)
rel_callback.args['nWR'] = select_wr
select_te = Select(options=["1","2"], value="1", title='nTE',callback=rel_callback)
rel_callback.args['nTE'] = select_te
select_flex = Select(options=["1","2","3","4"], value="2", title='nFlex',callback=rel_callback)
rel_callback.args['nFLEX'] = select_flex
select_teams = Select(options=["8","10","12","14"], value="10", title='nTeams',callback=rel_callback)
rel_callback.args['nTEAMS'] = select_teams

# Add select for choosing Y-axis
pts_callback = CustomJS.from_py_func(pts_callback)
pts_group = Select(title='Y-axis', options=['Raw Points', 'Relative Value','RB Baseline'], value='Raw Points', callback=pts_callback)
pts_callback.args['points'] = pts_group 
pts_callback.args['row_plot'] = row_plot

# Add widgetbox, define page layout
wbox = widgetbox(children=[pts_group, year_slider, checkboxes, select_qb, select_rb, 
                           select_wr, select_te, select_flex, select_teams], width=200)
l = row(children=[wbox, row_plot], sizing_mode='scale_height')

# Define output location
output_file('output/dashboard.html')

show(l)

