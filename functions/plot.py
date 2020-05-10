# Plot flattened graph and output HTML

from sklearn.manifold import MDS
import pandas as pd
from bokeh.models import ColumnDataSource, HoverTool, Legend
from bokeh.plotting import figure
from bokeh.embed import components

def plot_target(matrix, full_ids, boundaries, data, target, sim_threshold):

    # Flatten matrix
    flatten = MDS(n_components=2)
    flatten_matrix = flatten.fit_transform(matrix)

    # Make df for Bokeh
    bokeh_df = pd.DataFrame({'x': flatten_matrix[:, 0],
                             'y': flatten_matrix[:, 1],
                             'id': list(full_ids),
                             'name': pd.merge(full_ids, data, on='drugbank-id')['name'],
                             'desc': pd.merge(full_ids, data, on='drugbank-id')['description']
                            })

    # Plot
    p = figure(plot_width=1100, plot_height=800, \
        title=f"{bokeh_df[bokeh_df['id'] == target]['name'].item()}; Similarity threshold = {sim_threshold}")
    p.axis.visible = False
    p.grid.visible = False
    colors = ['black', 'green', 'red']
    labels = ['Target drug', 'Similar approved', 'Similar experimental']
    legend_items = []
    for i, (stop1, stop2) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        c = p.circle(x='x', y='y', color=colors[i], size=16, source=ColumnDataSource(bokeh_df[stop1:stop2]))
        legend_items.append((f'{labels[i]} ({stop2 - stop1})', [c]))

    p.add_tools(HoverTool(
        tooltips=[
            ('DrugBank ID', '@id'),
            ('Name', '@name'),
            ('Description', '@desc')
        ]))

    legend = Legend(items=legend_items)
    legend.click_policy='hide'

    p.add_layout(legend, 'right')

    return components(p)

def plot_category(matrix, full_ids, boundaries, joined_data, category, sim_threshold):

    # Flatten matrix
    flatten = MDS(n_components=2)
    flatten_matrix = flatten.fit_transform(matrix)

    # Make df for Bokeh
    bokeh_df = pd.DataFrame({'x': flatten_matrix[:, 0],
                             'y': flatten_matrix[:, 1],
                             'id': list(full_ids),
                             'name': pd.merge(full_ids, joined_data, on='drugbank-id')['name'],
                             'desc': pd.merge(full_ids, joined_data, on='drugbank-id')['description']
                            })

    # Plot
    p = figure(plot_width=1100, plot_height=800, title=f'{category}; Similarity threshold = {sim_threshold}')
    p.axis.visible = False
    p.grid.visible = False
    colors = ['blue', 'orange', 'green', 'red']
    labels = ['Approved', 'Experimental', 'Similar approved', 'Similar experimental']
    legend_items = []
    for i, (stop1, stop2) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        c = p.circle(x='x', y='y', color=colors[i], size=16, source=ColumnDataSource(bokeh_df[stop1:stop2]))
        legend_items.append((f'{labels[i]} ({stop2 - stop1})', [c]))

    p.add_tools(HoverTool(
        tooltips=[
            ('DrugBank ID', '@id'),
            ('Name', '@name'),
            ('Description', '@desc')
        ]))

    legend = Legend(items=legend_items)
    legend.click_policy='hide'

    p.add_layout(legend, 'right')

    return components(p)
