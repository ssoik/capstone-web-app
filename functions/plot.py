# Plot flattened graph and output HTML

from sklearn.manifold import MDS
import pandas as pd
import numpy as np
from bokeh.models import ColumnDataSource, HoverTool, Legend, CustomJS, Slider
from bokeh.plotting import figure
from bokeh.layouts import row, column
from bokeh.embed import components

def plot_target(matrix, full_ids, boundaries, data, target, sim_threshold):

    # Flatten matrix
    flatten = MDS(n_components=2)
    flatten_matrix = flatten.fit_transform(matrix)
    scaled_coords = flatten_matrix / np.sqrt(2) / (1 - sim_threshold)

    # Make df for Bokeh
    merged_data = pd.merge(full_ids, data, on='drugbank-id')
    source = pd.DataFrame({
        'x': flatten_matrix[:, 0],
        'y': flatten_matrix[:, 1],
        'id': np.array(full_ids),
        'index': np.array(merged_data['dice-similarity']),
        'name': np.array(merged_data['name']),
        'desc': np.array(merged_data['description'])
    })

    # Plot
    p = figure(plot_width=900, plot_height=650, title=f"{source[source['id'] == target]['name'].item()}; Similarity threshold = {sim_threshold}", \
        toolbar_location='below')
    p.axis.visible = False
    p.grid.visible = False
    colors = ['blue', 'green', 'red']
    labels = ['Target drug', 'Similar approved', 'Similar experimental']
    legend_items = []
    for i, (stop1, stop2) in enumerate(zip(boundaries[:-1][::-1], boundaries[1:][::-1])):
        c = p.circle('x', 'y', color=colors[::-1][i], size=16, source=source[stop1:stop2])
        legend_items = [(f"{labels[::-1][i]} ({stop2 - stop1})", [c])] + legend_items

    p.add_tools(HoverTool(
        tooltips=[
            ('DrugBank ID', '@id'),
            ('Name', '@name'),
            ('Similarity index', '@index'),
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
    merged_data = pd.merge(full_ids, joined_data, on='drugbank-id')
    source = pd.DataFrame({
        'x': flatten_matrix[:, 0],
        'y': flatten_matrix[:, 1],
        'id': np.array(full_ids),
        'name': np.array(merged_data['name']),
        'desc': np.array(merged_data['description'])
    })

    # Plot
    p = figure(plot_width=900, plot_height=650, title=f'{category}; Similarity threshold = {sim_threshold}', \
        toolbar_location='below')
    p.axis.visible = False
    p.grid.visible = False
    colors = ['blue', 'orange', 'green', 'red']
    labels = ['Approved', 'Experimental', 'Similar approved', 'Similar experimental']
    legend_items = []
    for i, (stop1, stop2) in enumerate(zip(boundaries[:-1][::-1], boundaries[1:][::-1])):
        c = p.circle('x', 'y', color=colors[::-1][i], size=16, source=source[stop1:stop2])
        legend_items = [(f'{labels[::-1][i]} ({stop2 - stop1})', [c])] + legend_items

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
