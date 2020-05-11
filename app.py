from flask import Flask, render_template, request, redirect
import os
import pandas as pd
from functions import find, matrix, plot

filepath = '~/Documents/TDI/capstone-project/data/'

app = Flask(__name__)

@app.route('/search')
def search():
    return render_template('search.html')

@app.route('/about')
def about():

    # Get components to plot stats
    target_stats = open('target_stats.txt', 'r')
    target_comps = target_stats.read().split(';;;;;')
    category_stats = open('category_stats.txt', 'r')
    category_comps = category_stats.read().split(';;;;;')
    return render_template('about.html', target_script=target_comps[0], target_div=target_comps[1], \
        category_script=category_comps[0], category_div=category_comps[1])

@app.route('/display_target', methods = ['POST'])
def display_target():
    
    # Request target and check if script, div exist
    target = request.form['target']
    sim_thresh = float(request.form['threshold'])
    fname = f"{'-'.join(target.split(' '))}.txt"
    """
    if os.path.exists(f'target-plots-cache/{fname}'):
        text = open(f'target-plots-cache/{fname}', 'r')
        text_parts = text.read().split(';;;;;')
        return render_template('display_target.html', script=text_parts[0], div=text_parts[1])
    """
    # Load data
    inds = pd.read_csv(f'{filepath}sim-data/{fname}')
    full_data = pd.read_csv(f'{filepath}full_database.csv')
    joined_data = pd.merge(inds, full_data, on='drugbank-id')

    # Select approved and experimental drugs
    approved_drugs = find.similar_approved([target], [], sim_thresh, joined_data, filepath)
    experimental_drugs = find.similar_experimental(pd.concat([pd.Series(target), approved_drugs]), [], \
        sim_thresh, joined_data, filepath)

    # Construct graph
    full_ids = pd.Series(pd.concat([pd.Series(target), approved_drugs, experimental_drugs]), name='drugbank-id')
    sim_matrix = matrix.similarity_matrix(full_ids, filepath)

    # Determine boundaries by status
    bounds = [0, 1]
    for df in [approved_drugs, experimental_drugs]:
        bounds.append(len(df) + bounds[-1])

    # Render HTML template and cache
    script, div = plot.plot_target(sim_matrix, full_ids, bounds, joined_data, target, sim_thresh)
    cache = open(f'target-plots-cache/{fname}', 'w')
    cache.write(f'{script};;;;;{div}')
    cache.close()
    return render_template('display_target.html', script=script, div=div)

@app.route('/display_category', methods = ['POST'])
def display_category():

    # Request category and check if script, div exist
    category = request.form['category']
    fname = f"{'-'.join(category.split(' '))}.txt"

    if os.path.exists(f'category-plots-cache/{fname}'):
        text = open(f'category-plots-cache/{fname}', 'r')
        text_parts = text.read().split(';;;;;')
        return render_template('display_category.html', script=text_parts[0], div=text_parts[1])

    # Load data
    full_data = pd.read_csv(f'{filepath}full_database.csv')
    struct_data = pd.read_csv(f'{filepath}structure_links_clean.csv')
    joined_data = pd.merge(full_data, struct_data, left_on='drugbank-id', right_on='DrugBank ID')

    # Select approved and experimental drugs
    approved_drugs = find.approved(category, joined_data)
    experimental_drugs = find.experimental(category, joined_data)

    # Set similarity threshold
    approved_sim_matrix = matrix.similarity_matrix(approved_drugs['drugbank-id'], filepath)
    sim_thresh = matrix.similarity_threshold(approved_sim_matrix, approved_drugs['drugbank-id'], filepath)

    # Find similar drugs
    sim_approved_drug_ids = find.similar_approved(approved_drugs['drugbank-id'], \
        experimental_drugs['drugbank-id'], sim_thresh, joined_data, filepath)
    sim_experimental_drug_ids = find.similar_experimental(approved_drugs['drugbank-id'], \
        experimental_drugs['drugbank-id'], sim_thresh, joined_data, filepath)

    # Construct full graph
    full_ids = pd.Series(pd.concat([approved_drugs['drugbank-id'], experimental_drugs['drugbank-id'], \
        sim_approved_drug_ids, sim_experimental_drug_ids]), name='drugbank-id')
    full_sim_matrix = matrix.similarity_matrix(full_ids, filepath)

    # Determine boundaries by status
    bounds = [0]
    for df in [approved_drugs, experimental_drugs, sim_approved_drug_ids, sim_experimental_drug_ids]:
        bounds.append(len(df) + bounds[-1])

    # Render HTML template and cache
    script, div = plot.plot_category(full_sim_matrix, full_ids, bounds, joined_data, category, sim_thresh)
    cache = open(f'category-plots-cache/{fname}', 'w')
    cache.write(f'{script};;;;;{div}')
    cache.close()
    return render_template('display_category.html', script=script, div=div)

if __name__ == '__main__':
    app.run(port=5000)
