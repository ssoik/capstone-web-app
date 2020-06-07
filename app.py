from flask import Flask, render_template, request, redirect
import requests
import re
import pandas as pd
from functions import find, matrix, plot

url = 'https://tdi-capstone-web-app-data.s3.us-east-2.amazonaws.com/'

app = Flask(__name__)

@app.route('/')
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
    
    # Request target
    target = request.form['target']
    sim_thresh = float(request.form['threshold'])
    fname = f'{target}.txt'
    
    # Load data
    inds = pd.read_csv(f'{url}sim-data/{fname}')
    full_data = pd.read_csv(f'{url}drugbank-data/full_database.csv')
    joined_data = pd.merge(inds, full_data, on='drugbank-id')

    # Select approved and experimental drugs
    approved_drugs = find.similar_approved([target], [], sim_thresh, joined_data, url)
    experimental_drugs = find.similar_experimental([target], [], sim_thresh, joined_data, url)

    # Construct graph
    full_ids = pd.Series(pd.concat([pd.Series(target), approved_drugs, experimental_drugs]), name='drugbank-id')
    sim_matrix = matrix.similarity_matrix(full_ids, url)

    # Determine boundaries by status
    bounds = [0, 1]
    for df in [approved_drugs, experimental_drugs]:
        bounds.append(len(df) + bounds[-1])

    # Render HTML template
    script, div = plot.plot_target(sim_matrix, full_ids, bounds, joined_data, target, sim_thresh)
    return render_template('display_target.html', script=script, div=div)

@app.route('/display_category', methods = ['POST'])
def display_category():

    # Request category
    category = request.form['category']
    fname = f"{'-'.join(re.split(' |/', category))}"

    # Get components from cache and render HTML template
    components = requests.get(f'{url}category-cache/{fname}.txt', 'r').text.split(';;;;;')
    return render_template('display_category.html', script=components[0], div=components[1])


    """The original code for these computations follows. Visual components were then cached to improve app performance.
    
    # Load data
    full_data = pd.read_csv(f'{url}drugbank-data/full_database.csv')
    struct_data = pd.read_csv(f'{url}drugbank-data/structure_links_clean.csv')
    joined_data = pd.merge(full_data, struct_data, left_on='drugbank-id', right_on='DrugBank ID')

    # Select approved and experimental drugs
    approved_drugs = find.approved(category, joined_data)
    experimental_drugs = find.experimental(category, joined_data)

    # Set similarity threshold
    approved_sim_matrix = matrix.similarity_matrix(approved_drugs['drugbank-id'], url)
    sim_thresh = matrix.similarity_threshold(approved_sim_matrix, approved_drugs['drugbank-id'], url)

    # Find similar drugs
    sim_approved_drug_ids = find.similar_approved(approved_drugs['drugbank-id'], \
        experimental_drugs['drugbank-id'], sim_thresh, joined_data, url)
    sim_experimental_drug_ids = find.similar_experimental(approved_drugs['drugbank-id'], \
        experimental_drugs['drugbank-id'], sim_thresh, joined_data, url)

    # Construct full graph
    full_ids = pd.Series(pd.concat([approved_drugs['drugbank-id'], experimental_drugs['drugbank-id'], \
        sim_approved_drug_ids, sim_experimental_drug_ids]), name='drugbank-id')
    full_sim_matrix = matrix.similarity_matrix(full_ids, url)

    # Determine boundaries by status
    bounds = [0]
    for df in [approved_drugs, experimental_drugs, sim_approved_drug_ids, sim_experimental_drug_ids]:
        bounds.append(len(df) + bounds[-1])

    # Get plot components and render HTML
    script, div = plot.plot_category(full_sim_matrix, full_ids, bounds, joined_data, category, sim_thresh)
    return render_template('display_category.html', script=script, div=div)
    """


if __name__ == '__main__':
    app.run(port=5000, debug=True)
