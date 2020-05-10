# Select appropriate rows from joined_data

import pandas as pd

def approved(category, joined_data):
    return joined_data[joined_data['categories'].str.contains(category, na=False) & \
                       joined_data['groups'].str.contains('approved', na=False)]

def experimental(category, joined_data):
    return joined_data[joined_data['categories'].str.contains(category, na=False) & \
                       ~joined_data['groups'].str.contains('approved', na=True)]

def similar_approved(approved_ids, experimental_ids, sim_threshold, joined_data, filepath):
    other_approved_ids = []
    for drug in approved_ids:
        inds = pd.read_csv(f'{filepath}sim-data/{drug}.txt')
        other_approved_ids += [db_id for db_id in pd.merge(inds['drugbank-id'][inds['dice-similarity'] > \
                                sim_threshold], \
                                joined_data[joined_data['groups'].str.contains('approved', na=False)], \
                                on='drugbank-id')['drugbank-id'] if \
                                db_id not in (list(approved_ids) + list(experimental_ids))]

    return pd.Series(list(set(other_approved_ids)), name='drugbank-id')

def similar_experimental(approved_ids, experimental_ids, sim_threshold, joined_data, filepath):
    other_experimental_ids = []
    for drug in approved_ids:
        inds = pd.read_csv(f'{filepath}sim-data/{drug}.txt')
        other_experimental_ids += [db_id for db_id in pd.merge(inds['drugbank-id'][inds['dice-similarity'] > \
                                    sim_threshold], \
                                    joined_data[~joined_data['groups'].str.contains('approved', na=True)], \
                                    on='drugbank-id')['drugbank-id'] if \
                                    db_id not in (list(approved_ids) + list(experimental_ids))]

    return pd.Series(list(set(other_experimental_ids)), name='drugbank-id')
