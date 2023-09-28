import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import RandomOverSampler
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
pd.options.display.width = 0

feat_imp_scores, species, feat_names = [], [], []
scores = []

for cls_task in ['Low-medium-High', 'Low-High']:
    mapped_read_counts = ['arabidopsis_counts.csv', 'zea_counts.csv', 'solanum_counts.csv', 'sbicolor_counts.csv']
    generic_feats = ['Arabidopsis_generated_features.csv', 'Zea_generated_features.csv',
                     'Solanum_generated_features.csv', 'Sorghum_generated_features.csv']

    for idx in range(4):
        tpm = pd.read_csv(f'tpm_counts/{mapped_read_counts[idx]}')
        predictors = pd.read_csv(generic_feats[idx], index_col=0)
        targets = []
        for log_count in tpm['logMaxTPM']:
            if log_count <= np.percentile(tpm['logMaxTPM'], 25):
                targets.append(0)
            elif log_count >= np.percentile(tpm['logMaxTPM'], 75):
                targets.append(1)
            else:
                targets.append(2)
        tpm['label'] = targets
        tpm = tpm[['gene_id', 'label']]
        data = predictors.merge(tpm, how='inner', on='gene_id')
        if cls_task == 'Low-High':
            data = data[data['label'] != 2]
        print(data.head())

        for chrom in data['Chromosome'].unique():
            data_train = data.copy()
            data_train = data_train[data_train['Chromosome'] != chrom]
            data_train.drop(columns=['gene_id', 'Chromosome'], inplace=True)
            data_test = data.copy()
            data_test = data_test[data_test['Chromosome'] == chrom]
            data_test.drop(columns=['gene_id', 'Chromosome'], inplace=True)

            # Balance data and  standardizing
            sampler = RandomOverSampler(random_state=42)
            x_train, y_train = data_train.values[:, :-1], data_train['label'].values
            x_train, y_train = sampler.fit_resample(x_train, y_train)
            x_test, y_test = data_test.values[:, :-1], data_test['label'].values
            scaler = StandardScaler()
            scaler.fit(x_train)
            x_train_std = scaler.transform(x_train)
            x_test_std = scaler.transform(x_test)
            random_forest = RandomForestClassifier(100)
            random_forest.fit(x_train_std, y_train)
            if cls_task == 'Low-High':
                feat_imp = random_forest.feature_importances_
                feat_imp_scores.extend(feat_imp)
                feat_names.extend(data_train.columns[:-1])
                species.extend([generic_feats[idx].split('_')[0]] * len(data_train.columns[:-1]))
            y_pred = random_forest.predict(x_test_std)
            acc = accuracy_score(y_test, y_pred)
            scores.append([generic_feats[idx].split('_')[0], acc, cls_task])
            print(generic_feats[idx].split('_')[0], random_forest.score(x_test_std, y_test), cls_task)

df_feat_imp = pd.DataFrame({'Importance scores': feat_imp_scores,
                            'Species': species,
                            'Feature': feat_names})
df_feat_imp.sort_values(by='Species', inplace=True)
df_performance = pd.DataFrame(scores, columns=['Species', 'accuracy',
                                               'task']).to_csv('../results/rand_for_perf.csv', index=False)
df_feat_imp.to_csv('../results/rand_for_feat_imp_SSR.csv', index=False)
