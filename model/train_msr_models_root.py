"""In this script we train models with several plant species and test on single species"""
import pandas as pd
import numpy as np
from utils import onehot
from pyfaidx import Fasta
from sklearn.utils import shuffle
from tensorflow.keras import Sequential, optimizers, backend, models
from tensorflow.keras.layers import Conv1D, Dense, MaxPool1D, Dropout, Flatten
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
from sklearn.metrics import accuracy_score, roc_auc_score
import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'


gene_model = pd.read_csv('gene_models/zea_sol_ara_sor_52.gtf', sep='\t')
gene_model['specie'] = [x.split('_')[0] for x in gene_model['Chromosome']]
print(gene_model.head())
genes_labels = pd.read_csv('tpm_counts/zea_sol_ara_sor_roots.csv', index_col=0)
print(genes_labels.head())
genome = Fasta('genomes/zea_sol_ara_sor_dna.fa',  as_raw=True, sequence_always_upper=True, read_ahead=10000)


def encode_seq(gene_models, fasta, labels, upstream=1000, downstream=500, training=False):
    encoded_seqs, true_labels = [], []

    for chrom, start, end, strand, gene_id, _ in gene_models.values:
        if gene_id in labels.index:
            if strand == '+':
                prom_start, prom_end = start - upstream, start + downstream
                term_start, term_end = end - downstream, end + upstream
                if prom_start > 0 and term_start > 0:
                    encoded_seq = np.concatenate([onehot(fasta[chrom][prom_start:prom_end]),
                                                  np.zeros(shape=(20, 4)),
                                                  onehot(fasta[chrom][term_start:term_end])])
                    if encoded_seq.shape[0] == 2 * (upstream + downstream) + 20:
                        encoded_seqs.append(encoded_seq)
                        true_labels.append(labels.loc[gene_id, 'true_target'])

            else:
                prom_start, prom_end = end - downstream, end + upstream
                term_start, term_end = start - upstream, start + downstream
                if prom_start > 0 and term_start > 0:
                    encoded_seq = np.concatenate([onehot(fasta[chrom][prom_start:prom_end])[::-1, ::-1],
                                                  np.zeros(shape=(20, 4)),
                                                  onehot(fasta[chrom][term_start:term_end])[::-1, ::-1]])

                    if encoded_seq.shape[0] == 2 * (upstream + downstream) + 20:
                        encoded_seqs.append(encoded_seq)
                        true_labels.append(labels.loc[gene_id, 'true_target'])

    encoded_seqs = np.array(encoded_seqs)
    true_labels = np.array(true_labels)

    if training:
        print('Encoding and balancing Training set')
        # Random down sampling to balance data
        low_train, high_train = np.where(true_labels == 0)[0], np.where(true_labels == 1)[0]
        min_class = min([len(low_train), len(high_train)])
        selected_low_train = np.random.choice(low_train, min_class, replace=False)
        selected_high_train = np.random.choice(high_train, min_class, replace=False)
        encoded_seqs = np.concatenate([
            np.take(encoded_seqs, selected_low_train, axis=0),
            np.take(encoded_seqs, selected_high_train, axis=0)
        ], axis=0)
        true_labels = np.concatenate([
            np.take(true_labels, selected_low_train, axis=0),
            np.take(true_labels, selected_high_train, axis=0)
        ], axis=0)
        encoded_seqs, true_labels = shuffle(encoded_seqs, true_labels, random_state=42)

    return encoded_seqs, true_labels


def build_network(x_train, x_val, y_train, y_val, specie_name, val_chrom):
    backend.clear_session()
    model = Sequential([
        # Conv Block 1
        Conv1D(64, kernel_size=8, activation='relu', padding='same',
               input_shape=(x_train.shape[1], x_train.shape[2])),
        Conv1D(64, kernel_size=8, activation='relu', padding='same'),
        MaxPool1D(8, padding='same'),
        Dropout(0.25),

        # Conv Block 2
        Conv1D(128, kernel_size=8, activation='relu', padding='same'),
        Conv1D(128, kernel_size=8, activation='relu', padding='same'),
        MaxPool1D(8, padding='same'),
        Dropout(0.25),

        # Conv Block 3
        Conv1D(64, kernel_size=8, activation='relu', padding='same'),
        Conv1D(64, kernel_size=8, activation='relu', padding='same'),
        MaxPool1D(8, padding='same'),
        Dropout(0.25),

        # Fully connected Block
        Flatten(),
        Dense(128, activation='relu'),
        Dropout(0.25),
        Dense(64, activation='relu'),
        Dense(1, activation='sigmoid')

    ])

    print(model.summary())

    model_save_name = f'saved_models/multi_specie_model_{specie_name}_{val_chrom}_root.h5'
    model_chkpt = ModelCheckpoint(model_save_name, save_best_only=True, verbose=1)
    early_stop = EarlyStopping(patience=10)
    reduce_lr = ReduceLROnPlateau(patience=5, factor=0.1)
    model.compile(loss='binary_crossentropy', optimizer=optimizers.Adam(0.0001), metrics=['accuracy'])
    model.fit(x_train, y_train, batch_size=64, epochs=100, validation_data=(x_val, y_val),
              callbacks=[early_stop, model_chkpt, reduce_lr])
    saved_model = models.load_model(model_save_name)
    predictions = saved_model.predict(x_val)
    val_auroc = roc_auc_score(y_val, predictions)
    predictions = predictions > 0.5
    val_acc = accuracy_score(y_val, predictions)
    print('Best model performance--------------------------\n')
    print(f'Accuracy: {val_acc}, auROC: {val_auroc}\n')
    print('------------------------------------------------')

    performace = [val_acc, val_auroc, specie_name]

    return performace


final_results = []
for specie in gene_model['specie'].unique():
    test_specie = gene_model.copy()
    test_specie = test_specie[test_specie['specie'] == specie]
    train_specie = gene_model.copy()
    train_specie = train_specie[train_specie['specie'] != specie]
    train_seqs, train_targets = encode_seq(train_specie, genome, genes_labels, training=True)

    for chrom in test_specie['Chromosome'].unique():
        test_specie_chrom = test_specie.copy()
        test_specie_chrom = test_specie_chrom[test_specie_chrom['Chromosome'] == chrom]
        test_seqs, test_targets = encode_seq(test_specie_chrom, genome, genes_labels)
        print(f'Training on {train_specie["specie"].unique()}')
        print(f'Testing on {test_specie["specie"].unique()}, chromosome {chrom}')
        result = build_network(train_seqs, test_seqs, train_targets, test_targets, specie, chrom)
        final_results.append(result)

final_results = pd.DataFrame(final_results, columns=['accuracy', 'auRoc', 'test specie'])
final_results.to_csv("../results/multi_specie_result_root.csv", index=False)







