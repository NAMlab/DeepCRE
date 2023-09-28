import os
import pandas as pd
import numpy as np
import pickle
import pyranges as pr
from pyfaidx import Fasta
from tensorflow.keras import Sequential, optimizers, backend, models
from tensorflow.keras.layers import Conv1D, Dense, MaxPool1D, Dropout, Flatten
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.utils import shuffle


def onehot(seq):
    code = {'A': [1, 0, 0, 0],
            'C': [0, 1, 0, 0],
            'G': [0, 0, 1, 0],
            'T': [0, 0, 0, 1],
            'unk': [0, 0, 0, 0]}
    encoded = np.zeros(shape=(len(seq), 4))
    for i, nt in enumerate(seq):
        if nt in ['A', 'C', 'G', 'T']:
            encoded[i, :] = code[nt]
        else:
            encoded[i, :] = code['unk']
    return encoded


class FastaSequenceLoader:
    def __init__(self, fasta, gtf, val_chromosome, pickled_val_ids, pickled_key, upstream=1000, downstream=500,
                 for_prediction=False):
        """
        :param fasta: path to reference genome
        :param gtf: path to gene models
        :param val_chromosome: test chromosome
        :param pickled_val_ids: gene_ids with no homologs in other chroms
        :param upstream: number of nucleotides upstream to extend anchor
        :param downstream: number of nucleotides downstream to extend anchor
        :param for_prediction: whether to use loader for training or prediction
        """
        self.fasta = Fasta(fasta, as_raw=True, sequence_always_upper=True, read_ahead=10000)
        self.upstream = upstream
        self.downstream = downstream
        self.val_chromosome = str(val_chromosome)
        self.pickled_val_ids = pickled_val_ids
        self.pickled_key = pickled_key

        gene_models = pr.read_gtf(gtf, as_df=True)
        gene_models = gene_models[gene_models['Feature'] == 'gene']
        gene_models = gene_models[gene_models['gene_biotype'] == 'protein_coding']
        gene_models = gene_models[['Chromosome', 'Start', 'End', 'Strand', 'gene_id']]
        if for_prediction:
            self.gtf = gene_models[gene_models['Chromosome'] == self.val_chromosome]
        else:
            self.gtf = gene_models

    def extract_seq(self):
        encoded_train_seqs, train_ids = [], []
        encoded_val_seqs, val_ids = [], []
        with open(self.pickled_val_ids, 'rb') as handle:
            validation_genes = pickle.load(handle)
        print(validation_genes.keys())
        for chrom, start, end, strand, gene_id in self.gtf.values:
            if strand == '+':
                prom_start, prom_end = start - self.upstream, start + self.downstream
                term_start, term_end = end - self.downstream, end + self.upstream
                if prom_start > 0 and term_start > 0:
                    encoded_seq = np.concatenate([onehot(self.fasta[chrom][prom_start:prom_end]),
                                                  np.zeros(shape=(20, 4)),
                                                  onehot(self.fasta[chrom][term_start:term_end])])
                    if encoded_seq.shape[0] == 2*(self.upstream + self.downstream) + 20:
                        if chrom == self.val_chromosome:
                            if gene_id in validation_genes[self.pickled_key]:
                                encoded_val_seqs.append(encoded_seq)
                                val_ids.append(gene_id)
                        else:
                            encoded_train_seqs.append(encoded_seq)
                            train_ids.append(gene_id)

            else:
                prom_start, prom_end = end - self.downstream, end + self.upstream
                term_start, term_end = start - self.upstream, start + self.downstream
                if prom_start > 0 and term_start > 0:
                    encoded_seq = np.concatenate([onehot(self.fasta[chrom][prom_start:prom_end])[::-1, ::-1],
                                                  np.zeros(shape=(20, 4)),
                                                  onehot(self.fasta[chrom][term_start:term_end])[::-1, ::-1]])

                    if encoded_seq.shape[0] == 2*(self.upstream + self.downstream) + 20:
                        if chrom == self.val_chromosome:
                            if gene_id in validation_genes[self.pickled_key]:
                                encoded_val_seqs.append(encoded_seq)
                                val_ids.append(gene_id)
                        else:
                            encoded_train_seqs.append(encoded_seq)
                            train_ids.append(gene_id)

        return encoded_train_seqs, encoded_val_seqs, train_ids, val_ids


class ConvNetwork:
    def __init__(self, encoded_train_seqs, encoded_val_seqs, train_ids, val_ids, val_chromosome,
                 tpm_counts, organism, case, outer_flank=1000, inner_flank=500, size_effect=False, tissue=""):
        """
        :param inner_flank: length of sequence inside gene
        :param encoded_train_seqs: onehot encoded train sequences from FastaSequenceLoader
        :param encoded_val_seqs: onehot encoded test sequences from FastaSequenceLoader
        :param train_ids: train gene ids to map true labels to encoded sequences
        :param val_ids: validation gene ids to map true labels to encoded sequences
        :param val_chromosome: chromosome to use for validation
        :param tpm_counts: mapped read counts with associated class
        :param organism: organism name
        :param case: control or treatment
        :param outer_flank: length of sequence upstream and downstream TSS and TTS respectively
        :param size_effect: whether we are investigating the effect of promoter and terminator lengths
        :param tissue: tissue containing expression profile, default is empty string to denote leaf tissue
        """
        self.train_seqs = encoded_train_seqs
        self.val_seqs = encoded_val_seqs
        self.train_ids = train_ids
        self.val_ids = val_ids
        self.val_chrom = val_chromosome
        self.counts = tpm_counts
        self.organism = organism
        self.case = case
        self.outer = outer_flank
        self.effect_size = size_effect
        self.tissue = tissue
        self.inner = inner_flank

    def build_network(self, x_train, x_val, y_train, y_val):
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
        if self.effect_size:
            model_save_name = f'saved_models/size_effect/{self.organism}_{self.val_chrom}_{self.outer}_{self.inner}.h5'
        else:
            model_save_name = f'saved_models/{self.organism}_model_{self.val_chrom}_{self.case}{self.tissue}.h5'
        print(f'save path: {model_save_name}')
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

        if self.effect_size:
            performace = [val_acc, val_auroc, self.organism, self.outer, x_train.shape[0]]
        else:
            performace = [val_acc, val_auroc, self.organism, self.case, x_train.shape[0]]

        if model_save_name.endswith(f'si-nucleotide_shuffle{self.tissue}.h5'):
            os.system(f'rm -rf {model_save_name}')
        return performace

    def train_network(self):
        train_labels, train_seqs = [], []
        val_labels, val_seqs = [], []

        for train_id, train_seq in zip(self.train_ids, self.train_seqs):
            train_labels.append(self.counts.loc[train_id, 'true_target'])
            train_seqs.append(train_seq)
        for val_id, val_seq in zip(self.val_ids, self.val_seqs):
            val_labels.append(self.counts.loc[val_id, 'true_target'])
            val_seqs.append(val_seq)

        train_labels, val_labels = np.array(train_labels), np.array(val_labels)
        train_seqs, val_seqs = np.array(train_seqs), np.array(val_seqs)

        # Random down sampling to balance data
        low_train, high_train = np.where(train_labels == 0)[0], np.where(train_labels == 1)[0]
        min_class = min([len(low_train), len(high_train)])
        selected_low_train = np.random.choice(low_train, min_class, replace=False)
        selected_high_train = np.random.choice(high_train, min_class, replace=False)
        x_train = np.concatenate([
            np.take(train_seqs, selected_low_train, axis=0),
            np.take(train_seqs, selected_high_train, axis=0)
        ], axis=0)
        y_train = np.concatenate([
            np.take(train_labels, selected_low_train, axis=0),
            np.take(train_labels, selected_high_train, axis=0)
        ], axis=0)
        x_train, y_train = shuffle(x_train, y_train, random_state=42)

        # Selecting validation sequences with label 1 and 0
        low_val, high_val = np.where(val_labels == 0)[0], np.where(val_labels == 1)[0]
        x_val = np.concatenate([
            np.take(val_seqs, low_val, axis=0),
            np.take(val_seqs, high_val, axis=0)
        ])
        y_val = np.concatenate([
            np.take(val_labels, low_val, axis=0),
            np.take(val_labels, high_val, axis=0)
        ])

        print(x_train.shape, x_val.shape)
        print(f'validation size: {x_val.shape[0]}')

        # Masking first three nucleotides after gene start and before gene end
        x_train[:, self.outer:self.outer + 3, :] = 0
        x_train[:, self.outer + (self.inner * 2) + 17:self.outer + (self.inner * 2) + 20, :] = 0
        x_val[:, self.outer:self.outer + 3, :] = 0
        x_val[:, self.outer + (self.inner * 2) + 17:self.outer + (self.inner * 2) + 20, :] = 0

        output = self.build_network(x_train, x_val, y_train, y_val)

        return output


def prepare_valid_seqs(fasta, gtf, tpms, val_chrom, pkey=False, upstream=1000, downstream=500):
    fasta = Fasta(f'genomes/{fasta}', as_raw=True, sequence_always_upper=True, read_ahead=10000)
    gene_models = pr.read_gtf(f'gene_models/{gtf}', as_df=True)
    gene_models = gene_models[gene_models['Feature'] == 'gene']
    gene_models = gene_models[gene_models['gene_biotype'] == 'protein_coding']
    gene_models = gene_models[['Chromosome', 'Start', 'End', 'Strand', 'gene_id']]
    if isinstance(val_chrom, list):
        gene_models = gene_models[gene_models['Chromosome'].isin(val_chrom)]
    else:
        gene_models = gene_models[gene_models['Chromosome'].isin([val_chrom])]
    if pkey:
        # Pickled validation IDs
        with open('validation_genes.pickle', 'rb') as handle:
            validation_genes = pickle.load(handle)
        gene_models = gene_models[gene_models['gene_id'].isin(validation_genes[pkey])]

    # Transcripts per Million
    tpm_counts = pd.read_csv(f'tpm_counts/{tpms}', index_col=0)
    true_targets = []

    for log_count in tpm_counts['logMaxTPM'].values:
        if log_count <= np.percentile(tpm_counts['logMaxTPM'], 25):
            true_targets.append(0)
        elif log_count >= np.percentile(tpm_counts['logMaxTPM'], 75):
            true_targets.append(1)
        else:
            true_targets.append(2)
    tpm_counts['true_target'] = true_targets

    encoded_val_seqs, labels, gene_ids = [], [], []
    for chrom, start, end, strand, gene_id in gene_models.values:
        if strand == '+':
            prom_start, prom_end = start - upstream, start + downstream
            term_start, term_end = end - downstream, end + upstream
            if prom_start > 0 and term_start > 0:
                encoded_seq = np.concatenate([onehot(fasta[chrom][prom_start:prom_end]),
                                              np.zeros(shape=(20, 4)),
                                              onehot(fasta[chrom][term_start:term_end])])
                if encoded_seq.shape[0] == 2 * (upstream + downstream) + 20:
                    encoded_val_seqs.append(encoded_seq)
                    labels.append(tpm_counts.loc[gene_id, 'true_target'])
                    gene_ids.append(gene_id)

        else:
            prom_start, prom_end = end - downstream, end + upstream
            term_start, term_end = start - upstream, start + downstream
            if prom_start > 0 and term_start > 0:
                encoded_seq = np.concatenate([onehot(fasta[chrom][prom_start:prom_end])[::-1, ::-1],
                                              np.zeros(shape=(20, 4)),
                                              onehot(fasta[chrom][term_start:term_end])[::-1, ::-1]])

                if encoded_seq.shape[0] == 2 * (upstream + downstream) + 20:
                    encoded_val_seqs.append(encoded_seq)
                    labels.append(tpm_counts.loc[gene_id, 'true_target'])
                    gene_ids.append(gene_id)

    # Selecting validation sequences with label 1 and 0
    labels, encoded_val_seqs, gene_ids = np.array(labels), np.array(encoded_val_seqs), np.array(gene_ids)
    low_val, high_val = np.where(labels == 0)[0], np.where(labels == 1)[0]
    x_val = np.concatenate([
        np.take(encoded_val_seqs, low_val, axis=0),
        np.take(encoded_val_seqs, high_val, axis=0)
    ])
    y_val = np.concatenate([
        np.take(labels, low_val, axis=0),
        np.take(labels, high_val, axis=0)
    ])
    gene_ids = np.concatenate([np.take(gene_ids, low_val, axis=0), np.take(gene_ids, high_val, axis=0)])
    x_val[:, upstream:upstream+3, :] = 0
    x_val[:, upstream+(downstream*2)+17:upstream+(downstream*2)+20, :] = 0
    return x_val, y_val, gene_ids

