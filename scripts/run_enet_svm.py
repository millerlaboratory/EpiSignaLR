#!/usr/bin/env python3
"""
Score sample(s) with the ElasticNet-selected-feature linear-SVM episignature
classifier and generate a PCA plot and a heatmap of the new sample(s) against
the training cohort.

Expects two model artifacts in --model-dir (produced by
nbs_epi_data/SETD1B_LRS/Corrected_BC/loo_SVM_enet_implementation.py):
    enet_svm_model.joblib      full pipeline: StandardScaler -> ElasticNetSelector -> linear SVM
    enet_svm_metadata.json
plus a training-data snapshot, --training-data (or <model-dir>/enet_svm_training_data.csv),
containing at minimum a sample-ID column, a label column, and EVERY DMR column
the model was trained on (e.g. a copy of the discovery-cohort CSV used to train
the model) -- not just the final selected panel. The pipeline's StandardScaler
and ElasticNetSelector were fit on the full DMR panel, so both the training-data
snapshot and the input region matrix must cover that same full panel; scoring
against only the final selected DMRs will fail (missing-column error) or, if it
didn't, would silently misalign the scaler/selector.

Input region matrix (--region-matrix) is the CSV produced by
10_region_average_methylation.py: rows are regions, columns are samples.

Usage:
    python run_enet_svm.py \
        --region-matrix sample.SETD1B_enet_svm_full_panel_DMPs.csv \
        --model-dir SETD1B_enet_svm_model/ \
        --out-prefix sample.svm_results
"""

import argparse
import json

import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression


# joblib pickled the fitted pipeline's 'select' step as '__main__.ElasticNetSelector'
# because loo_SVM_enet_implementation.py was run directly as a script. Unpickling
# looks up that class by module + name, so it must exist under __main__ here too --
# this definition must stay byte-for-byte in sync with the one in that script.
class ElasticNetSelector(BaseEstimator, TransformerMixin):
    def __init__(self, C=0.5, l1_ratio=0.5, max_iter=20000,
                 min_features=5, max_features=30, random_state=42):
        self.C = C
        self.l1_ratio = l1_ratio
        self.max_iter = max_iter
        self.min_features = min_features
        self.max_features = max_features
        self.random_state = random_state

    def fit(self, X, y):
        self.est_ = LogisticRegression(
            penalty='elasticnet', solver='saga',
            C=self.C, l1_ratio=self.l1_ratio,
            max_iter=self.max_iter, random_state=self.random_state,
        )
        self.est_.fit(X, y)
        coef = np.abs(self.est_.coef_[0])
        nonzero = np.where(coef > 1e-8)[0]

        if len(nonzero) < self.min_features:
            keep = np.argsort(coef)[::-1][:self.min_features]
        elif len(nonzero) > self.max_features:
            order = nonzero[np.argsort(coef[nonzero])[::-1]]
            keep = order[:self.max_features]
        else:
            keep = nonzero

        self.support_ = np.zeros(X.shape[1], dtype=bool)
        self.support_[keep] = True
        return self

    def transform(self, X):
        return X[:, self.support_]

    def get_support(self):
        return self.support_


# Validated palette (dataviz skill, references/palette.md) for the PCA plot
COLOR_CONTROL = '#2a78d6'   # categorical slot 1 (blue)
COLOR_CASE = '#e34948'      # categorical slot 6 (red)
COLOR_GRID = '#e1e0d9'

# pheatmap-style palette for the clustered heatmap
HEATMAP_CONTROL = '#2E86AB'
HEATMAP_CASE = '#E84855'
HEATMAP_NEW = '#4d4d4d'
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    'diverging_beta', ['#2166AC', '#F7F7F7', '#D6604D'])


def load_model(model_dir):
    pipeline = joblib.load(model_dir + 'enet_svm_model.joblib')
    with open(model_dir + 'enet_svm_metadata.json') as f:
        meta = json.load(f)
    print(f"Model loaded: {meta['classifier_type']} | "
          f"{meta['final_model']['n_features']} selected features "
          f"| best params={meta['final_model']['best_params']}")
    return pipeline, meta


def load_training_data(training_data_path, meta, sample_col, label_col):
    df = pd.read_csv(training_data_path)
    dmr_cols = [c for c in df.columns if c not in (sample_col, label_col)]

    train_ids = set(meta['training_cohort']['sample_ids'])
    df = df[df[sample_col].astype(str).isin(train_ids)].reset_index(drop=True)
    if len(df) != meta['training_cohort']['n_samples']:
        print(f"WARNING: training-data snapshot has {len(df)} matching samples, "
              f"metadata records {meta['training_cohort']['n_samples']}.")

    sample_ids = df[sample_col].astype(str).values
    y = df[label_col].astype(int).values
    X = df[dmr_cols].values.astype(float)
    return sample_ids, y, X, dmr_cols


def load_region_matrix(region_matrix_path, dmr_cols):
    """Region matrix is regions x samples; transpose to samples x DMRs and
    reindex to the FULL DMR panel the pipeline's scaler/selector were fit on
    (in that same column order), not just the final selected panel."""
    region_df = pd.read_csv(region_matrix_path, index_col=0)
    sample_df = region_df.T

    missing_dmrs = [d for d in dmr_cols if d not in sample_df.columns]
    if missing_dmrs:
        raise ValueError(
            f"{len(missing_dmrs)} required DMR(s) missing from region matrix "
            f"{region_matrix_path}: {missing_dmrs[:5]}\n"
            f"Ensure the region matrix covers the full DMR panel the ElasticNet "
            f"+ SVM pipeline was trained on (not just the final selected panel) "
            f"and that the same DMR BED file / pipeline version were used."
        )

    sample_ids = sample_df.index.astype(str).values
    X = sample_df[dmr_cols].values.astype(float)
    return sample_ids, X


def score_samples(pipeline, meta, X_new):
    if np.isnan(X_new).any():
        n_missing = int(np.isnan(X_new).sum())
        print(f"WARNING: {n_missing} missing DMR value(s) in new sample(s) "
              f"-- imputing with training-set means.")
        scaler = pipeline.named_steps['scaler']
        inds = np.where(np.isnan(X_new))
        X_new[inds] = np.take(scaler.mean_, inds[1])

    scores = pipeline.predict_proba(X_new)[:, 1]
    threshold = meta['scoring']['threshold']
    indet_lo, indet_hi = meta['scoring']['indeterminate_zone']
    preds = (scores >= threshold).astype(int)
    indet = (scores >= indet_lo) & (scores <= indet_hi)
    return scores, preds, indet


def write_results(sample_ids, scores, preds, indet, meta, out_prefix):
    results = pd.DataFrame({
        'sample': sample_ids,
        'score': np.round(scores, 4),
        'predicted_class': ['Case' if p == 1 else 'Control' for p in preds],
        'indeterminate': indet,
    })
    score_csv = f"{out_prefix}_score.csv"
    results.to_csv(score_csv, index=False)
    print(f"Saved: {score_csv}")

    indet_lo, indet_hi = meta['scoring']['indeterminate_zone']
    top_txt = f"{out_prefix}_top_result.txt"
    with open(top_txt, 'w') as f:
        for _, row in results.iterrows():
            f.write(f"Sample: {row['sample']}\n")
            f.write(f"Score: {row['score']}\n")
            f.write(f"Assigned Class: {row['predicted_class']}\n")
            if row['indeterminate']:
                f.write(f"NOTE: Score falls in the indeterminate zone "
                        f"({indet_lo}-{indet_hi}) -- expert review recommended.\n")
            f.write("\n")
    print(f"Saved: {top_txt}")

    if indet.any():
        print(f"WARNING: {int(indet.sum())} sample(s) in indeterminate zone "
              f"({indet_lo}-{indet_hi}) -- expert review recommended.")

    return results


def plot_pca(train_ids, y_train, X_train, new_ids, X_new, results, out_prefix):
    scaler_pca = PCA(n_components=2)
    X_train_pca = scaler_pca.fit_transform(X_train)
    X_new_pca = scaler_pca.transform(X_new)
    var_exp = scaler_pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(8, 6))
    colours = {0: COLOR_CONTROL, 1: COLOR_CASE}
    labels = {0: 'Training Control', 1: 'Training Case'}

    for label in [0, 1]:
        mask = y_train == label
        ax.scatter(X_train_pca[mask, 0], X_train_pca[mask, 1],
                   c=colours[label], s=90, label=labels[label],
                   edgecolors='white', linewidth=0.8, zorder=3)
        for xi, yi, sid in zip(X_train_pca[mask, 0], X_train_pca[mask, 1], train_ids[mask]):
            ax.annotate(str(sid), (xi, yi), textcoords='offset points',
                        xytext=(0, 7), fontsize=6, ha='center', color=colours[label])

    for i, (xi, yi, sid) in enumerate(zip(X_new_pca[:, 0], X_new_pca[:, 1], new_ids)):
        pred = results.loc[results['sample'] == sid, 'predicted_class'].iloc[0]
        indet = bool(results.loc[results['sample'] == sid, 'indeterminate'].iloc[0])
        pt_col = COLOR_CASE if pred == 'Case' else COLOR_CONTROL
        edge = '#eda100' if indet else '#0b0b0b'
        ax.scatter(xi, yi, c=pt_col, s=150, marker='D', zorder=4,
                   edgecolors=edge, linewidth=1.6)
        ax.annotate(str(sid), (xi, yi), textcoords='offset points',
                    xytext=(0, 9), fontsize=7, ha='center', fontweight='bold', color='#0b0b0b')

    ax.set_xlabel(f'PC1: {var_exp[0]:.1f}%', color='#52514e')
    ax.set_ylabel(f'PC2: {var_exp[1]:.1f}%', color='#52514e')
    ax.set_title('PCA -- Training Cohort + New Sample(s) (ElasticNet-Selected DMRs)', fontweight='bold')
    ax.grid(alpha=0.3, color=COLOR_GRID)

    new_sample_handle = Line2D([0], [0], marker='D', color='w', markerfacecolor='#898781',
                                markeredgecolor='#0b0b0b', markersize=9, label='New Sample')
    indet_handle = Line2D([0], [0], marker='D', color='w', markerfacecolor='#898781',
                           markeredgecolor='#eda100', markersize=9, label='New Sample (indeterminate)')
    handles, lbls = ax.get_legend_handles_labels()
    ax.legend(handles=handles + [new_sample_handle, indet_handle],
              labels=lbls + ['New Sample', 'New Sample (indeterminate)'], fontsize=8, loc='best')

    plt.tight_layout()
    out_svg = f"{out_prefix}_pca.svg"
    plt.savefig(out_svg, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out_svg}")


def plot_heatmap(train_ids, y_train, X_train, new_ids, X_new, selected_dmrs, out_prefix):
    """pheatmap-style clustered heatmap: DMRs as rows, samples as columns,
    hierarchical clustering on both axes, diverging colour scale, and a
    Condition column-annotation bar (Control / Case / New Sample). Expects
    raw (unscaled, 0-1 methylation fraction) values already restricted to the
    ElasticNet-selected panel, with any missing new-sample values imputed."""
    train_df = pd.DataFrame(X_train.T, index=selected_dmrs, columns=train_ids)
    new_df = pd.DataFrame(X_new.T, index=selected_dmrs, columns=new_ids)
    combined = pd.concat([train_df, new_df], axis=1)

    condition = pd.Series(['Case' if lab == 1 else 'Control' for lab in y_train], index=train_ids)
    condition = pd.concat([condition, pd.Series(['New Sample'] * len(new_ids), index=new_ids)])
    condition = condition.reindex(combined.columns)

    condition_colors = {'Control': HEATMAP_CONTROL, 'Case': HEATMAP_CASE, 'New Sample': HEATMAP_NEW}
    col_colors = condition.map(condition_colors)
    col_colors.name = 'Condition'

    g = sns.clustermap(
        combined,
        cmap=DIVERGING_CMAP, vmin=0, vmax=1, center=0.5,
        col_colors=col_colors,
        method='complete', metric='euclidean',
        row_cluster=True, col_cluster=True,
        xticklabels=False, yticklabels=False,
        figsize=(11, 9),
        dendrogram_ratio=(0.12, 0.12),
        cbar_kws={'label': 'Methylation fraction', 'ticks': [0, 0.25, 0.5, 0.75, 1]},
    )
    g.fig.suptitle('Training Cohort + New Sample(s) -- ElasticNet-Selected DMR Methylation',
                    fontweight='bold', y=1.02)

    legend_handles = [Patch(facecolor=c, label=l) for l, c in condition_colors.items()]
    g.ax_heatmap.legend(handles=legend_handles, title='Condition', fontsize=8,
                         title_fontsize=9, loc='upper left', bbox_to_anchor=(1.25, 1.0))

    out_svg = f"{out_prefix}_heatmap.svg"
    g.savefig(out_svg, bbox_inches='tight')
    plt.close(g.fig)
    print(f"Saved: {out_svg}")


def main():
    parser = argparse.ArgumentParser(
        description='Score sample(s) with the ElasticNet+linear-SVM episignature '
                    'classifier and plot a PCA and heatmap against the training cohort.')
    parser.add_argument('--region-matrix', required=True,
                        help='CSV from 10_region_average_methylation.py (regions x samples), '
                             'covering the FULL DMR panel the model was trained on')
    parser.add_argument('--model-dir', required=True,
                        help='Directory with enet_svm_model.joblib, enet_svm_metadata.json')
    parser.add_argument('--training-data', default=None,
                        help='CSV of training-cohort DMR values for PCA/heatmap reference '
                             '(default: <model-dir>/enet_svm_training_data.csv)')
    parser.add_argument('--sample-col', default='sample_id')
    parser.add_argument('--label-col', default='group_bin')
    parser.add_argument('--out-prefix', required=True)
    args = parser.parse_args()

    model_dir = args.model_dir if args.model_dir.endswith('/') else args.model_dir + '/'
    training_data_path = args.training_data or (model_dir + 'enet_svm_training_data.csv')

    pipeline, meta = load_model(model_dir)
    train_ids, y_train, X_train_full, dmr_cols = load_training_data(
        training_data_path, meta, args.sample_col, args.label_col)
    new_ids, X_new_full = load_region_matrix(args.region_matrix, dmr_cols)

    scores, preds, indet = score_samples(pipeline, meta, X_new_full.copy())
    results = write_results(new_ids, scores, preds, indet, meta, args.out_prefix)

    scaler = pipeline.named_steps['scaler']
    support = pipeline.named_steps['select'].get_support()
    selected_dmrs = [d for d, keep in zip(dmr_cols, support) if keep]

    X_new_filled = X_new_full.copy()
    if np.isnan(X_new_filled).any():
        inds = np.where(np.isnan(X_new_filled))
        X_new_filled[inds] = np.take(scaler.mean_, inds[1])

    X_train_sc = scaler.transform(X_train_full)[:, support]
    X_new_sc = scaler.transform(X_new_filled)[:, support]

    plot_pca(train_ids, y_train, X_train_sc, new_ids, X_new_sc, results, args.out_prefix)
    plot_heatmap(train_ids, y_train, X_train_full[:, support], new_ids,
                 X_new_filled[:, support], selected_dmrs, args.out_prefix)

    print("\nDone.")
    print(results.to_string(index=False))


if __name__ == '__main__':
    main()
