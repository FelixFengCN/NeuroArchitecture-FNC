#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCNtoolkit 0.30.post2: reproduce the original FNC normative-model workflow.

Inputs (original format):
  resp.txt            N_control x (ID + n_edges)   [whitespace, no header]
  cov.txt             N_control x (ID + age sex motion site)
  respdisorder.txt    N_patient x (ID + n_edges)
  covdisorder.txt     N_patient x (ID + age sex motion site)

The first column is a subject/observation ID by default (same as original index_col=0).
Use --no-id for files with no ID column. Rows must match exactly within each group.
All four input files must be whitespace-separated and have no header.

75-point normative age reference (6..80 inclusive):
  sex    = male code (default 1; change via --male-code if necessary)
  motion = average motion of healthy HCP-YA subjects (site ID 2 by default)
  site   = the HCP-YA site ID (default 2; change via --hcpya-site)

Models: one standard BLR per FNC edge, optimizer=powell,
        standardize=False, cubic B-spline p=3 nknots=5, age range [6,80].

Outputs:
  output/ROI_models_pos/edge####/Models/        serialized per-edge model
  output/ROI_models_pos/edge####/patient_results.npz
  output/ROI_models_pos/edge####/age75_results.npz
  output/fnc_normative_age75_and_patient_dscore.mat
  output/design/cov_bspline_{tr,di,age75}.txt
  output/design/design_information.json

The patient deviation score is PCNtoolkit's returned Z_di.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy.io import savemat


AGE_MIN_BASIS = 6.0
AGE_MAX_BASIS = 80.0
BSPLINE_DEGREE = 3
BSPLINE_NKNOTS = 5
AGE_GRID = np.arange(6, 81, dtype=np.float64)


def hash_file(path: Path):
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data-dir', type=Path,
                    default=Path('data/private_normative_inputs'))
    ap.add_argument('--healthy-fnc', type=Path, default=None, help='default: DATA_DIR/resp.txt')
    ap.add_argument('--healthy-cov', type=Path, default=None, help='default: DATA_DIR/cov.txt')
    ap.add_argument('--patient-fnc', type=Path, default=None, help='default: DATA_DIR/respdisorder.txt')
    ap.add_argument('--patient-cov', type=Path, default=None, help='default: DATA_DIR/covdisorder.txt')
    ap.add_argument('--hcpya-site', type=float, default=2.0,
                    help='Numeric site ID for HCP-YA (the second dataset; default 2.0)')
    ap.add_argument('--male-code', type=float, default=1.0,
                    help='Male encoding in cov.txt (default 1.0; use the coding of your input)')
    ap.add_argument('--out-dir', type=Path, default=None,
                    help='default: DATA_DIR/PCN_BLR_030_FNC_HCPYA_MALE_6_80')
    ap.add_argument('--no-id', action='store_true', help='Input matrices have no first ID column')
    ap.add_argument('--dry-run', action='store_true', help='Validate files and design only, no training')
    ap.add_argument('--limit-edges', type=int, default=None,
                    help='Test only first K FNC edges; save an explicitly partial MAT')
    ap.add_argument('--overwrite', action='store_true', help='Refit completed edges')
    ap.add_argument('--skip-in-sample-evaluation', action='store_true',
                    help='Fit without testresp/testcov; use if in-sample metrics are unnecessary')
    return ap.parse_args()


def read_numeric_matrix(path: Path, *, has_id: bool, label: str) -> tuple[np.ndarray, np.ndarray]:
    if not path.is_file():
        raise FileNotFoundError(f'{label}: {path}')
    # Original inputs are tab/space delimited; no header.
    data = pd.read_csv(path, sep=r'\s+', header=None, engine='c')
    if has_id:
        if data.shape[1] < 2:
            raise ValueError(f'{label} should have ID and numeric columns: {path}')
        ids = data.iloc[:, 0].astype(str).to_numpy()
        data = data.iloc[:, 1:]
    else:
        ids = np.arange(data.shape[0]).astype(str)
    if pd.Index(ids).has_duplicates:
        raise ValueError(f'{label} contains duplicated IDs; check input format.')
    try:
        arr = data.to_numpy(dtype=np.float64)
    except ValueError as exc:
        raise ValueError(f'{label} must contain only numeric values beyond ID: {path}') from exc
    if arr.ndim != 2 or not np.isfinite(arr).all():
        raise ValueError(f'{label} contains NaN, Inf, or is not a 2-D numeric matrix: {path}')
    return ids, arr


def check_same_rows(ids1, ids2, label):
    if len(ids1) != len(ids2) or not np.array_equal(ids1, ids2):
        raise ValueError(f'{label}: ID/order mismatch between response and covariates. '
                         'Stopping instead of silently misaligning subjects.')


def make_raw_covariates(arr, label):
    if arr.shape[1] != 4:
        raise ValueError(f'{label}: expected four numeric covariates: age, sex, motion, site; '
                         f'observed {arr.shape[1]} columns')
    df = pd.DataFrame(arr, columns=['age', 'sex', 'motion', 'site'])
    return df


def numeric_sites(site_series):
    return site_series.astype(float).to_numpy()


def design_from_cov(df, site_universe):

    cats = pd.Categorical(numeric_sites(df['site']), categories=site_universe)
    if np.any(pd.isna(cats)):
        raise ValueError('A site is missing from the fixed site universe')
    dummy = pd.get_dummies(cats, prefix='site', dtype=float)
    base = np.column_stack([
        df[['age', 'sex', 'motion']].to_numpy(dtype=float),
        dummy.to_numpy(dtype=float),
        np.ones(len(df), dtype=float)
    ])
    return base, [*['age', 'sex', 'motion'], *dummy.columns.tolist(), 'intercept']


def make_age_reference(healthy_cov, hcpya_site, male_code):
    hc_site_mask = np.isclose(healthy_cov['site'].to_numpy(dtype=float), hcpya_site)
    hc_hcpya = healthy_cov.loc[hc_site_mask]
    if hc_hcpya.empty:
        raise ValueError(f'HCP-YA site={hcpya_site} not present among healthy subjects. '
                         'Set --hcpya-site to the correct numeric ID.')
    motion_mean = float(hc_hcpya['motion'].mean())
    age_reference = pd.DataFrame({
        'age': AGE_GRID,
        'sex': np.repeat(float(male_code), len(AGE_GRID)),
        'motion': np.repeat(motion_mean, len(AGE_GRID)),
        'site': np.repeat(float(hcpya_site), len(AGE_GRID)),
    })
    metadata = {
        'sex_code_male': float(male_code),
        'hcpya_site_id': float(hcpya_site),
        'hcpya_motion_mean': motion_mean,
        'hcpya_n_healthy': int(hc_hcpya.shape[0]),
    }
    return age_reference, metadata


def build_designs(healthy_cov, patient_cov, age_reference_cov, B):
    sites = np.unique(np.concatenate([
        numeric_sites(healthy_cov['site']),
        numeric_sites(patient_cov['site']),
        numeric_sites(age_reference_cov['site']),
    ]))
    X_h_base, labels = design_from_cov(healthy_cov, sites)
    X_p_base, labels_p = design_from_cov(patient_cov, sites)
    X_age_base, labels_age = design_from_cov(age_reference_cov, sites)
    assert labels == labels_p == labels_age

    phi_h = np.asarray([B(x) for x in X_h_base[:, 0]], dtype=np.float64)
    phi_p = np.asarray([B(x) for x in X_p_base[:, 0]], dtype=np.float64)
    phi_age = np.asarray([B(x) for x in AGE_GRID], dtype=np.float64)
    for name, base, phi in [('healthy', X_h_base, phi_h),
                            ('patients', X_p_base, phi_p),
                            ('age75', X_age_base, phi_age)]:
        if len(base) != len(phi) or phi.ndim != 2:
            raise ValueError(f'{name} B-spline expansion unexpected shape {phi.shape}')
    X_h = np.column_stack((X_h_base, phi_h))
    X_p = np.column_stack((X_p_base, phi_p))
    X_age = np.column_stack((X_age_base, phi_age))
    if not (X_h.shape[1] == X_p.shape[1] == X_age.shape[1]):
        raise AssertionError('The three design matrices have different column counts')
    for mat in [X_h, X_p, X_age]:
        if not np.isfinite(mat).all():
            raise ValueError('Non-finite value in expanded covariate design matrix')
    full_labels = labels + [f'bspline_{i+1:02d}' for i in range(phi_h.shape[1])]
    return X_h, X_p, X_age, sites, full_labels


def as_vec(val, n_expected, label):
    v = np.asarray(val, dtype=np.float64).reshape(-1)
    if len(v) != n_expected:
        raise ValueError(f'{label}: expected {n_expected} values, got {len(v)}')
    if not np.isfinite(v).all():
        raise ValueError(f'{label}: found NaN/Inf from PCNtoolkit')
    return v


def main():
    args = parse_args()
    data_dir = args.data_dir.resolve()
    healthy_fnc = (args.healthy_fnc or data_dir / 'resp.txt').resolve()
    healthy_cov_path = (args.healthy_cov or data_dir / 'cov.txt').resolve()
    patient_fnc = (args.patient_fnc or data_dir / 'respdisorder.txt').resolve()
    patient_cov_path = (args.patient_cov or data_dir / 'covdisorder.txt').resolve()
    out_dir = (args.out_dir or data_dir / 'PCN_BLR_030_FNC_HCPYA_MALE_6_80').resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    has_id = not args.no_id

    id_h, Y_h = read_numeric_matrix(healthy_fnc, has_id=has_id, label='Healthy FNC')
    id_h_cov, cov_h_raw = read_numeric_matrix(healthy_cov_path, has_id=has_id,
                                             label='Healthy covariates')
    id_p, Y_p = read_numeric_matrix(patient_fnc, has_id=has_id, label='Patient FNC')
    id_p_cov, cov_p_raw = read_numeric_matrix(patient_cov_path, has_id=has_id,
                                             label='Patient covariates')
    check_same_rows(id_h, id_h_cov, 'Healthy')
    check_same_rows(id_p, id_p_cov, 'Patient')
    if Y_h.shape[1] != Y_p.shape[1]:
        raise ValueError(f'Healthy and patient FNC differ in number/order of edges: '
                         f'{Y_h.shape[1]} vs {Y_p.shape[1]}')
    healthy_cov = make_raw_covariates(cov_h_raw, 'healthy covariates')
    patient_cov = make_raw_covariates(cov_p_raw, 'patient covariates')
    age_reference, ref_details = make_age_reference(
        healthy_cov, hcpya_site=args.hcpya_site, male_code=args.male_code)
    ref_mode = 'HCPYA_male_site_and_HCPYA_mean_motion'

    try:
        import importlib.metadata as im
        version = im.version('pcntoolkit')
    except Exception as exc:
        raise RuntimeError('Use the fgz environment containing PCNtoolkit 0.30.post2') from exc
    if version != '0.30.post2':
        raise RuntimeError(f'Expected pcntoolkit 0.30.post2 but found {version}; '
                           'the exact old API and outputs must be preserved.')
    from pcntoolkit.normative import estimate, predict
    from pcntoolkit.util.utils import create_bspline_basis
    B = create_bspline_basis(AGE_MIN_BASIS, AGE_MAX_BASIS,
                             p=BSPLINE_DEGREE, nknots=BSPLINE_NKNOTS)
    X_h, X_p, X_age, site_categories, design_labels = build_designs(
        healthy_cov, patient_cov, age_reference, B)

    design_dir = out_dir / 'design'
    design_dir.mkdir(exist_ok=True)
    np.savetxt(design_dir / 'cov_bspline_tr.txt', X_h, fmt='%.16g')
    np.savetxt(design_dir / 'cov_bspline_di.txt', X_p, fmt='%.16g')
    np.savetxt(design_dir / 'cov_bspline_age75.txt', X_age, fmt='%.16g')
    manifest_file = design_dir / 'design_information.json'
    design_manifest = {
            'pcntoolkit_version': version, 'regression': 'blr', 'optimizer': 'powell',
            'warp': None, 'standardize': False,
            'basis_age_range': [AGE_MIN_BASIS, AGE_MAX_BASIS],
            'bspline_degree': BSPLINE_DEGREE, 'bspline_nknots': BSPLINE_NKNOTS,
            'age_grid': AGE_GRID.tolist(),
            'reference_mode': ref_mode,
            'reference_explanation': 'Age 6-80, male, HCP-YA mean motion, HCP-YA site',
            'reference_values': ref_details,
            'covariate_columns_in_order': design_labels,
            'site_categories_in_order': [str(x) for x in site_categories],
            'n_controls': len(id_h), 'n_patients': len(id_p),
            'n_edges': Y_h.shape[1],
            'input_files': {
                'healthy_fnc': str(healthy_fnc), 'healthy_cov': str(healthy_cov_path),
                'patient_fnc': str(patient_fnc), 'patient_cov': str(patient_cov_path),
            },
            'input_sha256': {
                'healthy_fnc': hash_file(healthy_fnc),
                'healthy_cov': hash_file(healthy_cov_path),
                'patient_fnc': hash_file(patient_fnc),
                'patient_cov': hash_file(patient_cov_path),
            },
            'in_sample_evaluation': not args.skip_in_sample_evaluation,
    }
    if manifest_file.exists() and not args.overwrite:
        with manifest_file.open('r', encoding='utf8') as f:
            old_manifest = json.load(f)
        if old_manifest != design_manifest:
            raise RuntimeError('Output folder already contains a DIFFERENT model/input configuration. '
                               'Choose a new --out-dir or explicitly pass --overwrite. '
                               'Cached model results must never be mixed across configurations.')
    with manifest_file.open('w', encoding='utf8') as f:
        json.dump(design_manifest, f, indent=2, ensure_ascii=False)

    n_edges = Y_h.shape[1]
    target_edges = n_edges if args.limit_edges is None else min(n_edges, args.limit_edges)
    if target_edges < 1:
        raise ValueError('No edges to run')
    print(f'Controls={len(id_h)}, patients={len(id_p)}, edges={n_edges}, '
          f'design_features={X_h.shape[1]}, age points=75')
    print('Expected PCNtoolkit:', version)
    print('Age reference:', ref_mode, ref_details)
    print('Design matrix shapes:', X_h.shape, X_p.shape, X_age.shape)
    if args.dry_run:
        print('DRY RUN DONE. Validated input matrices; trained no models.')
        return

    root = out_dir / 'ROI_models_pos'
    root.mkdir(exist_ok=True)
    patient_z = np.empty((len(id_p), target_edges), dtype=np.float64)
    patient_yhat = np.empty_like(patient_z)
    patient_s2 = np.empty_like(patient_z)
    normative_age75 = np.empty((75, target_edges), dtype=np.float64)
    normative_age75_s2 = np.empty_like(normative_age75)

    cov_tr = str(design_dir / 'cov_bspline_tr.txt')
    cov_di = str(design_dir / 'cov_bspline_di.txt')
    cov_age = str(design_dir / 'cov_bspline_age75.txt')
    original_cwd = os.getcwd()
    t0 = time.time()

    for j in range(target_edges):
        edge_name = f'edge{j+1:04d}'
        roi_dir = root / edge_name
        roi_dir.mkdir(parents=True, exist_ok=True)
        result_di = roi_dir / 'patient_results.npz'
        result_nm = roi_dir / 'age75_results.npz'
        if result_di.exists() and result_nm.exists() and not args.overwrite:
            with np.load(result_di) as s:
                z = as_vec(s['Z_di'], len(id_p), edge_name + ' Z_di cached')
                yhat = as_vec(s['yhat_di'], len(id_p), edge_name + ' yhat_di cached')
                s2 = as_vec(s['s2_di'], len(id_p), edge_name + ' s2_di cached')
            with np.load(result_nm) as s:
                yhat_nm = as_vec(s['yhat_nm'], 75, edge_name + ' yhat_nm cached')
                s2_nm = as_vec(s['s2_nm'], 75, edge_name + ' s2_nm cached')
        else:

            resp_tr = roi_dir / 'resp_tr.txt'
            resp_di = roi_dir / 'resp_di.txt'
            np.savetxt(resp_tr, Y_h[:, j], fmt='%.16g')
            np.savetxt(resp_di, Y_p[:, j], fmt='%.16g')
            os.chdir(roi_dir)
            try:
                estimate_kw = {
                    'alg': 'blr',
                    'optimizer': 'powell',
                    'savemodel': True,
                    'saveoutput': False,
                    'standardize': False,
                }
                if not args.skip_in_sample_evaluation:
                    estimate_kw['testresp'] = str(resp_tr)
                    estimate_kw['testcov'] = cov_tr
                estimate(cov_tr, str(resp_tr), **estimate_kw)
                model_path = str(roi_dir / 'Models')
                if not Path(model_path).is_dir():
                    raise FileNotFoundError(f'Model directory was not created: {model_path}')
                yhat_di_, s2_di_, z_di_ = predict(
                    cov_di, alg='blr', respfile=str(resp_di), model_path=model_path)
                yhat_nm_, s2_nm_ = predict(
                    cov_age, alg='blr', respfile=None, model_path=model_path)
                z = as_vec(z_di_, len(id_p), edge_name + ' Z_di')
                yhat = as_vec(yhat_di_, len(id_p), edge_name + ' yhat_di')
                s2 = as_vec(s2_di_, len(id_p), edge_name + ' s2_di')
                yhat_nm = as_vec(yhat_nm_, 75, edge_name + ' yhat_nm')
                s2_nm = as_vec(s2_nm_, 75, edge_name + ' s2_nm')
                np.savez_compressed(result_di, Z_di=z, yhat_di=yhat, s2_di=s2)
                np.savez_compressed(result_nm, yhat_nm=yhat_nm, s2_nm=s2_nm)
            finally:
                os.chdir(original_cwd)
        patient_z[:, j] = z
        patient_yhat[:, j] = yhat
        patient_s2[:, j] = s2
        normative_age75[:, j] = yhat_nm
        normative_age75_s2[:, j] = s2_nm
        if j == 0 or (j + 1) % 25 == 0 or (j + 1) == target_edges:
            print(f'[{j+1}/{target_edges}] {edge_name} completed. '
                  f'Elapsed: {(time.time()-t0)/60:.1f} min', flush=True)

    outfile = out_dir / ('fnc_normative_age75_and_patient_dscore.mat'
                         if target_edges == n_edges else
                         f'PARTIAL_{target_edges}_edges_fnc_normative_age75_and_patient_dscore.mat')
    mat_payload = {
        'ages': AGE_GRID.reshape(-1, 1),
        'FNC_normative_75': normative_age75,
        'FNC_normative_75_s2': normative_age75_s2,
        'Dscore': patient_z,
        'Z_di': patient_z,
        'yhat_di': patient_yhat,
        's2_di': patient_s2,
        'patient_FNC': Y_p[:, :target_edges],
        'patient_ids': np.asarray(id_p, dtype=object).reshape(-1, 1),
        'healthy_ids': np.asarray(id_h, dtype=object).reshape(-1, 1),
        'edge_names': np.asarray([f'edge{i+1:04d}' for i in range(target_edges)],
                                 dtype=object).reshape(1, -1),
        'pcntoolkit_version': version,
        'reference_mode': ref_mode,
        'reference_sex_male_code': float(args.male_code),
        'reference_hcpya_site': float(args.hcpya_site),
        'reference_hcpya_motion_mean': ref_details['hcpya_motion_mean'],
        'reference_hcpya_healthy_n': ref_details['hcpya_n_healthy'],
        'n_control': np.array([[len(id_h)]]),
        'n_patient': np.array([[len(id_p)]]),
        'total_edges_in_input': np.array([[n_edges]]),
        'edges_in_this_output': np.array([[target_edges]])
    }
    savemat(outfile, mat_payload, do_compression=True)
    if target_edges == n_edges:
        savemat(out_dir / 'normative_FNC_age6_80.mat', {
            'ages': mat_payload['ages'],
            'FNC_normative_75': normative_age75,
            'FNC_normative_75_s2': normative_age75_s2,
            'edge_names': mat_payload['edge_names'],
            'reference_mode': ref_mode,
            'reference_hcpya_motion_mean': ref_details['hcpya_motion_mean'],
            'reference_sex_male_code': float(args.male_code),
            'reference_hcpya_site': float(args.hcpya_site),
        }, do_compression=True)
        savemat(out_dir / 'patient_dscore.mat', {
            'Dscore': patient_z,
            'Z_di': patient_z,
            'yhat_di': patient_yhat,
            's2_di': patient_s2,
            'patient_ids': mat_payload['patient_ids'],
            'edge_names': mat_payload['edge_names'],
        }, do_compression=True)
    print('\nDONE. Final MATLAB file:', outfile)
    print('FNC_normative_75:', normative_age75.shape,
          'Dscore:', patient_z.shape, 'ages:', AGE_GRID.shape)
    print('Note: MATLAB rows of Dscore retain original patient input order.')


if __name__ == '__main__':
    main()

