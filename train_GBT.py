"""
train_GBT.py

Gradient boosting tree (GBT) como modelo substituto da turbina a gás com gás de
síntese.

O GBT constrói o modelo de forma aditiva,

    F_m(x) = F_{m-1}(x) + nu f_m(x),        m = 1, ..., M,

em que cada árvore f_m é ajustada ao gradiente negativo da perda quadrática
(os resíduos do estágio anterior), nu é a taxa de aprendizagem e M o número de
árvores. Com fração de amostragem < 1 cada árvore usa um subconjunto aleatório
do treino (boosting estocástico).

Divisão, normalização e etapas de execução em tcc_ml.py (70/15/15, min-max
ajustado no treino, divisão refeita em código a cada execução).

Uso
    python train_GBT.py --cv                 seleção por validação cruzada
    python train_GBT.py --train --sel best   treino final da melhor configuração

Saídas
    results/Reports/GBT_<resposta>/CV/summary_cv.csv   ranking por MSE_CV
    results/Reports/GBT_<resposta>/summary.csv         modelos treinados
    results/Models/GBT_<resposta>/*.joblib             modelo + normalizadores
    results/Models/GBT_<resposta>/*.npy                MSE x número de árvores

Autor  Hiago David Zogbi Silva Oliveira
"""

from __future__ import annotations

# CONFIGURAÇÃO


ALG = 'GBT'
ALVOS = ['eta_t', 'Sgen_kW_K', 'EINOx_g_kg']

CV = False          # True = só validação cruzada | False = treino final de SEL
SEL = 'best'        # 'best', k (as k melhores), 'all' ou {'n_estimators': 400, ...}
MODO = 'grid'       # 'grid' (todas as combinações) | 'random' (N_ITER sorteadas)
N_ITER = 100
N_PROC = 4          # processos simultâneos (1 = uma configuração por vez)

H = {
    'n_estimators':     [100, 200, 400],
    'learning_rate':    [0.01, 0.05, 0.10],
    'max_depth':        [2, 3, 4],
    'min_samples_leaf': [1, 3, 5],
    'subsample':        [0.7, 1.0],
}
CAST = {'n_estimators': int, 'learning_rate': float, 'max_depth': int,
        'min_samples_leaf': int, 'subsample': float}


import os
if N_PROC > 1:      # um núcleo por processo, o paralelismo vem dos processos
    for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        os.environ[v] = '1'

import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import GradientBoostingRegressor

import tcc_ml as tm


def build(h):
    return Pipeline([
        ('sx', MinMaxScaler()),
        ('m', GradientBoostingRegressor(loss='squared_error', **h, random_state=tm.SEED)),
    ])


def staged(model, D):
    """MSE normalizado em função do número de árvores (treino e validação)."""
    out = {}
    for k, nm in (('tr', 'loss'), ('va', 'val_loss')):
        xn = model['sx'].transform(D[f'X{k}'])
        y = D[f'y{k}_n'].ravel()
        out[nm] = np.array([np.mean((y - p) ** 2) for p in model['m'].staged_predict(xn)])
    return out


def fit(D, h, cv, rep, mod):
    """Avalia a configuração h (CV) ou treina e salva o modelo (treino final)."""
    tm.set_seed(tm.SEED)

    if cv:
        r, sc = tm.cv5(build(h), D['Xtr'], D['ytr_n'])
        sc.to_csv(rep / f"{tm.tag(h)}_r2_{r['cv_r2']:.4f}_mse_{r['cv_mse']:.4e}_cv_scores.csv", index=False)
        row = {**h, **r}
    else:
        model = build(h)
        t0 = time.time()
        model.fit(D['Xtr'], D['ytr_n'].ravel())
        dt = time.time() - t0
        met = tm.aval(model.predict, D)
        hist = staged(model, D)
        name = (f"{ALG}_{tm.tag(h)}_"
                f"{met['train_r2']:.4f}_{met['train_mse']:.4e}_"
                f"{met['val_r2']:.4f}_{met['val_mse']:.4e}_"
                f"{met['test_r2']:.4f}_{met['test_mse']:.4e}_{dt:.2f}")
        p = mod / f'{name}.joblib'
        tm.bundle(p, alg=ALG, kind='sklearn', t=D['t'], F=D['F'], h=h,
                  pipe=model, sx=D['sx'], sy=D['sy'], lb=D['lb'], ub=D['ub'],
                  met=met, hist=hist, imp=dict(zip(D['F'], model['m'].feature_importances_)))
        np.save(mod / f'{name}.npy', hist)
        row = {**h, **met, 'time_s': dt, 'file': tm.rel(p)}

    pd.DataFrame([row]).to_csv(rep / 'rows' / f'{tm.tag(h)}.csv', index=False)


if __name__ == '__main__':
    CV, SEL = tm.cli(CV, SEL)
    tm.run(ALG, H, CAST, fit, cv=CV, sel=SEL, modo=MODO, n_iter=N_ITER, n_proc=N_PROC, alvos=ALVOS)
