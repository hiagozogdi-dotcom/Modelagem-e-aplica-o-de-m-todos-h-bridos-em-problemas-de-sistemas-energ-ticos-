"""
train_RF.py

Floresta aleatória (RF) como modelo substituto da turbina a gás com gás de
síntese.

A RF é a média de B árvores de regressão,

    F(x) = (1/B) sum_b T_b(x),

cada uma ajustada a uma amostra bootstrap do treino e com um subconjunto
aleatório de atributos avaliado em cada divisão (max_features).

As configurações h do espaço H são avaliadas por validação cruzada de cinco
partições aplicada somente ao treino,

    MSE_CV(h) = (1/5) sum_k MSE_k(h),        h* = argmin_h MSE_CV(h),

com o normalizador min-max ajustado dentro de cada partição (Pipeline).
Divisão, normalização e etapas de execução em tcc_ml.py.

Uso
    python train_RF.py --cv                  seleção (busca aleatória em H)
    python train_RF.py --train --sel 5       treino final das 5 melhores

Autor  Hiago David Zogbi Silva Oliveira
"""

from __future__ import annotations

# CONFIGURAÇÃO

ALG = 'RF'
ALVOS = ['eta_t', 'Sgen_kW_K', 'EINOx_g_kg']

CV = False          # True = só validação cruzada | False = treino final de SEL
SEL = 5             # 'best', k (as k melhores), 'all' ou {'n_estimators': 300, ...}
MODO = 'random'     # 'random' (N_ITER sorteadas) | 'grid' (216 combinações)
N_ITER = 100
N_PROC = 4

H = {
    'n_estimators':      [100, 300, 500],
    'max_depth':         [10, 20, 30, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf':  [1, 2, 4],
    'max_features':      ['sqrt', 1.0],      # sqrt(nx) ou todos os atributos
}


import os
if N_PROC > 1:
    for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        os.environ[v] = '1'

import time
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestRegressor

import tcc_ml as tm


def _int(v):
    v = tm.num(v)
    return None if v is None else int(v)


CAST = {
    'n_estimators': _int,
    'max_depth': _int,                 # None = profundidade livre
    'min_samples_split': _int,
    'min_samples_leaf': _int,
    'max_features': tm.num,            # 'sqrt' ou fração (1.0 = todos)
}


def build(h, nj=1):
    return Pipeline([
        ('sx', MinMaxScaler()),
        ('m', RandomForestRegressor(**h, random_state=tm.SEED, n_jobs=nj)),
    ])


def fit(D, h, cv, rep, mod):
    tm.set_seed(tm.SEED)

    if cv:
        r, sc = tm.cv5(build(h), D['Xtr'], D['ytr_n'])
        sc.to_csv(rep / f"{tm.tag(h)}_r2_{r['cv_r2']:.4f}_mse_{r['cv_mse']:.4e}_cv_scores.csv", index=False)
        row = {**h, **r}
    else:
        model = build(h, nj=1 if N_PROC > 1 else -1)
        t0 = time.time()
        model.fit(D['Xtr'], D['ytr_n'].ravel())
        dt = time.time() - t0
        met = tm.aval(model.predict, D)
        name = (f"{ALG}_{tm.tag(h)}_"
                f"{met['train_r2']:.4f}_{met['train_mse']:.4e}_"
                f"{met['val_r2']:.4f}_{met['val_mse']:.4e}_"
                f"{met['test_r2']:.4f}_{met['test_mse']:.4e}_{dt:.2f}")
        p = mod / f'{name}.joblib'
        tm.bundle(p, alg=ALG, kind='sklearn', t=D['t'], F=D['F'], h=h,
                  pipe=model, sx=D['sx'], sy=D['sy'], lb=D['lb'], ub=D['ub'],
                  met=met, imp=dict(zip(D['F'], model['m'].feature_importances_)))
        row = {**h, **met, 'time_s': dt, 'file': tm.rel(p)}

    pd.DataFrame([row]).to_csv(rep / 'rows' / f'{tm.tag(h)}.csv', index=False)


if __name__ == '__main__':
    CV, SEL = tm.cli(CV, SEL)
    tm.run(ALG, H, CAST, fit, cv=CV, sel=SEL, modo=MODO, n_iter=N_ITER, n_proc=N_PROC, alvos=ALVOS)
