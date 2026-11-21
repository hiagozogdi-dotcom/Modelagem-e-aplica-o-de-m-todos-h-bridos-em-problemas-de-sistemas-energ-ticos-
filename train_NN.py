"""
train_NN.py

Rede neural artificial (perceptron multicamadas) como modelo substituto da
turbina a gás com gás de síntese.

Arquitetura (NL = 1)   x (nx) -> Dense(c1, ReLU) -> Dense(1, linear)
Arquitetura (NL = 2)   x (nx) -> Dense(c1, ReLU) -> Dense(c2, ReLU) -> Dense(1, linear)
Perda                  MSE, otimizador Adam (AMSGrad), interrupção antecipada

Etapas (tcc_ml.run)

  --cv     seleção. Cada h é avaliada por validação cruzada de cinco partições
           no treino, estratificadas pela resposta. Em cada partição o min-max
           de x é ajustado nas partições de ajuste e a interrupção antecipada usa
           a fração VSPLIT delas. Nenhum modelo final é treinado.
  --train  as configurações em SEL são treinadas no treino, com interrupção
           antecipada monitorando o erro na validação, e avaliadas em treino,
           validação e teste.

FRAC < 1 usa a fração FRAC de cada conjunto (subamostra estratificada, tcc_ml.sub)
e grava as saídas em pastas próprias, por exemplo NN1_f50_<resposta>.

Métricas r2 e mse na escala normalizada de y, rmse em unidades físicas. Junto
com cada modelo é gravado <modelo>_scalers.csv (min-max do treino), lido pelo
otimizador do Mathematica com os pesos do .keras.

Uso
    python train_NN.py --cv
    python train_NN.py --train --sel best

Autor  Hiago David Zogbi Silva Oliveira
"""

from __future__ import annotations

# CONFIGURAÇÃO

NL = 1              # camadas ocultas (1 ou 2)
ALG = f'NN{NL}'
ALVOS = ['eta_t', 'Sgen_kW_K', 'EINOx_g_kg']

CV = False          # True = só validação cruzada | False = treino final de SEL
SEL = 'best'        # 'best', k (as k melhores), 'all' ou {'c1': 46, 'lr': 0.003, 'batch': 32}
MODO = 'random'     # 'random' (N_ITER sorteadas) | 'grid' (todas), usado com --cv
N_ITER = 50
N_PROC = 1          # processos simultâneos (1 = uma configuração por vez)

EPOCH = 400
PAT = 15            # paciência da interrupção antecipada
VSPLIT = 0.15       # fração interna para a interrupção antecipada na CV
FRAC = 1            # fração de cada conjunto usada; 1 = todos os dados
NB = 10             # faixas da resposta para estratificar as partições da CV

# Espaço de busca. Os valores abaixo são a configuração selecionada no TCC;
# para refazer a busca, informe listas com os valores candidatos.
H = {
    'c1':    [46],
    'lr':    [0.003],
    'batch': [32],
}
if NL == 2:
    H = {'c1': [46], 'c2': [23], 'lr': [0.003], 'batch': [32]}
CAST = {'c1': int, 'c2': int, 'lr': float, 'batch': int}
CAST = {k: CAST[k] for k in H}


import os
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '0'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
if N_PROC > 1:
    for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
              'TF_NUM_INTRAOP_THREADS', 'TF_NUM_INTEROP_THREADS'):
        os.environ[v] = '1'

import gc
import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd

import keras
from keras import backend as K
from keras.optimizers import Adam
from keras.metrics import R2Score
from keras.models import Sequential
from keras.layers import Dense, Input
from keras.callbacks import EarlyStopping

from sklearn import metrics
from sklearn.preprocessing import MinMaxScaler

import tcc_ml as tm


def tag(h):
    c = f"{h['c1']}_{h['c2']}" if 'c2' in h else f"{h['c1']}"
    return f"{c}_{h['lr']:.6f}_{h['batch']}"


def set_seed(s):
    tm.set_seed(s)
    keras.utils.set_random_seed(s)


def build(nx, h):
    m = Sequential([Input(shape=(nx,)), Dense(h['c1'], activation='relu')]
                   + ([Dense(h['c2'], activation='relu')] if 'c2' in h else [])
                   + [Dense(1, activation='linear')])
    m.compile(optimizer=Adam(learning_rate=h['lr'], beta_1=0.855, beta_2=0.999, epsilon=1e-15, amsgrad=True),
              loss='mean_squared_error', metrics=[R2Score()])
    return m


def stop():
    return EarlyStopping(monitor='val_loss', mode='min', patience=PAT, restore_best_weights=True)


def cross_val(D, h, seed, rep):
    """Validação cruzada de cinco partições, com o min-max ajustado só nas
    partições de ajuste de cada rodada."""
    X, y = D['Xtr'], D['ytr_n'].ravel()
    kf, s = tm.kfold(y, 5, seed, nb=NB, how='quantile')
    rows = []
    for k, (ia, ib) in enumerate(kf.split(X, s)):
        set_seed(seed + k)
        sx = MinMaxScaler().fit(X[ia])
        m = build(X.shape[1], h)
        m.fit(sx.transform(X[ia]), y[ia], batch_size=h['batch'], epochs=EPOCH,
              validation_split=VSPLIT, callbacks=[stop()], shuffle=True, verbose=0)
        pa = np.ravel(m.predict(sx.transform(X[ia]), verbose=0))
        pb = np.ravel(m.predict(sx.transform(X[ib]), verbose=0))
        rows.append({'fold': k + 1,
                     'train_r2': metrics.r2_score(y[ia], pa),
                     'train_neg_mean_squared_error': -metrics.mean_squared_error(y[ia], pa),
                     'test_r2': metrics.r2_score(y[ib], pb),
                     'test_neg_mean_squared_error': -metrics.mean_squared_error(y[ib], pb)})
        K.clear_session()
        del m
        gc.collect()
    sc = pd.DataFrame(rows)
    mse = -sc['test_neg_mean_squared_error']
    r = {'cv_mse': mse.mean(), 'cv_mse_std': mse.std(ddof=1),
         'cv_r2': sc['test_r2'].mean(), 'cv_train_r2': sc['train_r2'].mean()}
    sc.to_csv(rep / f"{tag(h)}_r2_{r['cv_r2']:.4f}_mse_{r['cv_mse']:.4e}_cv_scores.csv", index=False)
    return r


def fit(D, h, cv, rep, mod):
    seed = 10 + h['c1'] * 100000 + h['batch'] * 1000 + int(h['lr'] * 1e6) + h.get('c2', 0) * 10**7
    set_seed(seed)

    if cv:
        row = {**h, **cross_val(D, h, seed, rep)}
    else:
        m = build(D['Xtr'].shape[1], h)
        t0 = time.time()
        hist = m.fit(D['Xtr_n'], D['ytr_n'], batch_size=h['batch'], epochs=EPOCH,
                     validation_data=(D['Xva_n'], D['yva_n']),
                     callbacks=[stop()], shuffle=True, verbose=0).history
        dt = time.time() - t0
        met = tm.aval(lambda Xn: m.predict(Xn, verbose=0), D, norm=True)
        name = (f"Arc_{'AB'[NL - 1]}_{tag(h)}_"
                f"{met['train_r2']:.4f}_{met['train_mse']:.4e}_"
                f"{met['val_r2']:.4f}_{met['val_mse']:.4e}_"
                f"{met['test_r2']:.4f}_{met['test_mse']:.4e}_{dt:.2f}")
        m.save(mod / f'{name}.keras')
        np.save(mod / f'{name}.npy', hist)
        p = mod / f'{name}.joblib'
        tm.bundle(p, alg=ALG, kind='keras', t=D['t'], F=D['F'], h=h, mfile=f'{name}.keras',
                  sx=D['sx'], sy=D['sy'], lb=D['lb'], ub=D['ub'], met=met, hist=hist)
        tm.scalers_csv({'F': D['F'], 't': D['t'], 'sx': D['sx'], 'sy': D['sy']},
                       mod / f'{name}_scalers.csv')
        row = {**h, **met, 'epochs': len(hist['loss']), 'time_s': dt, 'file': tm.rel(p)}

    pd.DataFrame([row]).to_csv(rep / 'rows' / f'{tag(h)}.csv', index=False)
    K.clear_session()
    gc.collect()


if __name__ == '__main__':
    CV, SEL = tm.cli(CV, SEL)
    tm.run(ALG, H, CAST, fit, cv=CV, sel=SEL, modo=MODO, n_iter=N_ITER, n_proc=N_PROC,
           alvos=ALVOS, frac=FRAC, tg=tag)
