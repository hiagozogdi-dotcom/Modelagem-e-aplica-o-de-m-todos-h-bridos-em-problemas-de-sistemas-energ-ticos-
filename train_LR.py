"""
train_LR.py

Regressão linear, modelo de referência para a turbina a gás com gás de síntese.

Sem busca de hiperparâmetros. O modelo

    y_n = b0 + sum_j b_j x_n,j

é ajustado por mínimos quadrados ordinários nas entradas normalizadas (min-max
ajustado no treino) e avaliado por validação cruzada de cinco partições no
treino e pelos conjuntos independentes de validação e teste. Com GRAU > 1 as
entradas são expandidas em termos polinomiais, o que mantém o modelo linear nos
coeficientes.

Os coeficientes também são escritos na base física (GRAU = 1). Como
x_n = (x - x_min)/(x_max - x_min) e y = y_min + (y_max - y_min) y_n,

    b_j,fis = b_j (y_max - y_min) / (x_max,j - x_min,j),
    b_0,fis = y_min + (y_max - y_min) b_0 - sum_j b_j,fis x_min,j.

Uso
    python train_LR.py

Saídas em results/Reports/LR_<resposta>/ (summary.csv, coeficientes.csv) e
results/Models/LR_<resposta>/.

Autor  Hiago David Zogbi Silva Oliveira
"""

from __future__ import annotations

import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler, PolynomialFeatures
from sklearn.linear_model import LinearRegression

import tcc_ml as tm


# CONFIGURAÇÃO

ALG = 'LR'
ALVOS = ['eta_t', 'Sgen_kW_K', 'EINOx_g_kg']
GRAU = 1            # 1 = regressão linear | 2 ou 3 = regressão polinomial



def build():
    st = [('sx', MinMaxScaler())]
    if GRAU > 1:
        st.append(('pf', PolynomialFeatures(degree=GRAU, include_bias=False)))
    st.append(('m', LinearRegression()))
    return Pipeline(st)


def coefs(pipe, D):
    """Coeficientes na base normalizada e, com GRAU = 1, na base física."""
    m = pipe['m']
    if GRAU > 1:
        nm = pipe['pf'].get_feature_names_out(D['F'])
        return pd.DataFrame({'termo': ['b0', *nm], 'b_n': [m.intercept_, *m.coef_]})
    dy = D['sy'].data_range_[0]
    dx = D['sx'].data_range_
    bf = m.coef_ * dy / dx
    b0 = D['sy'].data_min_[0] + dy * m.intercept_ - np.sum(bf * D['sx'].data_min_)
    return pd.DataFrame({'termo': ['b0', *D['F']],
                         'b_n': [m.intercept_, *m.coef_],
                         'b_fis': [b0, *bf]})


def fit(D, rep, mod):
    tm.set_seed(tm.SEED)

    cv, sc = tm.cv5(build(), D['Xtr'], D['ytr_n'])
    sc.to_csv(rep / f'GRAU_{GRAU}_cv_scores.csv', index=False)

    pipe = build()
    t0 = time.time()
    pipe.fit(D['Xtr'], D['ytr_n'].ravel())
    dt = time.time() - t0
    met = tm.aval(pipe.predict, D)

    name = (f"{ALG}_{GRAU}_"
            f"{met['train_r2']:.4f}_{met['train_mse']:.4e}_"
            f"{met['val_r2']:.4f}_{met['val_mse']:.4e}_"
            f"{met['test_r2']:.4f}_{met['test_mse']:.4e}_"
            f"{cv['cv_r2']:.4f}_{cv['cv_mse']:.4e}_{dt:.2f}")
    p = mod / f'{name}.joblib'
    tm.bundle(p, alg=ALG, kind='sklearn', t=D['t'], F=D['F'], h={'grau': GRAU},
              pipe=pipe, sx=D['sx'], sy=D['sy'], lb=D['lb'], ub=D['ub'], met={**cv, **met})

    coefs(pipe, D).to_csv(rep / 'coeficientes.csv', index=False)
    return {'rank': 1, 'grau': GRAU, **cv, **met, 'time_s': dt, 'file': tm.rel(p)}


if __name__ == '__main__':
    d, Fx = tm.read_data()
    print(f'{len(d)} pontos, {len(Fx)} entradas: {", ".join(Fx)}')

    for t in ALVOS:
        rep, mod = tm.dirs(ALG, t)
        D = tm.split(d, Fx, t)
        s = pd.DataFrame([fit(D, rep, mod)])
        s.to_csv(rep / 'summary.csv', index=False)
        print(f'\n[{ALG}] {t}')
        print(s[['cv_mse', 'cv_r2', 'train_r2', 'val_r2', 'test_r2']].to_string(index=False))
        print(pd.read_csv(rep / 'coeficientes.csv').to_string(index=False))
