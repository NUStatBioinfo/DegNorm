from sklearn.decomposition import TruncatedSVD
import numpy as np



def nmf(x, iter=100):

    svd = TruncatedSVD(n_components=1, algorithm='arpack')
    est = svd.fit_transform(x)

    lmbda = np.zeros(shape=x.shape)
    c = np.sqrt(iter)

    for _ in range(iter):
        res = x - est
        lmbda -= res * (1 / c)
        lmbda[lmbda < 0.] = 0.
        est = np.abs(svd.fit_transform(x + lmbda))

    return {'estimate': est
            , 'lambda': lmbda
            , 'residuals': res
            , 'svd': svd}