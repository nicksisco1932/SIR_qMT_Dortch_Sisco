import csv
import multiprocessing as mp

import numpy as np
from lmfit import minimize, Parameters
import julia

# TODO: is the data being indexed properly in the final data.
# TODO: send this to Richard.
# TODO: Clean up code
# TODO: Could do  this in a generator statement

X1 = []
Y1 = []
guesses = []
SMV = []
sigv = []

filename = 'xdata.txt'

with open('xdata.txt', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for col in plots:
        X1.append(col)
with open('ydata.txt', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for col in plots:
        Y1.append(col)
with open('guesses.txt', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for col in plots:
        guesses.append(col)
with open('estimated_saturated_fraction.txt', 'r') as csvfile:
# with open('estimated_saturated_fraction_800Lines.txt', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        SMV.append(row)
with open('measured_signals.txt', 'r') as csvfile:
# with open('measured_signals_800Lines.txt', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for col in plots:
        sigv.append(col)

xdata = X1[0]
ydata = Y1[0]  # named wrong

LB = np.array([0, 0.3, -1.05, 0])
UB = np.array([1, 3, 0, 10])


def fit_stuff(ti, td, M, p, sm):
    # % Get
    # pmf and Mfinf
    temp = p[0]
    a = np.asarray([1, float(temp[0])], dtype='double')
    M0 = a * float(temp[3])

    R1 = np.asarray([temp[1], temp[1]], dtype='double')
    S = np.asarray([temp[2], float(sm[0])], dtype='double')
    x = np.asarray(ti, dtype='double')
    delays = np.asarray(td, dtype='double')
    data = np.asarray(M, dtype='double')
    kmf = 14.5

    pmf = M0[1] / M0[0]
    Mfinf = M0[0]

    # % Get
    # R1f / R1m and Sf / Sm
    R1f = R1[0]
    R1m = R1[1]
    Sf = S[0]
    Sm = S[1]

    # % Calculate
    # R1 + / - in Eq.
    # 4
    R1diff = np.sqrt(np.square((R1f - R1m + (pmf - 1) * kmf)) + 4 * pmf * np.square(kmf))
    R1plus = (R1f + R1m + (1 + pmf) * kmf + R1diff) / 2
    R1minus = R1plus - R1diff

    # % Component
    # amplitude
    # terms
    # for td terms(Eq. 5)
    bftdplus = -(R1f - R1minus) / R1diff
    bftdminus = (R1f - R1plus) / R1diff
    bmtdplus = -(R1m - R1minus) / R1diff
    bmtdminus = (R1m - R1plus) / R1diff

    # % Signal
    # recovery
    # during
    # td(Eq.
    # 5)
    Mftd = bftdplus * np.exp(-R1plus * delays) + bftdminus * np.exp(-R1minus * delays) + 1
    Mmtd = bmtdplus * np.exp(-R1plus * delays) + bmtdminus * np.exp(-R1minus * delays) + 1

    # % Component
    # amplitude
    # terms
    # for ti terms(Eq. 5)
    bfplus = ((Sf * Mftd - 1) * (R1f - R1minus) + (Sf * Mftd - Sm * Mmtd) * pmf * kmf) / R1diff
    bfminus = -((Sf * Mftd - 1) * (R1f - R1plus) + (Sf * Mftd - Sm * Mmtd) * pmf * kmf) / R1diff
    # return (x[0]*t/(x[1]+t))-y
    model = (bfplus * np.exp(np.negative(R1plus) * x) + bfminus * np.exp(np.negative(R1minus) * x) + 1) * Mfinf

    params = Parameters()
    # params.add('bfplus', value=bfplus.all())
    # params.add('bfminus', value=bfminus.all())
    # params.add('R1plus', value=R1plus.all())
    # params.add('R1minus', value=R1minus.all())
    # params.add('Mfinf', value=Mfinf.all())
    params.add('Mm', value=model.all())

    return params


def residual(params, ti, data,i):
    # M = (bfplus. * exp(-R1plus * ti) + bfminus. * exp(-R1minus * ti) + 1) * Mfinf;
    x = np.asarray(ti, dtype='double')
    # bfplus = params['bfplus']
    # bfminus = params['bfminus']
    # R1plus = params['R1plus']
    # R1minus = params['R1minus']
    # Mfinf = params['Mfinf']
    params = fit_stuff(xdata, ydata, sigv[i], guesses, SMV[i])
    Mm = params['Mm']
    temp_m = np.abs(Mm)
    ddata = np.asarray(data, dtype='double')
    er = params['Mm'] - ddata
    params.add('P', value=er.all())
    # err = np.abs(Mm)-data
    return er


# M is the measured signal
# td are delay times
# ti is xdata
def work(i):
    # pool = mp.Pool(processes=12)
    # for i in range(0, len(sigv)):
    # print(sigv[i])
    #  i = 0
    print(i)
    #TODO: I think that this could be a generator instead so that it can be called in the minimize section instead of
    # here.
    #TODO: Use Julia

    params = fit_stuff(xdata, ydata, sigv[i], guesses, SMV[i])
    j = i
    result = minimize(residual, params, args=(xdata, sigv[i], j, ))
    final = result.residual
    return final

    # return result


def main():
    # pool = mp.Pool(processes=12)
    # for i in range(0, len(sigv)):
    # print(sigv[i])
    #  final_data = []
    pool = mp.Pool(processes=12)
    # for i in range(0,len(sigv)):
    # for i in range(0,800):
    # p = pool.apply_async(work, args=(i, ))
    # final_data.append(p.get())
    # print(p.get())
    # i = len(sigv)
    i = len(sigv)
    p = pool.map(work, range(0, i))
    print(np.asarray(p))
    pool.close()
    pool.join()
    # p = work(i-1)
    # final_data = np.array(p)
    # np.save('temp_data.txt',np.asarray(p))
    print('Processing has completed and your data has been saved as temp_Data.txt.gz')
    np.savetxt('temp_Data.txt.gz', np.asarray(p))
    # print(final_data)


if __name__ == '__main__':
    main()
