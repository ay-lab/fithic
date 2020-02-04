import gzip
import argparse
import sys
import math
import numpy as np
import scipy.sparse as sps
import time
import os

def parse_args(arguments):
    parser = argparse.ArgumentParser(description="Check help flag")
    parser.add_argument("-i", "--interactions", help="Path to the interactions file to generate bias values",required=True, type=str)
    parser.add_argument("-f", "--fragments", help="Path to the interactions file to generate bias values",required=True, type=str)
    parser.add_argument("-o", "--output", help="Full path to output the generated bias file to", required=True, type=str)
    parser.add_argument("-x", "--percentOfSparseToRemove", help="Percent of diagonal to remove", required=False, type=float, default=0.05)
    return parser.parse_args()

def loadfastfithicInteractions(interactionsFile, fragsFile):
    print("Creating sparse matrix...")
    startT = time.time()
    with gzip.open(fragsFile, 'rt') as frag:
        ctr = 0
        fragDic = {}
        revFrag = []
        for lines in frag:
            line = lines.rstrip().split()
            chrom = line[0]
            mid = int(line[2])
            if chrom not in fragDic:
                fragDic[chrom]={}
            fragDic[chrom][mid]=ctr
            revFrag.append((chrom,mid))
            ctr+=1
    x = []
    y = []
    z = []
    with gzip.open(interactionsFile, 'rt') as ints:
        for lines in ints:
            line = lines.rstrip().split()
            chrom1 = line[0]
            mid1 = int(line[1])
            chrom2 = line[2]
            mid2 = int(line[3])
            z.append(float(line[4]))
            x.append(fragDic[chrom1][mid1])
            y.append(fragDic[chrom2][mid2])
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    sparseMatrix = sps.coo_matrix((z, (x,y)), shape=(ctr,ctr))
    rawMatrix = sparseMatrix + sparseMatrix.T
    endT = time.time()
    print("Sparse matrix creation took %s seconds" % (endT-startT))
    return rawMatrix, revFrag

def returnBias(rawMatrix, perc):
    R = rawMatrix.sum()
    mtxAndRemoved = removeZeroDiagonalCSR(rawMatrix, perc)
    print("Sparse rows removed")
    initialSize = rawMatrix.shape
    print("Initial matrix size: %s rows and %s columns" % (initialSize[0], initialSize[1]))
    rawMatrix = mtxAndRemoved[0]
    removed = mtxAndRemoved[1]
    newSize = rawMatrix.shape
    print("New matrix size: %s rows and %s columns" % (newSize[0], newSize[1]))

    print("Normalizing with KR Algorithm")
    result = knightRuizAlg(rawMatrix)
    colVec = result[0]
    #x = sps.diags(colVec.flatten(), 0, format='csr')

    bias = computeBiasVector(colVec)
    biasWZeros = addZeroBiases(removed, bias)
    return biasWZeros

def removeZeroDiagonalCSR(mtx, perc):
    iteration = 0
    toRemove = []
    ctr = 0
    rowSums = mtx.sum(axis=0)
    rowSums = list(np.array(rowSums).reshape(-1,))
    rowSums = list(enumerate(rowSums))
    rowSums.sort(key=lambda tup: tup[1])
    size = len(rowSums)
    rem = int(perc * size)
    print("Removing %s percent of most sparse bins" % (perc))
    print("... corresponds to %s total rows" % (rem))
    valToRemove = rowSums[rem][1]
    #print(valToRemove)
    print("... corresponds to all bins with less than or equal to %s total interactions" % valToRemove)
    for value in rowSums:
        if value[1] <= valToRemove:
            toRemove.append(value[0])
    list(set(toRemove))
    toRemove.sort()
    mtx = dropcols_coo(mtx, toRemove)
    for num in toRemove:
        if iteration != 0:
            num -= iteration
        removeRowCSR(mtx,num)
        iteration +=1
    return [mtx, toRemove]

def computeBiasVector(x):
    one = np.ones((x.shape[0],1))
    x = one/x
    sums = np.sum(x)
    avg = (1.0*sums)/x.shape[0]
    bias = np.divide(x,avg)
    return bias

def addZeroBiases(lst, vctr):
    for values in lst:
        vctr = np.insert(vctr,values,-1,axis=0)
    return vctr

def dropcols_coo(M, idx_to_drop):
    idx_to_drop = np.unique(idx_to_drop)
    C = M.tocoo()
    keep = ~np.in1d(C.col, idx_to_drop)
    C.data, C.row, C.col = C.data[keep], C.row[keep], C.col[keep]
    C.col -= idx_to_drop.searchsorted(C.col) # decrement column indices
    C._shape = (C.shape[0], C.shape[1] - len(idx_to_drop))
    return C.tocsr()

def removeRowCSR(mat, i):
    if not isinstance(mat, sps.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    n = mat.indptr[i+1] - mat.indptr[i]
    if n > 0:
        mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        mat.data = mat.data[:-n]
        mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        mat.indices = mat.indices[:-n]
    mat.indptr[i:-1] = mat.indptr[i+1:]
    mat.indptr[i:] -= n
    mat.indptr = mat.indptr[:-1]
    mat._shape = (mat._shape[0]-1, mat._shape[1])

def knightRuizAlg(A, tol=1e-6, f1 = False):
    n = A.shape[0]
    e = np.ones((n,1), dtype = np.float64)
    res = []

    Delta = 3
    delta = 0.1
    x0 = np.copy(e)
    g = 0.9

    etamax = eta = 0.1
    stop_tol = tol*0.5
    x = np.copy(x0)

    rt = tol**2.0
    v = x * (A.dot(x))
    rk = 1.0 - v
    rho_km1 = ((rk.transpose()).dot(rk))[0,0]
    rho_km2 = rho_km1
    rout = rold = rho_km1

    MVP = 0 #we'll count matrix vector products
    i = 0 #outer iteration count

    if f1:
        print ("it        in. it      res\n"),

    while rout > rt: #outer iteration
        i += 1

        if i > 30:
            break

        k = 0
        y = np.copy(e)
        innertol = max(eta ** 2.0 * rout, rt)
        while rho_km1 > innertol: #inner iteration by CG
            k += 1
            if k == 1:
                Z = rk / v
                p = np.copy(Z)
                #rho_km1 = np.dot(rk.T, Z)
                rho_km1 = (rk.transpose()).dot(Z)
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p

            if k > 10:
                break

            #update search direction efficiently
            w = x * A.dot(x * p) + v * p
           # alpha = rho_km1 / np.dot(p.T, w)[0,0]
            alpha = rho_km1 / (((p.transpose()).dot(w))[0,0])
            ap = alpha * p
            #test distance to boundary of cone
            ynew = y + ap
            if np.amin(ynew) <= delta:
                if delta == 0:
                    break

                ind = np.where(ap < 0.0)[0]
                gamma = np.amin((delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            if np.amax(ynew) >= Delta:
                ind = np.where(ynew > Delta)[0]
                gamma = np.amin((Delta - y[ind]) / ap[ind])
                y += gamma * ap
                break

            y = np.copy(ynew)
            rk -= alpha * w
            rho_km2 = rho_km1
            Z = rk / v
            #rho_km1 = np.dot(rk.T, Z)[0,0]
            rho_km1 = ((rk.transpose()).dot(Z))[0,0]
        x *= y
        v = x * (A.dot(x))
        rk = 1.0 - v
        #rho_km1 = np.dot(rk.T, rk)[0,0]
        rho_km1 = ((rk.transpose()).dot(rk))[0,0]
        rout = rho_km1
        MVP += k + 1
        #update inner iteration stopping criterion
        rat = rout/rold
        rold = rout
        res_norm = rout ** 0.5
        eta_o = eta
        eta = g * rat
        if g * eta_o ** 2.0 > 0.1:
            eta = max(eta, g * eta_o ** 2.0)
        eta = max(min(eta, etamax), stop_tol / res_norm)
        if f1:
            print ("%03i %06i %03.3f %e %e \n") % \
                (i, k, res_norm, rt, rout),
            res.append(res_norm)
    if f1:
        print ("Matrix - vector products = %06i\n") % \
            (MVP),

    return [x,i,k]

def checkBias(biasvec):
    std = np.std(biasvec)
    mean = np.mean(biasvec)
    median = np.median(biasvec)
    if (mean < 0.5 or mean > 2): 
        print("WARNING... Bias vector has a mean outside of typical range (0.5, 2).")
        print("Consider running with a larger -x option if problems occur")
        print("Mean\t%s" % mean)
        print("Median\t%s" % median)
        print("Std. Dev.\t%s" % std)
    else:
        if (median<0.5 or median > 2):
            print("WARNING... Bias vector has a median outside of typical range (0.5, 2).")
            print("Consider running with a larger -x option if problems occur")
            print("Mean\t%s" % mean)
            print("Median\t%s" % median)
            print("Std. Dev.\t%s" % std)
    return

def outputBias(biasCol, revFrag, outputFilePath):
    bpath = outputFilePath
    with gzip.open(bpath,'wt') as biasFile:
        ctr = 0
        for values in np.nditer(biasCol):
            chrommidTup = revFrag[ctr]
            chrom = chrommidTup[0]
            mid = chrommidTup[1]
            biasFile.write("%s\t%s\t%s\n" % (chrom, mid, values))
            ctr += 1

def main():
    args = parse_args(sys.argv[3:])
    matrix,revFrag = loadfastfithicInteractions(args.interactions, args.fragments)
    bias = returnBias(matrix, args.percentOfSparseToRemove)
    checkBias(bias)
    outputBias(bias, revFrag, args.output)


if __name__=="__main__":
    main()
