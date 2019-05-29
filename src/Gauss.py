from copy import deepcopy


def troca(a, b, index):
    for i in range(index+1, len(b)):
        if a[i][index] != 0:
            auxa = a[i]
            auxb = b[i]
            a[i] = a[index]
            b[i] = b[index]
            a[index] = auxa
            b[index] = auxb


def zera(a, b, pivo):
    for i in range(pivo+1, len(a)):
        const = a[i][pivo]/a[pivo][pivo]

        b[i][1] = b[i][1] - (const * b[pivo][1])
        for j in range(len(a[0])):
            a[i][j] = a[i][j] - (const * a[pivo][j])


def solucao(a, b):
    x = [0] * len(b)

    for i in range(len(a)-1, -1, -1):
        soma = 0
        for j in range(len(a[0])):
            if i != j:
                soma += a[i][j] * x[j]

        x[i] = (b[i][1] - soma) / a[i][i]

    return x


def gauss(A, B):
    bAux = []
    a = deepcopy(A)
    b = deepcopy(B)

    for i in range(len(b)):
        bAux.append([i, b[i]])

    for i in range(len(a)):
        if a[i][i] == 0:
            troca(a, bAux, i)
        zera(a, bAux, i)

    return solucao(a, bAux)


def transposta(matriz):
    trans = []
    for i in range(len(matriz[0])):
        aux = []
        for j in range(len(matriz)):
            aux.append(matriz[j][i])
        trans.append(aux)

    return trans


def mult(a, b):
    result = 0
    if len(a) != len(b):
        print('Erro de dimensão')
        exit()

    for i in range(len(a)):
        result += a[i] * b[i]

    return result


def coluna(A, index):
    col = []
    for i in range(len(A)):
        col.append(A[i][index])
    return col


def print_matriz(mat, nome):
    print('\nmatriz {}:'.format(nome))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            print(mat[i][j], end='\t')
        print()
    print('\n')
