from Gauss import *


# adiciona a matriz identidade a matriz A
def adiciona_identidade(A, quant_rest):
    for i in range(quant_rest):
        for j in range(quant_rest):
            if i == j:
                A[i].append(1)
            else:
                A[i].append(0)


# apaga as variaveis artificiais da lista das nao basicas e apaga as colunas relacionadas as variaveis artificiais
def fimFase1(A, nVarArtificial, varNBasica):
    for i in range(nVarArtificial):
        varNBasica.remove(len(A[0]) - 1 - i)

    for i in range(len(A)):
        for j in range(nVarArtificial):
            A[i].pop(-1)


def forma_padrao(tipo_prob, funcao_obj, restricoes, operadores, b):
    fase1 = False
    A = []
    varBasicas = []
    varNBasicas = []
    nVarArtificial = 0

    if tipo_prob == 'max':  # transforma em min se for max
        funcao_obj = [(x * -1) for x in funcao_obj]

    for i in range(len(b)):  # verifica se tem restricao negativa e multiplica a linha por -19
        if b[i] < 0:
            b[i] *= -1
            for j in range(len(restricoes[0])):
                restricoes[i][j] *= -1

            if operadores[i] == '<=':
                operadores[i] = '>='
            elif operadores[i] == '>=':
                operadores[i] = '<='

    # criacao da matriz A
    for i in range(len(restricoes)):
        lin = []
        for j in range(len(restricoes[0])):
            lin.append(restricoes[i][j])
        A.append(lin)

    #  adiciona colunas das variavies de folga e excesso
    for i in range(len(operadores)):
        if operadores[i] == '<=':
            funcao_obj.append(0)  # adiciona o custo da variavel de folga a função objetivo
            for j in range(len(A)):
                if i == j:
                    A[j].append(1)
                else:
                    A[j].append(0)
        elif operadores[i] == '>=':
            funcao_obj.append(0)  # adiciona o custo da variavel de excesso a função objetivo
            for j in range(len(A)):
                if i == j:
                    A[j].append(-1)
                else:
                    A[j].append(0)

    #  verifica necessidade de fase 1
    if '>=' in operadores:
        fase1 = True

    #  seleciona as primeiras variaveis básicas de acordo com o número de variaveis
    for i in range(len(restricoes[0])):
        varNBasicas.append(i)

    #  seleciona as variaveis básicas e não basicas
    if fase1:
        adiciona_identidade(A, len(restricoes))
        nVarArtificial = len(restricoes)

        for i in range(len(restricoes)):
            varBasicas.append(len(A[0]) - len(restricoes) + i)
            funcao_obj.append(0)

        for i in range(len(varNBasicas), len(A[0]) - len(restricoes)):
            varNBasicas.append(i)
    else:
        for i in range(len(restricoes)):
            varBasicas.append(i + len(restricoes[0]))

    print_matriz(A, 'A')
    print('funcao objetivo: {}'.format(funcao_obj))
    print('vetor b: {}'.format(b))
    print('variaveis basicas: {}'.format(varBasicas))
    print('variaveis nao basicas: {}'.format(varNBasicas))
    print('Fase 1: {}'.format(fase1))
    print('Numero de variaveis artificiais {}'.format(nVarArtificial))

    return A, funcao_obj, b, varBasicas, varNBasicas, fase1, nVarArtificial


def atualiza_B(A, BS):
    B = []

    for i in range(len(A)):
        linha = []
        for j in range(len(BS)):
            linha.append(A[i][BS[j]])

        B.append(linha)

    return B


def passo1(B, b, varBasicas, n):
    XB = gauss(B, b)

    X = [0] * n
    for i in range(len(varBasicas)):
        X[varBasicas[i]] = XB[i]

    print('XB :', end='')
    for i in XB:
        print('{:.2f}'.format(i), end=' ')
    print()

    print('X:', end='')
    for i in X:
        print('{:.2f}'.format(i), end=' ')
    print()
    return XB, X


def passo2_1(B, varBasicas, funcao_obj):
    custos = []
    for i in varBasicas:
        custos.append(funcao_obj[i]) # pega os custos das variavies basicas na funcao objetivo
    lbda = gauss(transposta(B), custos) # realiza o gauss entre a matriz B transposta e os custos das variveis basicas

    print('lambda:', end='')
    for i in lbda:
        print('{:.2f}'.format(i), end=' ')
    print()
    # print('vetor lbda: {}'.format(lbda))
    return lbda


def passo2_2(funcao_obj, lbda, A, varNBasicas):
    vet_cnk = []
    for cn in range(len(varNBasicas)):
        cnk = funcao_obj[varNBasicas[cn]] - mult(lbda, coluna(A, varNBasicas[cn]))
        vet_cnk.append(cnk)

    print ('vet_cnk: ', end='')
    for i in vet_cnk:
        print('{:.2f}'.format(i), end=' ')
    print()
    return vet_cnk


# passo 2.3 e 3
def passo3(vet_cnk, passo1):
    fim = False
    menor = min(vet_cnk)
    k = vet_cnk.index(menor)
    if menor < 0:
        print('\033[35msolucao nao eh otima\033[m')
        print('entra: N{}'.format(k + 1))
    elif passo1 is True:
        fim = True
        print('\033[31mProblema infactível\033[m')
    else:
        print('\033[35msolucao otima encontrada\033[m')
        fim = True
    return fim, k


def passo4(B, A, varNBasicas, k):
    y = gauss(B, coluna(A, varNBasicas[k]))

    print('y: ', end='')
    for i in y:
        print('{:.2f}'.format(i), end=' ')
    print()

    # print('y: {}'.format(y))
    return y


def passo5(XB, y):
    if max(y) <= 0:
        print('Solução ótima ilimitada')
        return None, True

    aux = float('inf')
    sai = 0

    for i in range(len(XB)):
        if y[i] > 0 and XB[i] / y[i] < aux:
            aux = XB[i] / y[i]
            sai = i
    print('sai: B{}'.format(sai + 1))
    return sai, False


def passo5_fase1(XB, y):
    if max(y) <= 0:
        print('Problema Original Infactível')
        return None, True

    aux = float('inf')
    sai = 0

    for i in range(len(XB)):
        if y[i] > 0 and XB[i]/y[i] < aux:
            aux = XB[i]/y[i]
            sai = i

    print('sai: B{}'.format(sai + 1))
    return sai, False


def passo6(BS, N, entra, sai):
    aux = BS[sai]
    BS[sai] = N[entra]
    N[entra] = aux
    return BS, N


def passo6_fase1(BS, N, entra, sai, nVarArtificial):
    fim = True
    aux = BS[sai]
    BS[sai] = N[entra]
    N[entra] = aux

    totalVars = len(BS) + len(N)
    for i in range(len(BS)):
        for j in range(nVarArtificial):
            if BS[i] == totalVars - 1 - j:
                fim = False

    return BS, N, fim


def printInfo(varBasicas, varNBasicas, B, itr):
    print('\n')
    print('\033[34mIteração: {}\033[m'.format(itr))
    print_matriz(B, 'B')
    print('variaveis basicas: {}'.format(varBasicas))
    print('variaveis nao basicas: {}'.format(varNBasicas))
    print('\n')


def final(restricoes, funcao_obj, tipo_prob, X):
    print('\n\n')
    solucao_otima = 0
    for i in range(len(restricoes[0])):
        solucao_otima += funcao_obj[i] * X[i]

    for i in range(len(restricoes[0])):
        print('X{} = {}'.format(i + 1, X[i]))
    if tipo_prob == 'max':
        solucao_otima *= -1
    print('Solução Otima = {}'.format(solucao_otima))


def main():
    count = 0
    tipo_prob = 'max'
    funcao_obj = [2, 3, 1]
    restricoes = [[-1, -4, -2],
                  [3, 2, 0]]
    operadores = ['<=',
                '>=']

    b = [-8, 6]

    A, funcao_obj, b, varBasicas, varNBasicas, fase1, nVarArtificial = forma_padrao(tipo_prob, funcao_obj, restricoes, operadores, b)

    problema = False

    if fase1:
        itr = 1
        funcao_obj_fase1 = [0] * len(funcao_obj)
        for i in range(len(funcao_obj_fase1) - 1, len(funcao_obj_fase1) - 1 - nVarArtificial, -1): # atribui 1 para todos os custos das variaveis artificiais
            funcao_obj_fase1[i] = 1

        print('\n\n\033[33mFase I\033[m')
        print('funcao objetivo fase I: {}'.format(funcao_obj_fase1))

        while True:
            B = atualiza_B(A, varBasicas)
            printInfo(varBasicas, varNBasicas, B, itr)
            Xb, X = passo1(B, b, varBasicas, len(varBasicas) + len(varNBasicas))
            lbda = passo2_1(B, varBasicas, funcao_obj_fase1)
            vet_cnk = passo2_2(funcao_obj_fase1, lbda, A, varNBasicas)
            fim, k = passo3(vet_cnk, True)

            if fim:
                break

            y = passo4(B, A, varNBasicas, k)
            sai, problema = passo5_fase1(Xb, y)

            if problema:
                break

            varBasicas, varNBasicas, fim = passo6_fase1(varBasicas, varNBasicas, k, sai, nVarArtificial)

            if fim:
                break

            itr += 1

        fimFase1(A, nVarArtificial, varNBasicas)

    itr = 1
    X = []

    print('\n\n\033[33mFase II\033[m')
    while True:
        B = atualiza_B(A, varBasicas)
        printInfo(varBasicas, varNBasicas, B, itr)
        Xb, X = passo1(B, b, varBasicas, len(varBasicas) + len(varNBasicas))
        lbda = passo2_1(B, varBasicas, funcao_obj)
        vet_cnk = passo2_2(funcao_obj, lbda, A, varNBasicas)
        fim, k = passo3(vet_cnk, False)

        if fim:
            break

        y = passo4(B, A, varNBasicas, k)
        sai, problema = passo5(Xb, y)

        if problema:
            break

        varBasicas, varNBasicas = passo6(varBasicas, varNBasicas, k, sai)

        itr += 1

    if not problema:
        final(restricoes, funcao_obj, tipo_prob, X)


if __name__ == '__main__':
    main()
