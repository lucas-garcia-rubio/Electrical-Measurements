#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 14:29:15 2019

@author: lucas

Projeto da disciplina de instrumentação eletro-eletrônica

Projeto de um filtro Butterworth

Função de transferência: |H(jw)| =          1
                                    ------------------
                                    sqrt(1+(w/w0)**2n)
"""

import math
import numpy as np

#%% Obtendo parâmetros da projeto

#dB1 = float(input('Ganho na banda passante (dB): '))
dB1 = 0.5
#f1 = float(input('Frequencia para o ganho da banda passante (kHz): '))
f1 = 25
f1 = f1*1000
#dB2 = float(input('Atenuação na banda de corte (dB): '))
dB2 = 12
#f2 = float(input('Frequencia para a banda de corte (kHz): '))
f2 = 300
f2 = f2*1000

#%% Encontrando valores da função de transferência

pi = math.pi

w1 = 2*pi*f1     # transformando para rad/s
w2 = 2*pi*f2

g1 = 10**(-dB1/20)    # encontrando os ganhos  |H(w1)| = g1
g2 = 10**(-dB2/20)    #                        |H(w2)| = g2

#%% Encontrando (w2/w1)**2n para substituir e encontrar a ordem 'n' do filtro 

aux = ((1-g2**2)*(g1**2))/((1-g1**2)*(g2**2))
n = (math.log10(aux))/(math.log10(w2/w1))
n = n/2
n = math.ceil(n) # arredonda para o próximo inteiro
#%% Encontrando w0

w0 = (((g1**2)*(w1**(2*n)))/(1-(g1**2)))**(1/(2*n))

#%% Encontrando os polos
# pk = w0 * exp(j*pi/2) * exp(j * (2k-1)/2n * pi), k = polo em questão, 1, 2, 3... n

e = math.e

p = []
for i in range(n):
    k = i+1
    p.append(w0 * (np.exp(complex(0,pi/2))) * (np.exp(complex(0,(2*k-1)*pi/(2*n)))))
    
#%% Encontrando os resistores

def Sallen_Key(n, p, w0):
    if n%2 != 0: #se a ordem for par, necessita um filtro RC para o polo que não possui conjugado
        R = 1e3 #adota-se algum valor para calcular os capacitores
        n = n-1
        i = 0
        while i != n:
            wN = (p[i]*p[n])**(1/2)
            Q = -1*wN/(p[i]+p[n])
            C2 = 1/(2*wN*Q*R)
            C1 = ((2*Q)**2) * C2
            print("Sallen-Key ", i+1, ":")
            print("R1 = R2 = ", R, " ohms")
            print("C1 = ", C1.real, " F")
            print("C2 = ", C2.real, " F")
            print("--------------------------")
            i = i+1
            n = n-1
        Rc = 1e3 #resistor do filtro RC para ordens ímpares
        C3 = 1/(Rc*w0)
        print("Filtro RC:")
        print("R = ", Rc, " ohms")
        print("C = ", C3, " F")
        print("--------------------------")
    else:
        R = 1e3
        t = int(n/2) #para ir somente até a metade do vetor, já que ele será calculado com seu conjugado
        n = n-1
        for i in range(t+1):
            wN = (p[i]*p[n])**(1/2)
            Q = -1*wN/(p[i]+p[n])
            C2 = 1/(2*wN*Q*R)
            C1 = ((2*Q)**2) * C2
            print("Sallen-Key ", i+1, ":")
            print("R1 = R2 = ", R, " ohms")
            print("C1 = ", C1.real, " F")
            print("C2 = ", C2.real, " F")
            print("--------------------------")
Sallen_Key(n,p,w0)