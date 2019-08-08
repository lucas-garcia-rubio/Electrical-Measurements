#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 16:57:46 2019

@author: felipe arisi
"""

#%% Area de import
import numpy as np
import math
import matplotlib.pyplot as plt

#%% Função que recebe um resistor ou capacitor e retorna o valor comercial do mesmo 
def resistor_comercial (R):
    contador = 0
    while(R>9.1):
        R = R/10
        contador += 1
    
    Rcomercial = np.array ([1.0,1.1,1.2,1.3,1.5,1.6,1.8,2.0,2.2,2.4,2.7,3.0,3.3,3.6,3.9,4.3,4.7,5.1,5.6,6.2,6.8,7.5,8.2,9.1])
    return (Rcomercial[(np.abs(Rcomercial - R)).argmin()])*(10**(contador))
    
def capacitor_comercial (C):
    contador = 0
    while(C<1):
        C = C*10
        contador += 1
    
    Ccomercial = np.array ([1.0,1.1,1.2,1.3,1.5,1.6,1.8,2.0,2.2,3.0,3.3,3.6,3.9,4.7,5.6,6.2,6.8,7.5,8.2,9.1])
    return (Ccomercial[(np.abs(Ccomercial - C)).argmin()])*(10**(-contador))

#%% Funçao que gera o diagrama de bode e retorna as raizes da funçao de transferencia
    
def bode_raizes(n,f1,epslon,omega1):
    frequencias = np.arange(1e3,8e6,1e3)
    razao = frequencias/f1
    if(n == 1):
        Hpb = 1./((1+(((epslon)**2)*(razao**2)))*(1/2))
    elif(n == 2):
        Hpb = 1./((1+(((epslon)**2)*(((2*(razao**2))-1)**2)))**(1/2))
    elif(n == 3):
        Hpb = 1./((1+(((epslon)**2)*((((4*(razao**3))-(3*(razao)))**2))))**(1/2))
    elif(n == 4):
        Hpb = 1/((1+(((epslon)**2)*(((8*(razao**4))-(8*(razao**2))+1)**2)))**(1/2))
    elif(n == 5):
        Hpb = 1./((1+(((epslon)**2)*((((16*(razao**5))-(20*(razao**3))+(5*(razao)))**2))))**(1/2))
    elif(n == 6):
        Hpb = 1./((1+(((epslon)**2)*(((32*(razao**6))-(48*(razao**4))+(18*(razao**2))-(1))**2)))**(1/2))
    else:
        print("valor acima do implementado")
        return(-1)
    
    pk = np.ones(n,dtype = 'complex') 
    for i in range(1,n+1):
        pk[i-1] = complex(-(omega1*(math.sin(((2*i-1)*math.pi)/(2*n))*math.sinh((1/n)*math.asinh(1/epslon)))),(omega1*(math.cos(((2*i-1)*math.pi)/(2*n))*math.cosh((1/n)*math.acosh(1/epslon)))))
        
    
    
    plt.subplot(211)
    plt.semilogx(frequencias,Hpb)
    plt.hlines(0.95,frequencias[0],frequencias[len(frequencias)-1],linestyles=':',color = 'red')
    plt.xlabel('Frequência(Hz)')
    plt.title('Diagrama de bode')
    plt.ylabel('Hpb(B)')
    plt.axvline(f1,ls=':',color = 'red')
    plt.axvline(f2,ls=':',color = 'green')
    plt.grid(True)
    
    plt.subplot(212)
    plt.semilogx(frequencias,20*np.log10(Hpb))
    plt.xlabel('Frequência(Hz)')
    plt.ylabel('Hpb(dB)')
    plt.axvline(f1,ls=':',color = 'red')
    plt.axvline(f2,ls=':',color = 'green')
    plt.grid(True)
    
    plt.show()
    
    return pk
    
#%% especificações do projeto
#Projetar um filtro passa baixa, resposta Butterworth, com ganho na banda passante
#de 0.5dB(f1 = 15kHz), atenuação de banda de corte de 16dB (f2 = 25kHz)
#
#

#%% Definição de variaveis 
f1 = 17.5e3
f2 = 30e3
A1 = 0.5 #ganho na banda passante
B1 = 12 #atenuação de banda de corte 
omega1 = 2*math.pi*f1
omega2 = 2*math.pi*f2

#%% Calculos dos ganhos
Gbc = f2/f1
epslon = round(math.sqrt((10**(A1/10))-1),3)
Hcorte = round(1./((10**B1)**(1/20)),3) 

n = math.ceil((math.acosh((1-((Hcorte)**2))/(epslon*Hcorte)))/math.acosh(Gbc))

#%% Circuito 

pk = bode_raizes(n,f1,epslon,omega1)
denominadores = np.ones((math.floor(n/2),3),dtype='complex') 
omegaN = np.ones(math.floor(n/2),dtype='complex')
Q = np.ones(math.floor(n/2),dtype='float')
R1 = np.ones(math.floor(n/2),dtype='float')
R2 = np.ones(math.floor(n/2),dtype='float')
C1 = np.ones(math.floor(n/2),dtype='float')
C2 = np.ones(math.floor(n/2),dtype='float')

R1[0] = 1e3
R2[0] = 1e3
#R1[1] = 1e3
#R2[1] = 1e3
for i in range (0,math.floor(n/2)):
    denominadores[i] = np.convolve([1,-pk[i]],[1,-pk[n-(1+i)]])
    omegaN[i] = (denominadores[i][2]) ** (1/2)
    Q[i] = (omegaN[i].real)/ -(pk[i]+pk[n-(1+i)])
    C1[i] = (Q[i]* (R1[i]+R2[i]))/(omegaN[i].real*R1[i]*R2[i])
    C2[i] = 1./(R1[i]*R2[i]*C1[i]*omegaN[i].real*omegaN[i].real)
    #C1[i] = capacitor_comercial(C1[i])
    #C2[i] = capacitor_comercial(C2[i])
    
    
if (n % 2 != 0):
    R3 = 1e3
    C3 = -1./(R3*pk[1].real)