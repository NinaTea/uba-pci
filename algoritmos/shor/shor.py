
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer
from math import gcd
from numpy.random import randint, seed
from fractions import Fraction


"""
El algoritmo de Shor tiene dos partes -> Cuántica y clásica
1- Cuántica: usando la QFT busca el periodo de una funcion
2- Clásica: usa el periodo para factorizar el número


Se busca el periodo de esta funcion modular

f(x) = a^x mod N

Sea

a, N > INT_0
a < N y son coprimos

r es el periodo y es el entero mas chico que cumple:

a^r mod N = 1

Factorizamos N = 15

Comentarios despues de ejecutarlo varias veces:
O sea que,what, falopa pura
"""
N = 15
# Seguimos el algoritmo de la teorica
# 1. Elegimos un número "a" random entre 1 y N-1

seed(1)  # para poder repetir resultados
a = randint(2, N)
print(f'Se eligio a = {a}')

# 2. Verificamos que a y N son coprimos
if gcd(a, N) != 1:
    raise ValueError(f'No son coprimos: gcd({a}, {N}) = {gcd(a, N)}')


# 3. Buscamos el periodo r de la funcion a^{r} mod N = 1
# la fase que queremos es s/r y s es un numero random entre 0 y r-1

def c_amod15(a, power):
    # no hardcodeen
    if a not in [2, 4, 7, 8, 11, 13]:
        raise ValueError("'a' tiene que estar entre 2,4,7,8,11 or 13")
    U = QuantumCircuit(4)
    for _iteration in range(power):
        if a in [2, 13]:
            U.swap(2, 3)
            U.swap(1, 2)
            U.swap(0, 1)
        if a in [7, 8]:
            U.swap(0, 1)
            U.swap(1, 2)
            U.swap(2, 3)
        if a in [4, 11]:
            U.swap(1, 3)
            U.swap(0, 2)
        if a in [7, 11, 13]:
            for q in range(4):
                U.x(q)
    U = U.to_gate()
    U.name = f"{a}^{power} mod 15"
    c_U = U.control()
    return c_U


def qft_dagger(n):
    """Armamos QFT inversa"""
    qc = QuantumCircuit(n)

    for qubit in range(n//2):
        qc.swap(qubit, n-qubit-1)
    for j in range(n):
        for m in range(j):
            qc.cp(-np.pi/float(2**(j-m)), m, j)
        qc.h(j)
    qc.name = "QFT†"
    return qc


def qpe_amod15(a: int) -> float:
    """Estimamos la fase para N = 15.
    in: a
    out: fase estimada"""

    N_COUNT = 8
    qc = QuantumCircuit(4+N_COUNT, N_COUNT)

    for q in range(N_COUNT):
        qc.h(q)     # Inicializamos los qubits en estado |+>

    qc.x(N_COUNT)  # Seteamos un registro auxiliar en estado |1>

    for q in range(N_COUNT):  # Aplicamos la operacion controlada U
        qc.append(c_amod15(a, 2**q),
                  [q] + [i+N_COUNT for i in range(4)])

    # Aplicamos la inversa de QFT
    qc.append(qft_dagger(N_COUNT), range(N_COUNT))
    qc.measure(range(N_COUNT), range(N_COUNT))

    # simulamos el circuito
    aer_sim = Aer.get_backend('aer_simulator')
    # `memory=True` tells the backend to save each measurement in a list
    job = aer_sim.run(transpile(qc, aer_sim), shots=1, memory=True)
    lecturas = job.result().get_memory()

    print("Leemos el registro: " + lecturas[0])
    fase = int(lecturas[0], 2)/(2**N_COUNT)
    print(f"fase correspondiente: {fase}")
    return fase


fase = qpe_amod15(a)  # fase = s/r
frac = Fraction(fase).limit_denominator(N)
s, r = frac.numerator, frac.denominator
print(r)


posibles_factores_primos = [gcd(a**(r//2)-1, N), gcd(a**(r//2)+1, N)]
print(posibles_factores_primos)

# El algoritmo se repite hasta que encontrás factores no triviales de 15
# hay que meter un while
