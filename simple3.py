#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 11:54:20 2023

@author: user
"""

import randomga





# Constants
ncoins = 1000000  # number of coin flips
probability = 0.5  # probability of heads turning up

# Variables
heads_or_tails = [0, 0]  # heads/tails count
j=int()
toss=bool()

# Main program
randomga.randomize()  # Seed and warm up random number generator

for _ in range(ncoins):  # Coin toss loop
    xflip = 1 if randomga.flip(probability) else 0
    xrnd  = randomga.rnd(0,100)
    xrandom = randomga.random()
    print(f"sira: {_:9} xflip={xflip:2}  xrnd(0,100)={xrnd:8}  xrandom(0.0,1.0)={xrandom:20.18f}")

"""
    if toss:
        heads_or_tails[0] += 1
    else:
        heads_or_tails[1] += 1

print(f'In {ncoins} coin tosses there were {heads_or_tails[0]} heads and {heads_or_tails[1]} tails')
"""