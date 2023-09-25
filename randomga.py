# ============= begin random.scs ===============================

import math

# Contains random number generator and related utilities
# including advance_random, warmup_random, random, randomize,
# flip, rnd


# my additions - !!!

RND_FILE_NAME = 'RNDHAZIR.TXT';
RND_FROM_FILE = True;
global rndfile
global rndfileopened
rndfileopened = False


# Global variables - Don't use these names in other code
oldrand = [0.0] * 56  # Array of 55 random numbers
jrand = 0  # Current random
global randomseed

def advance_random():
    # Create the next batch of 55 random numbers
    global oldrand
    for j1 in range(1, 25):
        new_random = oldrand[j1] - oldrand[j1 + 31]
        if new_random < 0.0:
            new_random += 1.0
        oldrand[1] = new_random

    for j1 in range(25, 56):
        new_random = oldrand[j1] - oldrand[j1 - 24]
        if new_random < 0.0:
            new_random += 1.0
        oldrand[j1] = new_random

def warmup_random(random_seed):
    # Get random off and running
    global oldrand, jrand
    oldrand[55] = random_seed
    new_random = 1.0e-9
    prev_random = random_seed

    for j1 in range(1, 55):
        ii = 21 * j1 % 55
        oldrand[ii] = new_random
        new_random = prev_random - new_random
        if new_random < 0.0:
            new_random += 1.0
        prev_random = oldrand[ii]

    advance_random()
    advance_random()
    advance_random()
    jrand = 0

def random():
    # Fetch a single random number between 0.0 and 1.0 (Subtractive Method)
    # See Knuth, D. (1969), v. 2 for details
    global jrand
    global rndfile,rndfileopened
    tmpstr = ''

    if RND_FROM_FILE and rndfileopened:
        tmpstr=rndfile.readline().strip()
        if tmpstr =='':
            rndfile.close()
            rndfile=open(RND_FILE_NAME,"r")
            tmpstr=rndfile.readline().strip()
        return float(tmpstr)
    else:    
        jrand += 1
        if jrand > 55:
            jrand = 1
            advance_random()
        return oldrand[jrand]

def flip(probability):
    # Flip a biased coin, true if heads
    if probability == 1.0:
        return True
    else:
        return random() <= probability

def rnd(low, high):
    # Pick a random integer between low and high
    if low >= high:
        i = low
    else:
        i = math.floor(random() * (high - low + 1) + low)
        if i > high:
            i = high
    return i

def randomize():
    global randomseed
    global rndfile,rndfileopened

    if RND_FROM_FILE:
        rndfile = open(RND_FILE_NAME,"r")
        rndfileopened = True
    else:        # Get a seed number for random and start itup
        randomseed = 0.1  # You can change this seed value if desired
        warmup_random(randomseed)

# End of random.scs
