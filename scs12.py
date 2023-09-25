#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import math
import sys
import copy


# ============= myutil.py ===================


CLASS_FILE                 ='CLASS100.DTA'
BB_ENABLE                  = False
REPORT_PERIOD              = 500
CONSOLE_REPORT_PERIOD      = 2000
PLOT_REPORT_PERIOD         = 500
GA_PERIOD                  = 100
REPORT_FILE                ='REPO12.OUT'
PLOT_FILE                  ='PLOT12.OUT'
MAX_ITER                   = 10000000
RAND_FROM_PASCAL           = False
RND_FROM_FILE              = False
RAND_SEED_SCS              = 0.3333
RAND_SEED_PASCAL           = 2      # kullanilmayacak


# fncfile = 'MULTIRND.DTA'
fncfile = CLASS_FILE
# fncfile = 'PERFECT.DTA'
# fncfile = 'CLASS100.DTA'
fnefile = 'ENVIRON.DTA'
fnrfile = 'REINF.DTA'
fntfile = 'TIME.DTA'
fngfile = 'GA.DTA'
fnrep = REPORT_FILE
fnpfile = PLOT_FILE

# ========= end of myinit.py ================

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')


global cfile, efile, rfile, tfile, gfile, rep, pfile

current_dir = os.path.dirname("/home/hadi/PycharmProjects/pythonProjects/scs/")
# current_dir = os.path.dirname("/Users/user/PycharmProjects/scs/")




# *********************************************************************

# declarations for scs

MAX_POSITION = 10
MAX_CLASS = 101
WILDCARD = -1


class ClassType:
    def __init__(self):
        self.c = [0] * MAX_POSITION
        self.a = 0
        self.strength = 0.0
        self.bid = 0.0
        self.ebid = 0.0
        self.matchflag = False
        self.specificity = 0


class ClassArray:
    def __init__(self):
        self.classes = [ClassType() for _ in range(MAX_CLASS)]
        # self.classes = [ClassType()]*MAX_CLASS


class ClassList:
    def __init__(self):
        self.clist = [-1] * MAX_CLASS  # 0 olursa o hücrede 0. class oldugu dusunulur
        self.nactive = 0  # clist[0] henüz dolu olmadigi icin


class PopType:
    def __init__(self):
        self.classifier = ClassArray()
        self.nclassifier = 0
        self.nposition = 0
        self.pgeneral = 0.0
        self.cbid = 0.0
        self.bidsigma = 0.0
        self.bidtax = 0.0
        self.lifetax = 0.0
        self.bid1 = 0.0
        self.bid2 = 0.0
        self.ebid1 = 0.0
        self.ebid2 = 0.0
        self.sumstrength = 0.0
        self.maxstrength = 0.0
        self.avgstrength = 0.0
        self.minstrength = 0.0


population = PopType()  # population of classifiers

matchlist = ClassList()  # who matched

envmessage = [0] * MAX_POSITION  # environmental message

# The 'rep' variable can be represented using Python's file handling mechanisms.
# For example, to open a file named "report.txt" for writing:
# rep = open("report.txt", "w")


# ============= begin random.scs ===============================

# import math !!!

# Contains random number generator and related utilities
# including advance_random, warmup_random, random, randomize,
# flip, rnd

RND_FILE_NAME = 'RNDHAZIR.TXT';
#RND_FROM_FILE = False;
global rndfile
global rndfileopened
rndfileopened = False

# Global variables - Don't use these names in other code
oldrand = [0.0] * 56  # Array of 55 random numbers
jrand = 0  # Current random


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
    global rndfile, rndfileopened
    tmpstr = ''

    if RND_FROM_FILE and rndfileopened:
        tmpstr = rndfile.readline().strip()
        if tmpstr == '':
            rndfile.close()
            rndfile = open(RND_FILE_NAME, "r")
            tmpstr = rndfile.readline().strip()[:6]
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
    global rndfile, rndfileopened

    if RND_FROM_FILE:
        rndfile = open(RND_FILE_NAME, "r")
        rndfileopened = True
    else:
        # Get a seed number for random and start it up
        randomseed = RAND_SEED_SCS  # You can change this seed value if desired
        warmup_random(randomseed)


# End of random.scs


# io.py: IO Routines - File opening routines

from enum import Enum


class QueryType(Enum):
    INTERACTIVE = 0
    BATCH = 1


def page(rep):
    rep.write('\x0C')  # Page break character


def open_input(query_flag, message, filename):
    if query_flag == QueryType.BATCH:
        input_file = open(filename, 'r')
    else:
        filename = input(f'Enter {message} filename: ')
        input_file = open(filename, 'r')
    return input_file


def open_output(query_flag, message, filename):
    if query_flag == QueryType.BATCH:
        output_file = open(filename, 'w')
    else:
        filename = input(f'Enter {message} filename: ')
        output_file = open(filename, 'w')
    return output_file


"""
# Example usage
if __name__ == "__main__":
    qflag = QueryType.INTERACTIVE
    fn = ""

    input_file = open_input(None, qflag, "input", fn)
    output_file = open_output(None, qflag, "output", fn)

    # Use the input and output files as needed
    page(output_file)

    input_file.close()
    output_file.close()
"""


# ================      end of io.scs    =================


# utility.py: utility procedures and functions

# import math !!!
# import random !!!

def poweri(x, i):
    """
    Compute x raised to the power of i.
    
    Parameters:
    x (float): The base number.
    i (int): The exponent.
    
    Returns:
    float: The result of x raised to the power of i.
    """
    powertemp = 1.0
    if i == 0:
        powertemp = 1.0
    elif i > 0:
        while i > 0:
            powertemp *= x
            i -= 1
    else:
        while i < 0:
            powertemp /= x
            i += 1
    return powertemp


# Global variables for randomnormaldeviate - watch for conflicting names
rndx2 = 0.0
rndcalcflag = False


def initrandomnormaldeviate():
    """
    Initialization routine for randomnormaldeviate.
    """
    global rndcalcflag
    rndcalcflag = True


def randomnormaldeviate():
    """
    Generate a random normal deviate after ACM algorithm 267 / Box-Muller Method.
    
    Returns:
    float: A random normal deviate.
    """
    global rndx2, rndcalcflag
    if rndcalcflag:
        rndx1 = math.sqrt(-2.0 * math.log(random()))
        t = 6.2831853072 * random()
        rndx2 = rndx1 * math.sin(t)
        rndcalcflag = False
        return rndx1 * math.cos(t)
    else:
        rndcalcflag = True
        return rndx2


def noise(mu, sigma):
    """
    Generate normal noise with the specified mean and standard deviation.
    
    Parameters:
    mu (float): The mean of the noise.
    sigma (float): The standard deviation of the noise.
    
    Returns:
    float: The generated noise.
    """
    return randomnormaldeviate() * sigma + mu


def rndreal(lo, hi):
    """
    Generate a random real number between the specified limits.
    
    Parameters:
    lo (float): The lower limit.
    hi (float): The upper limit.
    
    Returns:
    float: A random real number between lo and hi.
    """
    return random() * (hi - lo) + lo


def max(x, y):
    """
    Return the maximum of two values.
    
    Parameters:
    x (float): The first value.
    y (float): The second value.
    
    Returns:
        object: 
    float: The maximum value between x and y.
    """
    return x if x > y else y


def min(x, y):
    """
    Return the minimum of two values.
    
    Parameters:
    x (float): The first value.
    y (float): The second value.
    
    Returns:
    float: The minimum value between x and y.
    """
    return x if x < y else y


def avg(x, y):
    """
    Return the average of two values.
    
    Parameters:
    x (float): The first value.
    y (float): The second value.
    
    Returns:
    float: The average value between x and y.
    """
    return 0.5 * (x + y)


def halt():
    """
    Test for key press and query for halt flag.
    
    Returns:
    bool: True if the program should halt, False otherwise.
    """
    times = 100
    j = 0
    # while not msvcrt.kbhit() and j < times:
    while j < times:
        j += 1
    temp = (j < times)
    if temp:
        ch = input('Halt (y/n)? > ')
        temp = ch.lower() == 'y'
    return temp


"""
# Example usage
if __name__ == "__main__":
    initrandomnormaldeviate()
    print(randomnormaldeviate())
    print(noise(0, 1))
    print(rndreal(0, 1))
    print(max(5, 10))
    print(min(5, 10))
    print(avg(5, 10))
    print(halt())
"""

# from pynput import keyboard

"""
def halt2():
    
    #Test for key press and query for halt flag.

    #Returns:
    #bool: True if the program should halt, False otherwise.
    
    times = 100
    j = 0
    listener = keyboard.Listener()

    def on_key_press(key):
        nonlocal j
        j += 1
        if j >= times:
            listener.stop()

    listener.on_press = on_key_press
    listener.start()

    if j < times:
        ch = input('Halt (y/n)? > ')
        temp = ch.lower() == 'y'
    else:
        temp = True

    return temp

"""


# =================== end of utility.scs ===================

# environ.py: multiplexer environment


# Environment declarations
class ERecord:
    def __init__(self):
        self.laddress = 0
        self.ldata = 0
        self.lsignal = 0
        self.address = 0
        self.output = 0
        self.classifieroutput = 0
        # self.signal = []
        self.signal = [0] * MAX_POSITION


environrec = ERecord()


def generatesignal(environrec):
    """
    Generate a random signal.
    
    Parameters:
    environrec (ERecord): The environment record.
    """
    for j in range(0, environrec.lsignal):
        environrec.signal[j] = 1 if flip(0.5) else 0
    # writemessage(rep,environrec.signal,environrec.lsignal)
    # rep.write("\n")


def decode(mess, start, length):
    """
    Decode a substring as an unsigned binary integer.
    
    Parameters:
    mess (list): The input message (signal).
    start (int): The starting index. 0'dan başlamalı !!!
    length (int): The length of the substring. 
    
    Returns:
    int: The decoded integer value.
    """
    accum = 0
    powerof2 = 1
    for j in range(start, start + length):
        accum += powerof2 * mess[j]
        powerof2 *= 2
    return accum


def multiplexeroutput(environrec):
    """
    Calculate the correct multiplexer output.
    
    Parameters:
    environrec (ERecord): The environment record.
    """
    environrec.address = decode(environrec.signal, 0, environrec.laddress)
    # output için, signal içerisindeki kaçıncı bit'in decode edilen adrese ait output olduğu bulunacak.
    # environrec.output = environrec.signal[environrec.laddress + environrec.address + 1]
    # ilk index'in sıfır olmadı gereği +1'i kaldırdım. Böylece 0.30'lardan 0.47'lere çıktı!
    # bakalım devamındaki hata nerde... çıktıyı kontrol ederek buraya kadar olanden emin olalım
    environrec.output = environrec.signal[environrec.laddress + environrec.address]


def environment(environrec):
    """
    Coordinate multiplexer calculations.
    
    Parameters:
    environrec (ERecord): The environment record.
    """
    generatesignal(environrec)
    # writesignal(rep,environrec.signal,environrec.lsignal)
    multiplexeroutput(environrec)


def initenvironment(efile, environrec):
    """
    Initialize the multiplexer environment.
    
    Parameters:
    efile (file): The environment file.
    environrec (ERecord): The environment record.
    """
    environrec.laddress = int(efile.readline()[:3])
    environrec.ldata = int(poweri(2.0, environrec.laddress))
    environrec.lsignal = environrec.laddress + environrec.ldata
    environrec.address = 0
    environrec.output = 0
    environrec.classifieroutput = 0
    environrec.signal = [0] * (environrec.lsignal + 1)


def initrepenvironment(rep, environrec):
    """
    Write the initial environmental report.
    
    Parameters:
    rep (file): The report file.
    environrec (ERecord): The environment record.
    """
    rep.write("\n")
    rep.write("Environmental Parameters (Multiplexer)\n")
    rep.write("--------------------------------------\n")
    rep.write(f"Number of address lines   = {environrec.laddress:8}\n")
    rep.write(f"Number of data lines      = {environrec.ldata:8}\n")
    rep.write(f"Total number of lines     = {environrec.lsignal:8}\n")


def writesignal(rep, signal, lsignal):
    """
    Write a signal in bit-reverse order.
    
    Parameters:
    rep (file): The report file.
    signal (list): The signal to write.
    lsignal (int): The length of the signal.
    """

    # writemessage(rep,environrec.signal,environrec.lsignal)
    # rep.write("\n")

    for j in range(lsignal - 1, -1, -1):
        rep.write(str(signal[j]))
    rep.write("\n")


def reportenvironment(rep, environrec):
    """
    Write the current multiplexer info to the report file.
    
    Parameters:
    rep (file): The report file.
    environrec (ERecord): The environment record.
    """
    rep.write("\n")
    rep.write("Current Multiplexer Status\n")
    rep.write("--------------------------\n")
    rep.write("Signal                 = ")
    writesignal(rep, environrec.signal, environrec.lsignal)
    rep.write(f"Decoded address        = {environrec.address:8}\n")
    rep.write(f"Multiplexer output     = {environrec.output:11}\n")
    rep.write(f"Classifier output      = {environrec.classifieroutput:8}\n")


"""
# Example usage
if __name__ == "__main__":
    environrec = ERecord()
    initenvironment(None, environrec)
    print(environrec.laddress)
    print(environrec.ldata)
    print(environrec.lsignal)
    initrepenvironment(None, environrec)
    reportenvironment(None, environrec)
"""


# ====================== end of environ.scs =================

# detector.py: convert environmental states to env. message

# Detector data declarations
class DRecord:
    pass  # For this problem, no detector record is required


detectrec = DRecord()  # Dummy detector record


def detectors(environrec, detectrec, envmessage):
    """
    Convert environmental state to env. message.
    
    Parameters:
    environrec (ERecord): The environment record.
    detectrec (DRecord): The detector record (dummy in this case).
    envmessage (list): The environmental message (signal).
    """
    writesignal(rep, envmessage, environrec.lsignal)
    # envmessage = environrec.signal  # Copy signal message to env. message
    # envmessage = environrec.signal.copy()  # Copy signal message to env. message
    # envmessage.messages = [Bit(bit.value) for bit in environrec.signal.messages]
    envmessage = copy.deepcopy(environrec.signal)
    writesignal(rep, envmessage, environrec.lsignal)


def writemessage(rep, mess, lmessage):
    """
    Write a message in bit-reverse order.
    
    Parameters:
    rep (file): The report file.
    mess (list): The message to write.
    lmessage (int): The length of the message.
    """
    for j in range(lmessage - 1, -1, -1):
        # rep.write(f"{mess[j]:1}")
        rep.write(str(mess[j]))


def reportdetectors(rep, envmessage, nposition):
    """
    Write out environmental message.
    
    Parameters:
    rep (file): The report file.
    envmessage (list): The environmental message.
    nposition (int): The length of the message.
    """
    rep.write("\n")
    rep.write("Environmental message:    ")
    writemessage(rep, envmessage, nposition)
    rep.write("\n")


def initdetectors(efile, detectrec):
    """
    Dummy detector initialization.
    
    Parameters:
    efile (file): The environment file (not used in this case).
    detectrec (DRecord): The detector record (dummy).
    """
    pass


def initrepdetectors(rep, detectrec):
    """
    Dummy initial detectors report.
    
    Parameters:
    rep (file): The report file.
    detectrec (DRecord): The detector record (dummy).
    """
    pass


"""
# Example usage
if __name__ == "__main__":
    environrec = ERecord()
    envmessage = [0] * (environrec.lsignal + 1)  # Assuming environrec is defined
    detectors(environrec, detectrec, envmessage)
    print(envmessage)
    print(environrec.lsignal)
    rep = open("report.txt", "w")  # You can replace "report.txt" with your report file name
    reportdetectors(rep, envmessage, environrec.lsignal)
    rep.close()
"""


# =================== end of detector.scs ==========================


# perform.py: performance system - classifier matching

# Performance declarations - most are in declare.py

# cfile = None  # Classifier file

def randomchar(pgeneral):
    """
    Set position at random with specified generality probability.
    
    Parameters:
    pgeneral (float): The generality probability.
    
    Returns:
    int: Randomly chosen character (-1, 0, or 1).
    """
    if flip(pgeneral):
        return WILDCARD
    elif flip(0.5):
        return 1
    else:
        return 0


def readcondition(cfile, c, pgeneral, nposition):
    """
    Read a single condition.
    
    Parameters:
    cfile (file): The classifier file.
    c (list): The condition array.
    pgeneral (float): The generality probability.
    nposition (int): The number of positions.
    """
    for j in range(nposition - 1, -1, -1):
        ch = cfile.read(1)
        if ch == '0':
            c[j] = 0
        elif ch == '1':
            c[j] = 1
        elif ch == '#':
            c[j] = WILDCARD
        elif ch == 'R':
            c[j] = randomchar(pgeneral)
        # print(c[j])


def readclassifier(cfile, _class, pgeneral, nposition):
    """
    Read a single classifier.
    
    Parameters:
    cfile (file): The classifier file.
    _class (classtype): The classifier object.
    pgeneral (float): The generality probability.
    nposition (int): The number of positions.
    """
    # ch = cfile.read(1) #!!!

    readcondition(cfile, _class.c, pgeneral, nposition)  # Read condition
    cfile.read(1)  # read ":" & ignore
    _class.a = int(cfile.read(1))  # Read action (a single trit)
    _class.strength = float(cfile.readline()[:15])  # Read strength
    _class.bid = 0.0
    _class.ebid = 0.0
    _class.matchflag = False  # Initialization


def countspecificity(c, nposition):
    """
    Count condition specificity.
    
    Parameters:
    c (list): The condition array.
    nposition (int): The number of positions.
    
    Returns:
    int: The count of specific conditions.
    """

    temp = nposition - c[:nposition].count(WILDCARD)
    return temp

    """  hallettim gitti basbelasini
    temp = 0
    while nposition >= 0:
        if c[nposition-1] != WILDCARD:
            temp += 1
        nposition -= 1
    return temp
    """


def initclassifiers(cfile, population):
    """
    Initialize classifiers.
    
    Parameters:
    cfile (file): The classifier file.
    population (poptype): The population object.
    """
    nposition = int(cfile.readline()[:15])
    population.nposition = nposition
    population.nclassifier = int(cfile.readline()[:15])
    population.pgeneral = float(cfile.readline()[:15])
    population.cbid = float(cfile.readline()[:15])
    population.bidsigma = float(cfile.readline()[:15])
    population.bidtax = float(cfile.readline()[:15])
    population.lifetax = float(cfile.readline()[:15])
    population.bid1 = float(cfile.readline()[:15])
    population.bid2 = float(cfile.readline()[:15])
    population.ebid1 = float(cfile.readline()[:15])
    population.ebid2 = float(cfile.readline()[:15])
    for j in range(0, population.nclassifier):
        readclassifier(cfile, population.classifier.classes[j], population.pgeneral, nposition)
        population.classifier.classes[j].specificity = countspecificity(population.classifier.classes[j].c, nposition)


def initrepclassifiers(rep, population):
    """
    Initial report on population parameters.
    
    Parameters:
    rep (file): The report file.
    population (poptype): The population object.
    """
    rep.write("\n")
    rep.write("Population Parameters\n")
    rep.write("---------------------\n")
    rep.write(f"Number of classifiers  = {population.nclassifier:8}\n")
    rep.write(f"Number of positions    = {population.nposition:8}\n")
    rep.write(f"Bid coefficient        = {population.cbid:8.4f}\n")
    rep.write(f"Bid spread             = {population.bidsigma:8.4f}\n")
    rep.write(f"Bidding tax            = {population.bidtax:8.4f}\n")
    rep.write(f"Existence tax          = {population.lifetax:8.4f}\n")
    rep.write(f"Generality probability = {population.pgeneral:8.4f}\n")
    rep.write(f"Bid specificity base   = {population.bid1:8.4f}\n")
    rep.write(f"Bid specificity mult.  = {population.bid2:8.4f}\n")
    rep.write(f"ebid specificity base  = {population.ebid1:8.4f}\n")
    rep.write(f"ebid specificity mult. = {population.ebid2:8.4f}\n")


def writecondition(rep, c, nposition):
    """
    Convert internal condition format to external format and write to file/device.
    
    Parameters:
    rep (file): The report file.
    c (list): The condition array.
    nposition (int): The number of positions.
    """
    for j in range(nposition - 1, -1, -1):
        if c[j] == 1:
            rep.write("1")
        elif c[j] == 0:
            rep.write("0")
        elif c[j] == WILDCARD:
            rep.write("#")


def writeclassifier(rep, _class, number, nposition):
    """
    Write a single classifier.
    
    Parameters:
    rep (file): The report file.
    _class (classtype): The classifier object.
    number (int): The classifier number.
    nposition (int): The number of positions.
    """
    rep.write(f"{number + 1:5} {_class.strength:8.2f} {_class.bid:8.2f} {_class.ebid:8.2f}")
    if _class.matchflag:
        rep.write(" X ")
    else:
        rep.write("   ")
    writecondition(rep, _class.c, nposition)
    rep.write(f": [{_class.a}]\n")


def reportclassifiers(rep, population):
    """
    Generate classifiers report.
    
    Parameters:
    rep (file): The report file.
    population (poptype): The population object.
    """
    rep.write("\n")
    rep.write("No. Strength bid ebid M Classifier\n")
    rep.write("---------------------------------\n")
    rep.write("\n")
    for j in range(0, population.nclassifier):
        writeclassifier(rep, population.classifier.classes[j], j, population.nposition)


def match(c, m, nposition):
    """
    Match a single condition to a single message.
    
    Parameters:
    c (list): The condition array.
    m (list): The message (signal) array.
    nposition (int): The number of positions.
    
    Returns:
    bool: True if the condition matches the message, False otherwise.
    """

    matchtemp = True
    while matchtemp and nposition >= 0:
        matchtemp = (c[nposition - 1] == WILDCARD) or (c[nposition - 1] == m[nposition - 1])
        nposition -= 1
    return matchtemp


def matchclassifiers(population, emess, matchlist):
    """
    Match all classifiers against environmental message and create match list.
    
    Parameters:
    population (poptype): The population object.
    emess (message): The environmental message.
    matchlist (classlist): The match list object.
    
    Burada j: matchlist.clist dizisi içerisindeki indexler kaç classifier'in emess ile match olduğunu 
    gösterdiği için 0. index anlamsız olacak ve kullanılmayacak. 
    matchlist.clist dizi içeriği ise emess'e match olan son classifier'in numarası (0-nclassifier-1)
    olacak. 

    """
    matchlist.nactive = 0  # !!!
    for j in range(0, population.nclassifier):
        population.classifier.classes[j].matchflag = match(population.classifier.classes[j].c, emess,
                                                           population.nposition)
        try:
            if population.classifier.classes[j].matchflag:
                matchlist.nactive += 1
                matchlist.clist[matchlist.nactive - 1] = j  # !!!
                # matchlist.clist[matchlist.nactive] = j  #!!! Pascal icin orijinal satir (index uyarlandi)
        except:
            print("hata olustu...")


"""
# Example usage
if __name__ == "__main__":
    environrec = ERecord()  # Assuming environrec is defined
    emess = [0] * (environrec.lsignal + 1)  # Assuming emess is defined
    population = Poptype()  # Assuming population is defined
    matchlist = Classlist()  # Assuming matchlist is defined
    
    matchclassifiers(population, emess, matchlist)
    print(matchlist.nactive)
    rep = open("report.txt", "w")  # You can replace "report.txt" with your report file name
    reportclassifiers(rep, population)
    rep.close()
"""


# ================= end of perform.scs ==================================

class CRecord:
    def __init__(self):
        self.winner = -1  # class indexleri 0'dan basladigi icin
        self.oldwinner = -1
        self.bucketbrigadeflag = False


clearingrec = CRecord()


def initaoc(cfile, clearingrec):
    """
    Initialize clearinghouse record.
    
    Parameters:
    clearingrec (CRecord): The clearinghouse record object.
    """
    # print(cfile)
    ch = cfile.read(1)
    clearingrec.bucketbrigadeflag = ch.lower() == 'y'
    clearingrec.bucketbrigadeflag = BB_ENABLE
    
    clearingrec.winner = 0
    clearingrec.oldwinner = 0


def initrepaoc(rep, clearingrec):
    """
    Initial report of clearinghouse parameters.
    
    Parameters:
    rep (file): The report file object.
    clearingrec (CRecord): The clearinghouse record object.
    """
    rep.write("\n")
    rep.write("Apportionment of Credit Parameters\n")
    rep.write("----------------------------------\n")
    rep.write(f"Bucket brigade flag       = {clearingrec.bucketbrigadeflag}\n")


def auction(population, matchlist, oldwinner):
    """
    Auction among currently matched classifiers - return winner.
    
    Parameters:
    population (Poptype): The population object.
    matchlist (Classlist): The match list object.
    oldwinner (int): The old winner classifier index.
    
    Returns:
    int: The index of the winner classifier.

    !!! matchlist.nactive[0] ASLA kullanılmayacak-> fikir degisti!!!

    """
    bidmaximum = 0.0
    winner = oldwinner  # If no match, oldwinner wins again
    if matchlist.nactive > 0:
        for j in range(0, matchlist.nactive):
            k = matchlist.clist[j]
            # with population.classifier.classes[k] as _class:
            population.classifier.classes[k].bid = population.cbid * (
                        population.bid1 + population.bid2 * population.classifier.classes[k].specificity) * \
                                                   population.classifier.classes[k].strength
            population.classifier.classes[k].ebid = population.cbid * (
                        population.ebid1 + population.ebid2 * population.classifier.classes[k].specificity) * \
                                                    population.classifier.classes[k].strength + noise(0.0,
                                                                                                      population.bidsigma)

            if population.classifier.classes[k].ebid > bidmaximum:
                winner = k
                bidmaximum = population.classifier.classes[k].ebid

    return winner


def clearinghouse(population, clearingrec):
    """
    Distribute payment from recent winner to oldwinner.
    
    Parameters:
    population (Poptype): The population object.
    clearingrec (CRecord): The clearinghouse record object.
    """
    # with population.classifier.classes[clearingrec.winner] as winner_classifier:
    payment = population.classifier.classes[clearingrec.winner].bid
    population.classifier.classes[clearingrec.winner].strength -= payment

    if clearingrec.bucketbrigadeflag:
        population.classifier.classes[clearingrec.oldwinner].strength += payment


def taxcollector(population):
    """
    Collect existence and bidding taxes from population members.
    
    Parameters:
    population (Poptype): The population object.
    """
    bidtaxswitch = 0  # define first to remember

    if (population.lifetax != 0.0) or (population.bidtax != 0.0):
        for j in range(0, population.nclassifier):
            # with population.classifier.classes[j] as _class:

            if population.classifier.classes[j].matchflag:
                bidtaxswitch = 1.0
            else:
                bidtaxswitch = 0.0

            population.classifier.classes[j].strength += -population.lifetax * population.classifier.classes[
                j].strength - population.bidtax * bidtaxswitch * population.classifier.classes[j].strength


def reportaoc(rep, clearingrec):
    """
    Report who pays to whom.
    
    Parameters:
    rep (file): The report file object.
    clearingrec (CRecord): The clearinghouse record object.
    """
    rep.write("\n")
    rep.write(f"New winner ({clearingrec.winner + 1}) : Old winner ({clearingrec.oldwinner + 1})\n")


def aoc(population, matchlist, clearingrec):
    """
    Apportionment of credit coordinator.
    
    Parameters:
    population (Poptype): The population object.
    matchlist (Classlist): The match list object.
    clearingrec (CRecord): The clearinghouse record object.
    """
    # with clearingrec:
    clearingrec.winner = auction(population, matchlist, clearingrec.oldwinner)
    taxcollector(population)
    clearinghouse(population, clearingrec)


"""
# Example usage
if __name__ == "__main__":
    population = Poptype()  # Assuming population is defined
    matchlist = Classlist()  # Assuming matchlist is defined
    clearingrec = CRecord()  # Create an instance of CRecord
    
    aoc(population, matchlist, clearingrec)
    
    rep = open("report.txt", "w")  # Replace with your report file name
    reportaoc(rep, clearingrec)
    rep.close()
"""


# ===================== end of aoc.scs ===================================

def effector(population, clearingrec, environrec):
    """
    Set action in object as dictated by auction winner.
    
    Parameters:
    population (Poptype): The population object.
    clearingrec (CRecord): The clearinghouse record object.
    environrec (ERecord): The environmental record object.
    """

    # Set classifieroutput as the action of the winner classifier
    environrec.classifieroutput = population.classifier.classes[clearingrec.winner].a


"""
# Example usage
if __name__ == "__main__":
    population = Poptype()  # Assuming population is defined
    clearingrec = CRecord()  # Assuming clearingrec is defined
    environrec = ERecord()  # Assuming environrec is defined
    
    effector(population, clearingrec, environrec)
    
    # Now the classifieroutput field in population object holds the action of the winner classifier

"""


# ======================== end of effector.scs ===========================

# reinforc.py: reinforcement and criterion procedures

# rfile=None !!! #

# Reinforcement data declarations
class RRecord:
    def __init__(self):
        self.reward = 0.0
        self.rewardcount = 0.0
        self.totalcount = 0.0
        self.count50 = 0.0
        self.rewardcount50 = 0.0
        self.proportionreward = 0.0
        self.proportionreward50 = 0.0
        self.lastwinner = -1


reinforcementrec = RRecord()


# Initialize reinforcement parameters
def initreinforcement(rfile, reinforcementrec):
    reinforcementrec.reward = float(rfile.readline()[:5])
    reinforcementrec.rewardcount = 0.0
    reinforcementrec.rewardcount50 = 0.0
    reinforcementrec.totalcount = 0.0
    reinforcementrec.count50 = 0.0
    reinforcementrec.proportionreward = 0.0
    reinforcementrec.proportionreward50 = 0.0
    reinforcementrec.lastwinner = -1


# Initial reinforcement report
def initrepreinforcement(rep, reinforcementrec):
    rep.write('\n')
    rep.write('Reinforcement parameters\n')
    rep.write('------------------------\n')
    rep.write('Reinforcement reward = {:8.1f}\n'.format(reinforcementrec.reward))


# Return true if criterion is achieved
def criterion(rrec, environrec):
    tempflag = (environrec.output == environrec.classifieroutput)
    rrec.totalcount += 1
    rrec.count50 += 1

    # Increment reward counters
    if tempflag:
        rrec.rewardcount += 1
        rrec.rewardcount50 += 1

    # Calculate reward proportions: running & last 50
    rrec.proportionreward = rrec.rewardcount / rrec.totalcount
    if round(rrec.count50 - 50.0) == 0:
        rrec.proportionreward50 = rrec.rewardcount50 / 50.0
        rrec.rewardcount50 = 0.0
        rrec.count50 = 0.0  # Reset

    return tempflag


# Pay reward to appropriate individual
def payreward(population, rrec, clearingrec):
    # with population, rrec, clearingrec:
    #    with population.classifier.classes[winner]:
    population.classifier.classes[clearingrec.winner].strength += rrec.reward
    rrec.lastwinner = clearingrec.winner


# Report reinforcement award
def reportreinforcement(rep, rrec):
    rep.write('\n')
    rep.write('Reinforcement Report\n')
    rep.write('-------------------\n')
    rep.write('Proportion Correct (from start) = {:8.4f}\n'.format(rrec.proportionreward))
    rep.write('Proportion Correct (last fifty) = {:8.4f}\n'.format(rrec.proportionreward50))
    rep.write('Last winning classifier number  = {:8}\n'.format(rrec.lastwinner + 1))


# Make payment if criterion satisfied
def reinforcement(reinforcementrec, population, clearingrec, environrec):
    if criterion(reinforcementrec, environrec):
        payreward(population, reinforcementrec, clearingrec)


"""
# Example usage
if __name__ == "__main__":
    reinforcementrec = RRecord()  # Assuming reinforcementrec is defined
    population = Poptype()  # Assuming population is defined
    clearingrec = CRecord()  # Assuming clearingrec is defined
    environrec = ERecord()  # Assuming environrec is defined
    rfile = open("reinforcement_data.txt", "r")  # Open reinforcement data file

    init_reinforcement(rfile, reinforcementrec)
    
    # Rest of the code: perform reinforcement and report

    rfile.close()  # Close reinforcement data file

"""
# =================== end of reinforc.scs =============================

# timekeep.py: timekeeper routines

# Data declarations

# tfile=None

# Number of iterations per block
# global iterationsperblock, carryflag, dummyflag
iterationsperblock = 10000


# Timekeeper record type
class TRecord:
    """Timekeeper record type."""

    def __init__(self):
        self.initialiteration = 0
        self.initialblock = 0
        self.iteration = 0
        self.block = 0
        self.reportperiod = 0
        self.gaperiod = 0
        self.consolereportperiod = 0
        self.plotreportperiod = 0
        self.nextplotreport = 0
        self.nextconsolereport = 0
        self.nextreport = 0
        self.nextga = 0
        self.reportflag = False
        self.gaflag = False
        self.consolereportflag = False
        self.plotreportflag = False


timekeeprec = TRecord()


# Add time function
# def addtime(t, dt, carryflag):
def addtime(t, dt, carry):
    """Increment iterations counter and set carry flag if necessary."""
    tempadd = t + dt
    carry[0] = tempadd >= iterationsperblock
    if carry[0]:
        tempadd = tempadd % iterationsperblock
    return tempadd


# Initialize timekeeper
def inittimekeeper(tfile, timekeeprec):
    """Initialize timekeeper."""
    dummyflag = [False]
    timekeeprec.iteration = 0
    timekeeprec.block = 0
    timekeeprec.initialiteration = int(tfile.readline()[:15])
    timekeeprec.initialblock = int(tfile.readline()[:15])
    timekeeprec.reportperiod = int(tfile.readline()[:15])
    timekeeprec.reportperiod = REPORT_PERIOD
    timekeeprec.consolereportperiod = int(tfile.readline()[:15])
    timekeeprec.consolereportperiod = CONSOLE_REPORT_PERIOD
    timekeeprec.plotreportperiod = int(tfile.readline()[:15])
    timekeeprec.plotreportperiod = PLOT_REPORT_PERIOD
    timekeeprec.gaperiod = int(tfile.readline()[:15])
    timekeeprec.gaperiod = GA_PERIOD
    timekeeprec.iteration = timekeeprec.initialiteration
    timekeeprec.block = timekeeprec.initialblock
    timekeeprec.nextga = addtime(timekeeprec.iteration, timekeeprec.gaperiod, dummyflag)
    timekeeprec.nextreport = addtime(timekeeprec.iteration, timekeeprec.reportperiod, dummyflag)
    timekeeprec.nextconsolereport = addtime(timekeeprec.iteration, timekeeprec.consolereportperiod, dummyflag)
    timekeeprec.nextplotreport = addtime(timekeeprec.iteration, timekeeprec.plotreportperiod, dummyflag)


# Initialize timekeeper report
def initreptimekeeper(rep, timekeeprec):
    """Initial timekeeper report."""
    rep.write("\n")
    rep.write("Timekeeper Parameters\n")
    rep.write("---------------------\n")
    rep.write(f"Initial iteration         ={timekeeprec.initialiteration:8}\n")
    rep.write(f"Initial block             ={timekeeprec.initialblock:8}\n")
    rep.write(f"Report period             ={timekeeprec.reportperiod:8}\n")
    rep.write(f"Console report period     ={timekeeprec.consolereportperiod:8}\n")
    rep.write(f"Plot report period        ={timekeeprec.plotreportperiod:8}\n")
    rep.write(f"Genetic algorithm period  ={timekeeprec.gaperiod:8}\n")


# Timekeeper function
def timekeeper(timekeeprec):
    """Keep time and set flags for time-driven events."""

    carryflag = [False]
    dummyflag = [False]

    timekeeprec.iteration = addtime(timekeeprec.iteration, 1, carryflag)
    if carryflag[0]:
        timekeeprec.block += 1

    timekeeprec.reportflag = (timekeeprec.nextreport == timekeeprec.iteration)
    if timekeeprec.reportflag:
        timekeeprec.nextreport = addtime(timekeeprec.iteration, timekeeprec.reportperiod, dummyflag)

    timekeeprec.consolereportflag = (timekeeprec.nextconsolereport == timekeeprec.iteration)
    if timekeeprec.consolereportflag:
        timekeeprec.nextconsolereport = addtime(timekeeprec.iteration, timekeeprec.consolereportperiod, dummyflag)

    timekeeprec.plotreportflag = (timekeeprec.nextplotreport == timekeeprec.iteration)
    if timekeeprec.plotreportflag:
        timekeeprec.nextplotreport = addtime(timekeeprec.iteration, timekeeprec.plotreportperiod, dummyflag)

    timekeeprec.gaflag = (timekeeprec.nextga == timekeeprec.iteration)
    if timekeeprec.gaflag:
        timekeeprec.nextga = addtime(timekeeprec.iteration, timekeeprec.gaperiod, dummyflag)


# Report time function
def reporttime(rep, timekeeprec):
    """Print out block and iteration number."""
    rep.write(f"[ Block:iteration ] - [ {timekeeprec.block}:{timekeeprec.iteration} ]\n")


"""
# Example usage
if __name__ == "__main__":
    # Open tfile and rep as needed
    with open("tfile.txt", "r") as tfile, open("report.txt", "w") as rep:
        timekeeprec = TRecord()
        inittimekeeper(tfile, timekeeprec)
        initreptimekeeper(rep, timekeeprec)
        
        # Perform timekeeping and report
        for _ in range(100000):
            timekeeper(timekeeprec)
            if timekeeprec.reportflag:
                reporttime(rep, timekeeprec)
"""


# ========================= end of TIMEKEEP.PAS =========================


# advance.py: Advance variables for the next time step

def advance(clearingrec):
    """
    Advance the winner for the next time step.
    
    Args:
        clearingrec (crecord): The clearinghouse record containing winner information.
    """
    clearingrec.oldwinner = clearingrec.winner


# ========================= end of ADVANCE.PAS =========================

"""
Example usage:
# Create a sample clearingrec object
class ClearingRecord:
    def __init__(self, winner, oldwinner):
        self.winner = winner
        self.oldwinner = oldwinner

# Initialize a clearingrec object
sample_clearingrec = ClearingRecord(winner=2, oldwinner=1)

# Display initial values
print("Initial winner:", sample_clearingrec.winner)
print("Initial oldwinner:", sample_clearingrec.oldwinner)

# Advance winner
advance(sample_clearingrec)

# Display updated values
print("Updated winner:", sample_clearingrec.winner)
print("Updated oldwinner:", sample_clearingrec.oldwinner)
"""
# =============== end of advance.scs ===========


# ========================== begin of ga.py =======================


maxmating = 10


class MRecord:
    def __init__(self):
        self.mate1 = 0
        self.mate2 = 0
        self.mort1 = 0
        self.mort2 = 0
        self.sitecross = 0


class MArray:
    def __init__(self):
        self.records = [MRecord() for _ in range(maxmating)]
    def __getitem__(self, item):
        return self.records[item]


class GRecord:
    def __init__(self):
        self.proportionselect = 0.0
        self.pmutation = 0.0
        self.pcrossover = 0.0
        self.ncrossover = 0
        self.nmutation = 0
        self.crowdingfactor = 0
        self.crowdingsubpop = 0
        self.nselect = 0
        self.mating = MArray()  # mating records for ga report


garec = GRecord()


def initga(gfile, garec, population):
    garec.proportionselect = float(gfile.readline()[:15])
    garec.pmutation = float(gfile.readline()[:15])
    garec.pcrossover = float(gfile.readline()[:15])
    garec.crowdingfactor = int(gfile.readline()[:15])
    garec.crowdingsubpop = int(gfile.readline()[:15])
    garec.nselect = round(garec.proportionselect * population.nclassifier * 0.5)
    garec.nmutation = 0
    garec.ncrossover = 0


def initrepga(rep, garec):
    rep.write("\n")
    rep.write("Genetic Algorithm Parameters\n")
    rep.write("----------------------------\n")
    rep.write(f"Proportion to select/gen  = {garec.proportionselect:8.4f}\n")
    rep.write(f"Number to select          = {garec.nselect:8}\n")
    rep.write(f"Mutation probability      = {garec.pmutation:8.4f}\n")
    rep.write(f"Crossover probability     = {garec.pcrossover:8.4f}\n")
    rep.write(f"Crowding crowdingfactor   = {garec.crowdingfactor:8}\n")
    rep.write(f"Crowding subpopulation    = {garec.crowdingsubpop:8}\n")


def select(population):
    partsum = 0.0
    j = 0
    rand = random() * population.sumstrength
    while partsum < rand and j < population.nclassifier:
        j += 1
        partsum += population.classifier.classes[j].strength
    return j - 1  # -1


def mutation(positionvalue, garec):
    # mutate a single position with specified probability 
    tempmutation = int()
    if flip(garec.pmutation):
        tempmutation = (positionvalue + rnd(1, 2) + 1) % 3 - 1
        garec.nmutation += 1
    else:
        tempmutation = positionvalue
    return tempmutation


def bmutation(positionvalue, pmutation, nmutation):
    # mutate a single bit with specified probability  
    tempmutation = int(0)
    if flip(garec.pmutation):
        tempmutation = (positionvalue + 1) % 2
        garec.nmutation += 1
    else:
        tempmutation = positionvalue
    return tempmutation


#def crossover(parent1, parent2, child1, child2,
#              pcrossover, pmutation, nposition, ncrossover, nmutation):
def crossover(parent1, parent2, child1, child2,
              nposition, garec):

    inheritance = float(0)
    j = int(0)
    sitecross = int(0) # locak 

    if flip(garec.pcrossover):
        sitecross = rnd(0, nposition-1)
        # sitecross = rnd(1, nposition)
        garec.ncrossover += 1
    else:
        #sitecross = nposition + 1  # transfer, but no cross
        sitecross = nposition   # transfer, but no cross
        

    # transfer action part regardless of sitecross
    child1.a = bmutation(parent1.a, garec.pmutation, garec.nmutation)
    child2.a = bmutation(parent2.a, garec.pmutation, garec.nmutation)

    # transfer and cross above cross site
    j = sitecross
    #while j <= nposition:
    while j <= nposition-1:
        child2.c[j] = mutation(parent1.c[j], garec)
        child1.c[j] = mutation(parent2.c[j], garec)
        j += 1

    #j = 1
    j = 0

    # transfer only below cross site
    while j < sitecross:
        child1.c[j] = mutation(parent1.c[j], garec)
        child2.c[j] = mutation(parent2.c[j], garec)
        j += 1

    inheritance = avg(parent1.strength, parent2.strength)
    child1.strength = inheritance
    child1.matchflag = False
    child1.ebid = 0
    child1.bid = 0
    child1.specificity = countspecificity(child1.c, nposition)

    child2.strength = inheritance
    child2.matchflag = False
    child2.ebid = 0
    child2.bid = 0
    child2.specificity = countspecificity(child2.c, nposition)

    return sitecross


def worstofn(population, n):
    # select worst individual from random subpopulation of size n
    j = 0
    worst = 0
    candidate = 0
    worststrength = 0.0

    # Initialize with random selection
    worst = rnd(1, population.nclassifier)-1  # !!!
    worststrength = population.classifier.classes[worst].strength
    # select and compare from remaining subpopulation
    if n > 1:
        for j in range(2, n + 1):
            candidate = rnd(1, population.nclassifier)-1
            if worststrength > population.classifier.classes[candidate].strength:
                worst = candidate
                worststrength = population.classifier.classes[worst].strength
    return worst


def matchcount(classifier1, classifier2, nposition):
    # count number of positions of similarity
    tempcount = 1 if classifier1.a == classifier2.a else 0
    for j in range(0, nposition):
        if classifier1.c[j] == classifier2.c[j]:
            tempcount += 1
    return tempcount


def crowding(child, population, crowdingfactor, crowdingsubpop):
    # replacement using modified De Jong crowding
    matchmax = -1
    mostsimilar = 0
    if crowdingfactor < 1:
        crowdingfactor = 1
    for _ in range(0, crowdingfactor):
        popmember = worstofn(population, crowdingsubpop)
        # pick worst of n
        match = matchcount(child, population.classifier.classes[popmember], population.nposition)
        if match > matchmax:
            matchmax = match
            mostsimilar = popmember
    return mostsimilar


def statistics(population):
    j = 0
    population.maxstrength = population.classifier.classes[0].strength
    population.minstrength = population.classifier.classes[0].strength
    population.sumstrength = population.classifier.classes[0].strength
    for j in range(1, population.nclassifier):  # 0 indisi yukarida ilk deger olaran alindigi icin 1
        population.maxstrength = max(population.maxstrength, population.classifier.classes[j].strength)
        population.minstrength = min(population.minstrength, population.classifier.classes[j].strength)
        population.sumstrength += population.classifier.classes[j].strength
    population.avgstrength = population.sumstrength / population.nclassifier


def ga(garec, population):
    # coordinate selection, mating, crossover, mutation, & replacement
    j = 0
    child1 = ClassType()
    child2 = ClassType()
    statistics(population)  # get average, max, min, sumstrength
    for j in range(0, garec.nselect):
        garec.mating[j].mate1 = select(population)
        garec.mating[j].mate2 = select(population)
        garec.mating[j].sitecross = crossover(
            population.classifier.classes[garec.mating[j].mate1],
            population.classifier.classes[garec.mating[j].mate2],
            child1,
            child2,
            population.nposition,
            garec
        )  # cross 4 mutate
        garec.mating[j].mort1 = crowding(child1, population, garec.crowdingfactor, garec.crowdingsubpop)
        population.sumstrength = population.sumstrength - population.classifier.classes[garec.mating[j].mort1].strength + child1.strength  # update sumstrength
        population.classifier.classes[garec.mating[j].mort1] = child1  # insert child in morel's place

        garec.mating[j].mort2 = crowding(child2, population, garec.crowdingfactor, garec.crowdingsubpop)
        population.sumstrength = population.sumstrength - population.classifier.classes[garec.mating[j].mort2].strength + child2.strength  # update sumstrength
        population.classifier.classes[garec.mating[j].mort2] = child2  # insert child in morel's place


def reportga(rep, garec, population):
    j = 0
    page(rep)
    rep.write("\n")
    rep.write("Genetic Algorithm Report\n")
    rep.write("------------------------\n")
    rep.write("\n")
    rep.write("Pair (Mate1,Mate2) SiteCross Mort1 Mort2\n")
    rep.write("--------------------------------------\n")
    for j in range(0, garec.nselect):
        rep.write(
            f"{j+1:3})  ({garec.mating[j].mate1+1:3},{garec.mating[j].mate2+1:3}) >{garec.mating[j].sitecross+1:3}< ({garec.mating[j].mort1+1:3},{garec.mating[j].mort2+1:3})\n"
        )
    rep.write("\n")
    rep.write(" Statistics Report       \n")
    rep.write("-------------------------\n")
    rep.write(f" Average    strength   = {population.avgstrength:8.2f}\n")
    rep.write(f" Maximum    strength   = {population.maxstrength:8.2f}\n")
    rep.write(f" Minimum    strength   = {population.minstrength:8.2f}\n")
    rep.write(f" Sum of     strength   = {population.sumstrength:8.2f}\n")
    rep.write(f" Number of crossings   = {garec.ncrossover:8}\n")
    rep.write(f" Number of mutations   = {garec.nmutation:8}\n")


# ========================== End of ga.py =======================

# report.py: Report coordination routines

# Data declarations

# Function to send report header to specified file/dev.
def reportheader(rep):
    """
    Send the report header to the specified file/device.
    
    Args:
        rep (TextIO): The report file object.
    """
    page(rep)
    rep.write("Snapshot Report\n")
    rep.write("---------------\n\n")


# Report coordination routine
def report(rep):
    """
    Generate a comprehensive report containing various information.
    
    Args:
        rep (TextIO): The report file object.
    """
    reportheader(rep)
    reporttime(rep, timekeeprec)
    reportenvironment(rep, environrec)
    reportdetectors(rep, envmessage, population.nposition)
    reportclassifiers(rep, population)
    reportaoc(rep, clearingrec)
    reportreinforcement(rep, reinforcementrec)


# Write console report
def consolereport(reinforcementrec):
    """
    Display a console report containing reinforcement information.
    
    Args:
        reinforcementrec (rrecord): The reinforcement record containing information.
    """
    clear_console()  # Clear the console screen
    print(f"|----------------- = -----------------------------|")
    print(f"         iteration = {reinforcementrec.totalcount: 8}")
    print(f"         P correct = {reinforcementrec.proportionreward:8.6f}")
    print(f"       P50 correct = {reinforcementrec.proportionreward50:8.6f}")
    print(f"|----------------- = -----------------------------|")
    print(f"   pop.minstrength = {population.minstrength:8.6f}")
    print(f"   pop.avgstrength = {population.avgstrength:8.6f}")
    print(f"   pop.maxstrength = {population.maxstrength:8.6f}")
    print(f"   pop.sumstrength = {population.sumstrength:8.6f}")
    print(f"|----------------- = -----------------------------|")
    print(f" matchlist.nactive = {matchlist.nactive}")
    print(f"    matchlist.clist= {[j+1 for j in matchlist.clist[0:10]]}")
    print(f"|----------------- = -----------------------------|")


# Write plot report to pfile
def plotreport(pfile, reinforcementrec):
    """
    Write a plot report containing reinforcement information to a file.
    
    Args:
        pfile (TextIO): The plot file object.
        reinforcementrec (rrecord): The reinforcement record containing information.
    """
    pfile.write(f"{reinforcementrec.totalcount:8} "
                f"{reinforcementrec.proportionreward:8.6f} "
                f"{reinforcementrec.proportionreward50:8.6f}\n")


# =============== End of report.py =========================


# initial.py: Initialization coordination


# Function to write a header to the specified file/device.
def initrepheader(rep):
    """
    Write a header to the specified file/device.
    
    Args:
        rep (TextIO): The file object for writing.
    """

    rep.write("********************************************\n")
    rep.write("   A Simple Classifier System - SCS\n")
    rep.write("   (C) David E. Goldberg, 1987\n")
    rep.write("   All Rights Reserved\n")
    rep.write("         in Python    \n")
    rep.write("********************************************\n")
    rep.write("\n\n")


# Clear screen and print interactive header
def interactiveheader():
    """
    Clear the screen and print an interactive header.
    """
    clear_console()  # Clear the console screen
    initrepheader(sys.stdout)


# Initialization coordination
def initialization():
    global cfile, efile, rfile, tfile, gfile, rep, pfile
    global fncfile, fnefile, fnrfile, fntfile, fngfile, fnrep, fnpfile

    """
    Coordinate input and initialization.
    """
    interactiveheader()
    randomize()  # Random number and normal initialization
    initrandomnormaldeviate()

    qtype = QueryType.BATCH  # Use batch input
    cfile = open_input(qtype, '   classifier   ', fncfile)
    efile = open_input(qtype, '   environment  ', fnefile)
    rfile = open_input(qtype, '   reinforcement', fnrfile)
    tfile = open_input(qtype, '   timekeeper   ', fntfile)
    gfile = open_input(qtype, '   gen algorithm', fngfile)
    rep = open_output(qtype, '   report       ', fnrep)
    pfile = open_output(qtype, '   plot file    ', fnpfile)

    # Segment initialization: classifiers, environment, detectors, aoc, reinforcement, timekeeper, ga
    initrepheader(rep)
    initclassifiers(cfile, population)
    initrepclassifiers(rep, population)
    initenvironment(efile, environrec)
    initrepenvironment(rep, environrec)
    initdetectors(efile, detectrec)
    initrepdetectors(rep, detectrec)
    initaoc(cfile, clearingrec)
    initrepaoc(rep, clearingrec)
    initreinforcement(rfile, reinforcementrec)
    initrepreinforcement(rep, reinforcementrec)
    inittimekeeper(tfile, timekeeprec)
    initreptimekeeper(rep, timekeeprec)
    initga(gfile, garec, population)
    initrepga(rep, garec)


# ======================= End of initial.py ======================

# Main function
if __name__ == "__main__":
    # breakpoint()
    initialization()  # Initialize the system
    detectors(environrec, detectrec, envmessage)
    report(rep)
    try:
        sayac = 0
        while True:
            sayac += 1
            if (sayac == 49):
                pass
                # print("Başlıyoruz...")

            timekeeper(timekeeprec)
            environment(environrec)

            # detectors(environrec, detectrec, envmessage)
            envmessage = copy.deepcopy(environrec.signal)

            matchclassifiers(population, envmessage, matchlist)
            aoc(population, matchlist, clearingrec)
            effector(population, clearingrec, environrec)
            reinforcement(reinforcementrec, population, clearingrec, environrec)

            if timekeeprec.reportflag:
                report(rep)

            if timekeeprec.consolereportflag:
                consolereport(reinforcementrec)

            if timekeeprec.plotreportflag:
                plotreport(pfile, reinforcementrec)

            advance(clearingrec)

            if timekeeprec.gaflag:
                ga(garec, population)
                if timekeeprec.reportflag:
                    reportga(rep, garec, population)

            # Halt condition
            if halt() or sayac > MAX_ITER:
                break

    except KeyboardInterrupt:
        print("Kod Ctrl+C ile kesildi.")
    finally:
        # Final report
        print("rapor yapılacak.")
        report(rep)
        print("rapor yapıldı.")
        print("dosyalar kapatılacak.")
        if rndfileopened:
            rndfile.close()

        rep.close()
        pfile.close()  # Close plot file
        cfile.close()
        efile.close()
        rfile.close()
        tfile.close()
        gfile.close()
        print("dosyalar kapatıldı.")


# kenarda dursun
# from os.path import dirname, join
"""
from random import random
from classes import *
from initial import initialization
from environment import environment
from detectors import detectors
from report import report, consolereport, plotreport, reportga
from timekeep import timekeeper, reporttime, advance
"""






