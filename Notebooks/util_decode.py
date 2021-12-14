import random
import math
import numpy as np
import functools
import reedsolo

def CheckBiochemicalRequirements(s, max_length=3, check=False):
    """ Check if a DNA-strand meets the biochemical requirements
    Input: 
        s: a DNA strand (string) of A,G,C,T 
        max_length: a max_length (int) which is the allowed length for a homopolymer
    Output:
        check: Boolean, True or False
    """
    if s.find((max_length+1)*'C') ==-1 & s.find((max_length+1)*'G') ==-1 & s.find((max_length+1)*'T')==-1 & s.find((max_length+1)*'A')==-1:
#         print('yes, no homopolymers')
        nr_CG = s.count('C')+s.count('G')
        per_CG=nr_CG/len(s)
        if (per_CG < 0.55) & (per_CG > 0.45):
#             print('Suitable CG content')
            check = True
#         else:
#             print('CG content not in accepted ranges')
#     else:
#         print('not a valid sequence: it contains homopolymers longer than three bases')
    return check


def check_oligo_length(s, length=152, check=False):
    """ Check if Dna strand has expected length
    Input:
        s: a DNA strand (string) of A,G,C,T
        length: expected length, normally set at 152
    Output:
        check: Boolean, True or False
    """
    if len(s) == length:
        check = True
    else:
        if len(s) < length:
            print('oligo length is too short')
        else:
            print('oligo length is too long')
    return check


def DNA2Binary(s):
    """ convert DNA strand (string) to binary (string)
    Input:
        s: DNA strand
    Output:
        bs: binary counterpart of s
    """
    bs = ""
    for x in s:
        if x == "A":
            bs = bs + "00"
        if x == "C":
            bs = bs + "01"
        if x == "G":
            bs = bs + "10"
        if x == "T":
            bs = bs + "11"
    return bs

##weghalen
def string_to_bytearray(bs):
    """
    INPUTS: 
    bs: string in utf-8 encoding)
    OUTPUT:
    ba: bytearray 
    Converts string to a bytearray
    """
    ba = bytarray(bs,'utf-8')
    return ba


def Reed_solo(ba):
    """"
    INPUT:
    ba: droplet bytearray with rs encoding.
    OUTPUT:
    decoded_rs: decoded (so errors corrected) droplet 
    Uses RS bits to solve errors in a droplet.
    """
    rsc = reedsolo.RSCodec(2)
    if rsc.check(ba):
        decoded_rs = rsc.decode(ba)[0]
    return decoded_rs

def ideal_soliton(K) -> list: 
    """
    INPUT:
    K: length of list
    OUTPUT:
    probabilities: list of probabilities following ideal soliton distribution
    Generates a list of probalities of length K, following an ideal soliton distribution
    """
    # initialize with the first two values, p(0) = 0 and p(1) = 1/K
    probabilities = [0, 1/K]
    # calculate the rest of the values using p(i) = 1/(i*(i-1))
    probabilities += [1/(i*(i-1)) for i in range(2, K+1)]
    return probabilities

def robust_soliton(K,c,delta) -> list:
    """
    INPUT:
    K: length of list
    c: value of c variable in distribution
    delta: value of delta variable in distribution
    OUTPUT:
    probabilities: list of probabilities following ideal soliton distribution
    Generates a list of probalities of length K, following an robust soliton distribution with variables c and delta.
    """
    #initialize with the ideal distribution
    probabilities = ideal_soliton(K)
    # Define R
    R = c*(math.log(K/delta))*math.sqrt(K)
    # calculate the additional probabilities
    robust_probabilities = [0] + [R/(i*K) for i in range(1, int(K/R)-1)]
    robust_probabilities += [(R*math.log(R/delta))/K]
    robust_probabilities += [0 for i in range(int(K/R),K+1)]
    # add together
    probabilities = np.add(robust_probabilities, probabilities)
    #normalize 
    probabilities /= sum(probabilities)
    return probabilities


def recover_seed(decoded_rs, total_segments):
    """
    INPUT: 
    decoded_rs: decoded droplet, in bytearray format
    total_segments: total number of segments in input file.
    OUTPUT: 
    amount_recovery: amount of segments in the droplet
    segment_indices: indices of the segments in the droplet
    takes a droplet and uses it's seed to determine which segments were XORd into it.
    """
    prng = random.Random()
    droplet_seed = decoded_rs[0:32]
    prng.seed(droplet_seed)
    amount_recovery = prng.choices(range(0,101), robust_soliton(100,0.1,0.05), k = 1)[0]
    segment_indices = prng.sample(range(total_segments), k = amount_recovery)
    return (amount_recovery, segment_indices)

# def input_data(segment_indices, decoded_rs):
#     segment = decoded_rs[32:len(decoded_rs)]
#     segmentlist += [segment_indices, segment]
#     return input_data

def Decoding_step1and2(input_data, output_data):
    """ Decode the segments from the segment_indices and payload
    input:
        inputdata: a list of tuples ([list of segment_indices], payload)
    output:
        output_data: a dictionary with {nr of segment : value of segment}
    """
    remaining_droplets = []
    for droplet in input_data:
        segment_indices = droplet[0]
        XOR = droplet[1]
        remaining_segments = []
        for i in range(len(segment_indices)):
            if segment_indices[i] in output_data:
                dif = int(XOR,2) ^ int(output_data[segment_indices[i]],2)
                XOR = '{0:0{1}b}'.format(dif,len(droplet[1]))
            else:
                remaining_segments.append(segment_indices[i])    

        if len(remaining_segments)==1:
            output_data[remaining_segments[0]] = XOR
            newsolves=1
            while newsolves>0: 
                startsolves = len(output_data)
                output_data = Decoding_step1and2(remaining_droplets, output_data)
                endsolves = len(output_data)
                newsolves = endsolves-startsolves
        else:
            remaining_droplets.append((remaining_segments, XOR))
    return output_data


def unrandomize(solution, keygenseed = 7) -> bytearray:
    """
    INPUT: 
    solution: randomized file, in bytearray format
    keygenseed: integer value of the seed to be used, default 7 
    OUTPUT:
    Unrandomize the message one byte at a time, using a known keygenseed.
    """
    unrandomized = bytearray(len(solution))
    for i in range(len(solution)):
        key = generate_key(keygenseed)
        keygenseed += 1
        unrandomized[i] = key ^ solution[i]
    return unrandomized

def original_data(unrandomized):
    original_data = unrandomized.decode()
    return original_data