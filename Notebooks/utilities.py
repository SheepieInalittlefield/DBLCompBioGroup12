import random
import math
import numpy as np
import functools
import reedsolo

def OpenFile(file_path): 
    """ PROBABLY UNNECESSARY
    INPUT:
        filepath: path to the file
    OUTPUT:
        message: string containing the data of the file
    
    """
    with open(file_path, 'rb') as file:
        message = file.read()
    return message

def RandomizeMessage(segment, keygenseed = 7) -> bytearray:
    """Randomize all bytes to reduce homopolymers in the final DNA segments
    INPUT: 
        segment: original file, in bytearray format
        keygenseed: integer value of the seed to be used, default 7 
    OUTPUT:
        Randomize the message one byte at a time, using a known keygenseed.
    """
    randomized_segment = bytearray(len(segment))
    for i in range(len(segment)):
        random.seed(keygenseed)
        key = random.randint(0,255)
        keygenseed += 1 
        randomized_segment[i] = key ^ segment[i]
    return randomized_segment

def Segment(message) -> list:
    """ Segment the data into segments of 32 bytes, using padding to make sure the last segment has the correct lenght
    INPUT:
        message: string containing the data of the file
    OUTPUT:
        segment_list: list of binary strings
    """
    n=256
    ## add zero padding and amount of zero padding
    bits_to_add = n - (len(message) % n)
    message += '0'*(bits_to_add-8)
    padding_length = format(bits_to_add-1,'08b')
    message += padding_length
    
    ## segment data
    segments = [message[i:i+n] for i in range(0, len(message), n)]
    return segments

def IdealSoliton(K) -> list: 
    """ Generate a list of probalities of length K, following ideal soliton distribution
    INPUT:
        K: length of list
    OUTPUT:
        probabilities: list of probabilities following ideal soliton distribution
    """
    # initialize with the first two values, p(0) = 0 and p(1) = 1/K
    probabilities = [0, 1/K]
    # calculate the rest of the values using p(i) = 1/(i*(i-1))
    probabilities += [1/(i*(i-1)) for i in range(2, K+1)]
    return probabilities 

def RobustSoliton(K,c,delta) -> list:
    """ Generates a list of probalities of length K, following an robust soliton distribution with variables c and delta.
    INPUT:
        K: length of list
        c: value of c variable in distribution
        delta: value of delta variable in distribution
    OUTPUT:
        probabilities: list of probabilities following robust soliton distribution
    """
    #initialize with the ideal distribution
    probabilities = IdealSoliton(K)
    # Define R
    R = c*(math.log(K/delta)**2)*math.sqrt(K)
    # calculate the additional probabilities
    pivot = int(math.floor(K/R))
    robust_probabilities = [0] + [R/(i*K) for i in range(1, pivot)]
    robust_probabilities += [(R*math.log(R/delta))/K]
    robust_probabilities += [0 for i in range(pivot,K)]
    # add together
    probabilities = np.add(robust_probabilities, probabilities)
    #normalize 
    probabilities /= sum(probabilities)
    return probabilities

def Binary2Bytearray(binary_string): 
    """ Convert a string of binary to bytearray.
    INPUT:
        binary_string: string of binary, so length multiple of 8
    OUTPUT 
        bytearray: bytearray conversion of the string
    """
    return bytearray(int(binary_string,2).to_bytes((len(binary_string)+7) // 8, byteorder='big'))

def PrepareSeed(seed):
    """ Convert seed into a 4 byte long bytearray for packing into droplet
    INPUT:
        seed: integer seed
    OUTPUT 
        seed_index: seed converted to bytearray 
    """
    seed_array = bytearray(Binary2Bytearray(bin(seed)))
    seed_index = bytearray()
    if len(seed_array) < 4:  
        seed_index = bytearray(4-len(seed_array)) + seed_array
    elif len(seed_array) > 4:
        if seed_array[0] == 0:
            del seed_array[0]
            seed_index = seed_array
        else:
            raise Exception("seed too big")
    else: 
        seed_index = seed_array
    return seed_index

def MakeDroplet(randomized_segments, segment_seed, prng, nr_droplets_probabilities) -> bytearray:
    """ Create droplet from seed and segments
    INPUT:
        randomized_segments: randomized and segmented data.
        segment_seed: integer seed for the creation of the droplet.
        prng: random.Random() random number generator.
    OUTPUT:
        droplet_rs: droplet encoded using Reed Solomon
    """
    prng.seed(segment_seed)
    seed_index = PrepareSeed(segment_seed)
    amount = prng.choices(range(0,len(nr_droplets_probabilities)), nr_droplets_probabilities, k = 1)[0]
    segment_indices = prng.sample(range(len(randomized_segments)),k = amount)
    segments = [randomized_segments[i] for i in segment_indices]
    droplet = seed_index + bytearray(functools.reduce(lambda i, j: bytes(a^b for (a, b) in zip(i,j)), segments))
    # prepare reedsolomon
    rsc = reedsolo.RSCodec(2)
    # create the encoded droplet (what will eventually be stored in DNA)
    droplet_rs = rsc.encode(droplet)
    return droplet_rs

def LFSR():
    """Linear forward shifting register based on Galois fields
    The result is an integer between 0 and 32^2-1 which will not repeat until all possible values have been passed
    """
    mask=0b100000000000000000000000011000101 # Do not change this one
    result = 0b101011100 # Set the first state of the register
    nbits = mask.bit_length()-1
    while True:
        result = (result << 1) #Shift the register left once
        xor = result >> nbits #Shift the register right by the amount of bits in the mask -1
        if xor != 0: #XOR is useless if it is 0
            result ^= mask #XOR the state of the register with the mask

        yield result

def Bytearray2Binary(Bytearray) -> str:
    """ Convert bytearray to binary (string)
    INPUT:
        Bytearray: bytearray to be converted
    OUTPUT:
        binary_string: string of converted bytearray
    """
    binary_string = ""
    for i in Bytearray:
        binary_string += format(i, '08b')
    return binary_string

def Binary2DNA(binary_string):
    """ Convert binary (string) to DNA strand (string)
    Input:
        binary_string: binary strand
    Output:
        dna_string: DNA counterpart
    """
    dna_string = ""
    for x in [binary_string[i:i+2] for i in range(0,len(binary_string),2)]:
        if x == "00":
            dna_string += "A"
        if x == "01":
            dna_string += "C"
        if x == "10":
            dna_string += "G"
        if x == "11":
            dna_string += "T"
    return dna_string

def DNA2Binary(dna_string):
    """ Convert DNA strand (string) to binary (string)
    Input:
        dna_string: DNA strand
    Output:
        binary_string: binary counterpart of s
    """
    binary_string = ""
    for x in dna_string:
        if x == "A":
            binary_string += "00"
        if x == "C":
            binary_string += "01"
        if x == "G":
            binary_string += "10"
        if x == "T":
            binary_string += "11"
    return binary_string

def CheckBiochemicalRequirements(s, max_length=3, check=False):
    """ Check if DNA strand has expected structure
    Input:
        s: DNA strand (string) of A,G,C,T
        length: expected length, normally set at 152
    Output:
        check: Boolean, True or False
    """
    if s.find((max_length+1)*'C') ==-1 & s.find((max_length+1)*'G') ==-1 & s.find((max_length+1)*'T')==-1 & s.find((max_length+1)*'A')==-1:
        nr_CG = s.count('C')+s.count('G')
        per_CG=nr_CG/len(s)
        if (per_CG < 0.55) & (per_CG > 0.45):
            check = True
    return check

def CheckOligoLength(s, length=152, check=False):
    """ Check if DNA strand has expected length
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


def RecoverSeed(droplet_seed, total_segments, distribution_size=1000):
    """ Determine segments XORd into this droplet
    INPUT: 
        droplet_seed: seed (integer) of this droplet
        total_segments: total number of segments in input file.
    OUTPUT: 
        amount_recovery: amount of segments in the droplet
        segment_indices: indices of the segments in the droplet
    """
    prng = random.Random()
    prng.seed(droplet_seed)
    amount_recovery = prng.choices(range(0,distribution_size+1), RobustSoliton(distribution_size,0.001,0.025), k = 1)[0]
    segment_indices = prng.sample(range(total_segments), k = amount_recovery)
    return (amount_recovery, segment_indices)

def Decode(input_data, output_data):
    """ Decode the payloads which have one segment not in output_data
    Input:
        inputdata: a list of tuples ([list of segment_indices], payload)
    Output:
        output_data: a dictionary with {nr of segment : value of segment}
    """
    remaining_droplets = []
    newsolves = 0
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
        else:
            remaining_droplets.append((remaining_segments, XOR))
    return remaining_droplets, output_data, newsolves

def RemovePadding(solution):
    """ Remove added padding
    Input:
        solution: Binary string of all solved droplets combined
    Output:
        solution: Binary string of solution without paddding
    """
    bits_padded = int(solution[-8:],2)+1
    solution = solution[:-bits_padded]
    return solution