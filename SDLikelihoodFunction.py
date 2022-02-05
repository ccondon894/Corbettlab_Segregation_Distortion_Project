import math
import sys


class SDData:

    def __init__(self, inFile):
        
        self.inFile = inFile
        with open(inFile, 'r') as f:

            self.SD_data = f.readlines()

        

class SDAnalyzer:

    def __init__(self, k, sdPos, data):

        self.data = data
        self.k = k
        self.sdPos = sdPos
        self.data = data
        self.numReads = len(data)

    def sumLikelihood(self, iterations):

        likelihood = 0
        for i in range(iterations):
            likelihood += self.calcLikelihood(self.data[i])

        return likelihood

    def calcLikelihood(self, line):

        # get data from line
        tabSplit = line.split("\t")
        rPos = float(tabSplit[2])
        a_count = int(tabSplit[3])
        A_count = int(tabSplit[4])


        # recombination probability from Haldane's mapping distance
        p_r = 0.5*(1-math.e**(-2 * abs(rPos - self.sdPos)))
        # probability of sequencing error (just an arbitrary number rn)
        p_e = 0.0001

        # likelihood equations for allele A
        e1 = (1 - p_e) * (1 - p_r) * self.k
        e2 = p_e * p_r * self.k
        e3 = (1 - p_e) * p_r * (1 - self.k)
        e4 = p_e * (1 - p_r) * (1 - self.k)
        # likelihood equations for allele a
        e5 = (1 - p_e) * (1 - p_r) * (1 - self.k)
        e6 = (1 - p_e) * p_r * self.k
        e7 = p_e * (1 - p_r) * self.k
        e8 = p_e * p_r * (1 - self.k)

        # calc likelihoods
        A_likelihood = math.log(e1 + e2 + e3 + e4) * A_count
        a_likelihood = math.log(e5 + e6 + e7 + e8) * a_count

        return A_likelihood + a_likelihood

def main():

    dataset_1 = SDData("k0.5.txt").SD_data

    ds1_analysis1 = SDAnalyzer(0.5, 1, dataset_1)
    print(ds1_analysis1.sumLikelihood(ds1_analysis1.numReads))

    ds1_analysis2 = SDAnalyzer(0.25, 1, dataset_1)
    print(ds1_analysis2.sumLikelihood(ds1_analysis2.numReads))


    


if __name__ == "__main__":
    main()




        
