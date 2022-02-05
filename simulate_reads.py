import numpy as np
import argparse

### get args
parser = argparse.ArgumentParser(description='Create SD data.')
optional = parser.add_argument_group('optional arguments')
optional.add_argument("--depth", type=int, default = 100, help="mean sequencing depth per site")
optional.add_argument("--position", type=int, default = 5000, help="position of distorting gene")
optional.add_argument("--morgans", type=float, default = 2, help="total recombination disotance along scaffold")
optional.add_argument("-k", type=float, default = 0.5, help="distortion effect")
optional.add_argument("--error", type=float, default = 0.01, help="sequencing error per site per read")
optional.add_argument("--nsites", type = int, default = 10000, help="total number of sites")
args = parser.parse_args()

### iterate across all sites create read 
for i in range(0, args.nsites) :
	
    ### sequencing depth at the site
    site_depth = np.random.poisson(args.depth)

    ### recombination rate to the distorting position 
    recombination = ( 1 - np.exp( -2 * abs( ( ( args.position / args.nsites ) - ( i / args.nsites ) ) ) * args.morgans ) ) / 2 

    ### sequencing depth of A prior to recombination etc. 
    drawA = np.random.binomial(site_depth, args.k)

    ### depth swaps due to recombination 
    swapA = np.random.binomial( drawA, recombination ) 
    swapa = np.random.binomial( site_depth - drawA, recombination ) 
    drawA = drawA - swapA + swapa

    ### swapping for sequencing error 
    swapA = np.random.binomial( drawA, args.error )
    swapa = np.random.binomial( site_depth - drawA, args.error )
    drawA = drawA - swapA + swapa 

    ### print the counts 
    position_in_morgans = i/args.nsites * args.morgans
    drawa = site_depth - drawA
    print ("scaffold_X", i, position_in_morgans, drawA, drawa, sep="\t")


	

