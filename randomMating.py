import copy
import numpy as np
import random
import sys
import unittest
from collections import defaultdict
from bisect import bisect

def weighted_choice(choices):
    ''' Returns a randomly chosen value from the weighted distribution, choices.
    Choices is an iterable containing tuples of values and integer weights,
    ( (v1,w1), (v2,w2), ...) '''
    values, weights = zip(*choices)
    total = 0
    cum_weights = []
    for w in weights:
        total += w
        cum_weights.append(total)
    x = random.random() * total
    i = bisect(cum_weights, x)
    return values[i]

class Gene(tuple):
    ''' diploid gene (genotype) implemented as a 2-tuple
        g = Gene(5,4) --> (4,5)
    '''
    def __new__(cls, *alleles):
        ''' In order to implement Gene(5,4) constructor form with tuple,
        need to use __new__. If willing to let constructor be of form Gene((5,4)),
        can implement more typically using __init__
        '''
        assert len(alleles), 2
        assert all(isinstance(al, int) for al in alleles) # alleles are integers
        # sorted returns a list; force it to a tuple.  Rest is just magic
        # cls plays similar role to self here
        return super(Gene, cls).__new__(cls, tuple(sorted(alleles)))
    def allele_count(self):
        return { k:self.count(k) for k in self }
    def mate(self, other):
        '''
        Given genes (a1, a2) and (a3, a4) randomly return one of the four alternatives:
            (a1, a3), (a1, a4), (a2, a3), (a2, a4)
        '''
        # sorted returns a list -- need to force a tuple
        return Gene(random.choice(self), random.choice(other))

class Genome(tuple):
    '''
    A Genome object is a tuple of gene objects (genotypes).  The index of each gene
    in the Genome tuple is its locus
    '''
    def __new__(cls, *genotypes):
        '''
        Genome((1,1),(2,2),(3,3))
        '''
        genes = [Gene(g[0],g[1]) for g in genotypes]
        return super(Genome, cls).__new__(cls, tuple(genes))
    def allele_count(self):
        return [ g.allele_count() for g in self ]
    def mate(self, other):
        ''' Returns a diploid Genome from the random mating of self with other.
        '''
        assert isinstance(other, Genome)
        assert len(self) == len(other)
        # Note: the * is required to turn the list into a list of args
        return Genome(*[s.mate(o) for (s,o) in zip(self,other)])

class Species(defaultdict):
    '''
    Represents an hermaphroditic species as a dict whose keys are Genome objects and
    whose values are integers representing the population of that Genome in the species
    '''
    def __init__(self, genome_dict):
        ''' 
        genome_dict keys are genome objects; values are integers, representing the number
        of individuals in the species with that genome.  Internally, store as defaultdict
        in order to allow easy incrementing as new genomes are created through mating
        '''
        assert all(isinstance(g, Genome) for g in genome_dict)
        # Pass genome_dict to the defaultdict __init__ method specifying int factory
        super(Species, self).__init__(int, genome_dict)
    def population(self):
        return sum(self.values())
    def genome_frequencies(self):
        total = float(sum(self.values()))
        return { g:cnt/total for g,cnt in self.items() }
    def allele_counts(self, locus=0):
        ''' Return the distribution of alleles at *locus* across all genomes in the species '''
        # first collect the allele dists at locus for each genome
        alleles_by_genome = {g:g.allele_count()[locus] for g in self}
        # next sum up number of each allele in each genome weighted
        # by pop of each g in species 
        allele_sum = defaultdict(int)
        for g in alleles_by_genome:
            for allele in alleles_by_genome[g]:
                allele_sum[allele] += alleles_by_genome[g][allele] * self[g] 
        return allele_sum
    def allele_frequencies(self, locus=0):
        total = float(sum(self.allele_counts(locus).values()))
        return { al:cnt/total for al,cnt in self.allele_counts(locus).items() }
    def mate(self):
        # make a working copy of the species object
        s = Species(self) 
        while s.population() >= 2:
            # randomly choose one genome from copy and decrement population in copy
            g1 = weighted_choice(s.items())
            s[g1] -= 1
            # randomly choose another genome from copy and decrement copy
            g2 = weighted_choice(s.items())
            s[g2] -= 1
            # mate g1 and g2 and add offspring to species population
            offspring = g1.mate(g2)
            self[offspring] += 1  # defaultdict used here

class TestRM(unittest.TestCase):
    def setUp(self):
        self.g1 = Gene(1,1)
        self.g2 = Gene(2,2)
        self.genome1 = Genome(self.g1)
        self.genome2 = Genome(self.g2)
        self.sp0 = Species({self.genome1:1, self.genome2:1})
    def test_2genotype_2allele(self):
        sp0 = self.sp0
        genome1 = self.genome1
        genome2 = self.genome2
        self.assertEqual(sp0.population(), 2)
        self.assertEqual(sp0.genome_frequencies(), {genome1: 0.5, genome2: 0.5})
        self.assertEqual(sp0.allele_counts(), defaultdict(int, {1: 2, 2: 2}))
        self.assertEqual(sp0.allele_frequencies(), {1: 0.5, 2: 0.5})
        sp0.mate()
        self.assertEqual(sp0.population(), 3)
        self.assertEqual(sp0.genome_frequencies(), {genome1: 1./3,
                                        Genome((1, 2)): 1./3,
                                        genome2: 1./3}
                         )
        self.assertEqual(sp0.allele_frequencies(), {1: 0.5, 2: 0.5})
        self.assertEqual(sp0.allele_counts(), defaultdict(int, {1: 3, 2: 3}))


# Tests
if __name__ == '__main__':
    unittest.main()
    trials = 10
    generations = 10
    for n in range(trials):
        sp0 = Species({genome1:23, genome2:11})
        for i in range(generations):
            sp0.mate()
        print(sp0.population())
        print(sp0.allele_frequencies())
        print(sp0.genome_frequencies())
