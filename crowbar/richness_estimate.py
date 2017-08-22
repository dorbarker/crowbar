import random
from collections import Counter
import attr
import pandas as pd

@attr.s
class Population(object):

    alleles = attr.ib(validator=attr.validators.instance_of(pd.Series))

    def monte_carlo(self, replicates: int, seed: int = 1) -> pd.Series:
        '''Monte Carlo estimation of the probability of finding a new class
        on the next sampling.

        Returns a pandas Series of the proportion replicates for which each
        observation yielded a new class
        '''

        random.seed(seed)  # for reproducibility

        new_allele = pd.Series([0 for _ in self.alleles])

        individuals = list(self.alleles)
        number_obs = len(individuals)

        for replicate in range(replicates):

            already_observed = set()

            individuals = random.sample(individuals, number_obs)

            for idx, observation in enumerate(individuals):

                # 1 if True, 0 if False
                new_allele[idx] += observation not in already_observed

                already_observed.add(observation)

        percent_new_allele = new_allele / replicates

        return percent_new_allele

    def proportion_successes(self) -> float:
        '''Returns the probability of an individual representing a new class
        given previous discovery rate of new classes
        '''
        return len(set(self.alleles)) / len(self.alleles)

    def chao1(self) -> float:
        '''Implementation of the Chao 1 estimator

        Nonparametric Estimation of the Number of Classes in a Population
        Anne Chao, 1984

        Built with the help of:
        https://www.uvm.edu/~ngotelli/manuscriptpdfs/Chapter%204.pdf

        Returns the estimated number of classes in a population
        '''

        n_observed = len(set(self.alleles))

        abundances = Counter(self.alleles)
        abundance_counts = list(abundances.values())

        f_1 = abundance_counts.count(1)

        f_2 = abundance_counts.count(2)

        return n_observed + ((f_1 * (f_1 - 1)) / (2 * (f_2 + 1)))

    def ace(self) -> float:
        '''Implementation of the abundace-based coverage estimator (ACE)
        of species richness

        Nonparametric Prediction in Species Sampling
        Anne Chao and Tsung-Jen Shen, 2004

        Built with the help of:
        https://www.uvm.edu/~ngotelli/manuscriptpdfs/Chapter%204.pdf

        Returns the estimated number of classes in a population
        '''

        n_observed = len(set(self.alleles))

        abundances = Counter(self.alleles)
        abundance_counts = pd.Series(list(abundances.values()))

        rare_range = list(range(1, 11))  # 1-10
        abund_range = list(range(11, max(abundance_counts) + 1))

        s_rare = sum(abundance_counts <= 10)

        s_abund = sum(abundance_counts > 10)

        n_rare = sum(abundance_counts[abundance_counts <= 10])

        if n_rare == sum(abundance_counts[abundance_counts == 1]):
            msg = 'ACE is undefined when all rare species are singletons'
            raise ZeroDivisionError(msg)

        singletons = sum(abundance_counts == 1)
        c_ace = 1 - (singletons / n_rare)

        k2 = ((k * (k - 1)) * sum(abundance_counts == k)
              for k in rare_range)


        y_ace = ((s_rare / c_ace) * (sum(k2) / (n_rare * (n_rare - 1)))) - 1

        y_ace = max(y_ace, 0)

        s_ace = s_abund + (s_rare / c_ace) + \
                ((singletons / c_ace) * y_ace)

        return s_ace

