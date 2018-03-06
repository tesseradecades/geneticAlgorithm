import random

class Gene:
    """
    A class that encapsulates a boolean value which can be output in a digital
        representation("0" or "1"), and easily flipped(negated)
    @param bit: str - a string of "0" or "1" used to set the bit
    """
    def __init__(self, bit: str):    
        self.setBit(bit)
    
    """
    getter for the digital representation of the bit "0" or "1"
    """
    def bit(self) -> str:
        return str(self)
    """
    getter for the boolean representation of the bit False or True
    """
    def trueBit(self) -> bool:
        return self._bit
    
    """
    @param bit: str - a string of either "0" or "1" used to set the bit to
        False or True respectively
    """
    def setBit(self, bit: str):
        self._bit = (bit == "1")
    def __str__(self) -> str:
        if(self._bit == True):
            return "1"
        else:
            return "0"

    """
    Negates the gene's underlying boolean value, setting it from True to False,
        or from False to True
    """
    def flip(self):
        self._bit = not self._bit

    """
    Method to create a Gene from another gene. Necessary for creating snapshots
        of Genes, Chromosomes, and Populations
    """
    def initFromGene(gene):
        return Gene(gene.bit())

class Chromosome:
    """
    A Chromosome wraps a list of genes.
    @param genes: list - a collection of genes that this class wraps
    """
    def __init__(self, genes: list):
        self.setGenes(genes)
    """
    Initializes a chromosome from a binary string meant to represent the 
    chromosome's genes.
    @param genes: str - a binary string meant to represent a chromosome's genes.
    @return chromosome: Chromosome
    """
    def initFromString(genes: str):
        ret = []
        for g in genes:
            ret.append(Gene(g))
        return Chromosome(ret)
    """
    Getter for the list of Genes wrapped by this class
    """
    def genes(self) -> list:
        return self._genes
    """
    Setter for the list of genes wrapped by this class
    """
    def setGenes(self, genes: list):
        self._genes = genes
    """
    Flips the bit of a gene at a particular index of the list wrapped by this
    class.
    @param bit: int - the index of the gene that need s its bit flipped
    """
    def flipBit(self, bit: int):
        self._genes[bit].flip()
    def __repr__(self):
        return str(self)
    def __str__(self) -> str:
        ret = ""
        for g in self._genes:
            ret += str(g)
        return ret
    """
    Initializes a Chromosome as a copy of an existing Chromosome
    """
    def initFromChromosome(chromosome):
        ret = []
        for g in chromosome.genes():
            ret.append(Gene.initFromGene(g))
        return Chromosome(ret)

class Population:
    """
    This class wraps a list of chromosomes, meant to be used as the population
    in a genetic algorithm.
    @param chromosomes: list - the list of Chromosomes wrapped by this class.
    @return Population
    """
    def __init__(self, chromosomes: list):
        self.setChromosomes(chromosomes)
    """
    Getter for the list of chromosomes wrapped by this class
    @return chromosomes: list
    """
    def chromosomes(self) -> list:
        return self._chromosomes
    """
    Setter for the list of chromosomes wrapped by this class
    @param chromosomes: list
    """
    def setChromosomes(self, chromosomes: list):
        self._chromosomes = chromosomes
    def __repr__(self):
        return str(self)
    def __str__(self) -> str:
        ret = ""
        for c in self._chromosomes:
            ret += str(c)
        return ret
    """
    Returns a list of the string representations of the chromosomes in the list
    wrapped by this class
    @return ret: list
    """
    def asList(self) -> list:
        ret = []
        for c in self._chromosomes:
            ret.append(str(c))
        return ret

    """
    Creates and returns one Population as an exact copy of another Population
    down to the Gene level.
    @param population: Population - the Population to be copied
    @return Population - an exact copy of the population param
    """
    def initFromPopulation(population):
        ret = []
        for c in population.chromosomes():
            ret.append(Chromosome.initFromChromosome(c))
        return Population(ret)


class Context:
    """
    Contextually, for a genetic algorithm to work, it needs an initial 
    population, and a fitness function.
    @param fitnessFunction - the function that is used to calculate the fitness
        of a particular chromosome
    @param initialPopulation - the Population of chromosomes that the algorithm
        will start from
    @return context        
    """
    def __init__(self, fitnessFunction, initialPopulation: Population):
        self.setFitnessFunction(fitnessFunction)
        self.setInitialPopulation(initialPopulation)
    
    #getter and setter for the fitness function
    def fitnessFunction(self):
        return self._fitnessFunction
    def setFitnessFunction(self, fitnessFunction):
        self._fitnessFunction = fitnessFunction
    
    #getter and setter for the initial population
    def initialPopulation(self) -> Population:
        return self._initialPopulation
    def setInitialPopulation(self,initialPopulation: Population):
        self._initialPopulation = initialPopulation

"""
Method to pick the two most fit chromosomes and two least fit chromosomes out
of a dictionary.
@param fitness: dict - A dictionary with Chromosomes as keys, and floats as values
@return :list - a list containing the two Chromosomes that map to the two 
    highest values in @fitness @return[0,1], as well as the two Chromosomes 
    that map to the two lowest values in @fitness @return[2,3].
"""
def selection(fitnesses: dict) -> list:
    chromosomes = list(fitnesses.keys())

    #get the most and least fit
    fittest = chromosomes[0]
    leastFit = chromosomes[1]
    for chromosome in chromosomes:
        fit = fitnesses[chromosome]
        if(fit > fitnesses[fittest]):
            fittest = chromosome
        elif(fit < fitnesses[leastFit]):
            leastFit = chromosome
    chromosomes.remove(fittest)
    chromosomes.remove(leastFit)

    #get the 2nd most and 2nd least fit
    fittest2 = chromosomes[1]
    leastFit2 = chromosomes[0]
    for chromosome in chromosomes:
        fit = fitnesses[chromosome]
        if(fit > fitnesses[fittest2]):
            fittest2 = chromosome
        elif(fit < fitnesses[leastFit2]):
            leastFit2 = chromosome

    return [fittest,fittest2,leastFit2,leastFit]

"""
Method to simulate the mating of two chromosomes. At a random point in the 
Chromosome pair, each pair is cut in half, and the end halves of the 
Chromosomes are swapped, creating two "child" Chromosomes. Each child 
Chromosome then has a chance to undergo mutation, before being returned by the
method.
@param parent1: Chromosome - A Chromosome representing one of the parents to be
    mated.
@param parent2: Chromosome - A Chromosome representing the other parent to be
    mated.
@param rateOfMutation: float - The probability of each gene mutating in a child
    Chromosome
@return : tuple - the two child Chromosomes
"""
def crossover(parent1: Chromosome, parent2: Chromosome, rateOfMutation: float
    ) -> tuple:
    parent1String = str(parent1)
    parent2String = str(parent2)
    pointOfCrossingOver = random.randint(1,len(parent1String))
    
    child1 = Chromosome.initFromString(parent1String[:pointOfCrossingOver] +
    parent2String[pointOfCrossingOver:])
    
    mutation(child1, rateOfMutation)
    
    child2 = Chromosome.initFromString(parent2String[:pointOfCrossingOver] +
    parent1String[pointOfCrossingOver:])

    mutation(child2,rateOfMutation)
    return (child1,child2)

"""
A method representing the act of mutating genes in a Chromosome.
@param chromosome: Chromosome - the Chromosome to be mutated
@param rateOfMutation: float - The probability of each gene mutating in 
    chromosome
"""
def mutation(chromosome: Chromosome, rateOfMutation: float):
    i = 0
    while i < len(str(chromosome)):
        if(random.randint(1,101) < 100*rateOfMutation):
            chromosome.flipBit(i)
        i+=1

"""
A method used to determine when a Population's genetic diversity has reached an 
acceptable level of convergence.
@param parentGeneration: Population - the generation before the most recent 
    generation.
@param childGeneration: Population - the most recent generation to be compared
@param targetHammingDistance: int - An integer representing the number of 
    genetic differences between the parentGeneration and childGeneration 
    necessary for the populations to qualify as having converged.
@return : bool - A boolean that is True if the Population has converged, and
    False otherwise
"""
def termination(parentGeneration: Population, childGeneration: Population, 
    targetHammingDistance: int) -> bool:
    hammingDistance = 0
    i = 0
    while(i < len(str(parentGeneration))):
        if(str(parentGeneration)[i] != str(childGeneration)[i]):
            hammingDistance += 1
        i+=1
    return hammingDistance <= targetHammingDistance

"""
A method to compute the Hamming Distance necessary for termination. The current
formula returns the average number of mutations expected in a child Population,
relative to its parent Population. This is working under the assumption that 
the only genetic changes between the Populations are due to mutation.
"""
def getTargetHammingDistance(rateOfMutation: float, chromosomeLength: int
    ) -> int:
    return rateOfMutation*chromosomeLength*2

"""
A method to run the Genetic Algorithm after computing a Target Hamming Distance
based on the rateOfMutation and the context
@param rateOfMutation: float - The probability of each gene mutating in 
    chromosome
@param context: Context - The context in which the algorithm is to be run,
    which should contain the starting Population, and the FitnessFunction
    that will be used to judge the fitness of Chromosomes.
@return : Chromosome - a the most fit Chromosome at the time the Population
    reaches satisfactory convergence.
"""
def geneticAlgorithm(rateOfMutation: float, context: Context) -> Chromosome:
    return geneticAlgorithmHelper(
        rateOfMutation, 
        getTargetHammingDistance(
            rateOfMutation, 
            len(str(context.initialPopulation().chromosomes()[0]))), 
        context)

"""
A method to run the Genetic Algorithm based on the rateOfMutation 
targetHammingDistance, and the context. 
1) Make copy of initial Population, known as the current generation.
    1a) and a snapshot of the current generation.
2) Compute fitness for each Chromosome in the Population
3) Select the two most fit Chromosomes, and two least fit Chromosomes
    3a) Remove the two least fit Chromosomes from the Population
    3b) Generate two Child Chromosomes from the two most fit and add them to
        the Population
4) Repeat steps 1a to 3b until the genetic diversity has converged.
5) Return the most fit individual
@param rateOfMutation: float - The probability of each gene mutating in 
    chromosome
@param context: Context - The context in which the algorithm is to be run,
    which should contain the starting Population, and the FitnessFunction
    that will be used to judge the fitness of Chromosomes.
@return : Chromosome - a the most fit Chromosome at the time the Population
    reaches satisfactory convergence.
"""
def geneticAlgorithmHelper(rateOfMutation: float, targetHammingDistance: int, 
    context: Context) -> Chromosome:
    generation = Population.initFromPopulation(context.initialPopulation())
    snapshot = Population.initFromPopulation(generation)
    fitnessDict = {}
    for chromosome in generation.chromosomes():
        fitnessDict[chromosome] = context.fitnessFunction()(chromosome)
    selected = selection(fitnessDict)
    #kill the least fit
    generation.chromosomes().remove(selected[2])
    generation.chromosomes().remove(selected[3])

    #mate the most fit
    children = crossover(selected[0], selected[1], rateOfMutation)
    generation.chromosomes().append(children[0])
    generation.chromosomes().append(children[1])

    while(not termination(snapshot, generation,targetHammingDistance)):
        snapshot = Population.initFromPopulation(generation)
        fitnessDict = {}
        for chromosome in generation.chromosomes():
            fitnessDict[chromosome] = context.fitnessFunction()(chromosome)
        selected = selection(fitnessDict)
        #kill the least fit
        generation.chromosomes().remove(selected[2])
        generation.chromosomes().remove(selected[3])

        #mate the most fit
        children = crossover(selected[0], selected[1], rateOfMutation)
        generation.chromosomes().append(children[0])
        generation.chromosomes().append(children[1])
    return selected[0]

def main():
    g = Gene("1")
    c = Chromosome([g])
    p = Population([c])
    print(Gene.initFromGene(g))
    print(Chromosome.initFromChromosome(c))
    print(Population.initFromPopulation(p))
    return

if __name__ == "__main__":
    main()
