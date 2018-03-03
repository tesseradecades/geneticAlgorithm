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
    def __init__(self, genes: list):
        self.setGenes(genes)
    def initFromString(genes: str):
        ret = []
        for g in genes:
            ret.append(Gene(g))
        return Chromosome(ret)
    
    """def __init__(self, genes: int):
        self.__init__(str(genes))"""

    def genes(self) -> list:
        return self._genes
    def setGenes(self, genes: list):
        self._genes = genes
    def flipBit(self, bit: int):
        self._genes[bit].flip()
    def __repr__(self):
        return str(self)
    def __str__(self) -> str:
        ret = ""
        for g in self._genes:
            ret += str(g)
        return ret

    def initFromChromosome(chromosome):
        ret = []
        for g in chromosome.genes():
            ret.append(Gene.initFromGene(g))
        return Chromosome(ret)

class Population:
    def __init__(self, chromosomes: list):
        self.setChromosomes(chromosomes)
    def chromosomes(self) -> list:
        return self._chromosomes
    def setChromosomes(self, chromosomes: list):
        self._chromosomes = chromosomes
    def __repr__(self):
        return str(self)
    def __str__(self) -> str:
        ret = ""
        for c in self._chromosomes:
            ret += str(c)
        return ret
    
    def asList(self) -> list:
        ret = []
        for c in self._chromosomes:
            ret.append(str(c))
        return ret

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

def selection(fitnesses: dict) -> list:
    chromosomes = list(fitnesses.keys())
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

    fittest2 = chromosomes[1]
    leastFit2 = chromosomes[0]
    for chromosome in chromosomes:
        fit = fitnesses[chromosome]
        if(fit > fitnesses[fittest2]):
            fittest2 = chromosome
        elif(fit < fitnesses[leastFit2]):
            leastFit2 = chromosome

    return [fittest,fittest2,leastFit2,leastFit]

def crossover(parent1: Chromosome, parent2: Chromosome, rateOfMutation: float) -> tuple:
    parent1String = str(parent1)
    parent2String = str(parent2)
    pointOfCrossingOver = random.randint(1,len(parent1String))
    child1 = Chromosome.initFromString(parent1String[:pointOfCrossingOver]+parent2String[pointOfCrossingOver:])
    mutation(child1, rateOfMutation)
    child2 = Chromosome.initFromString(parent2String[:pointOfCrossingOver]+parent1String[pointOfCrossingOver:])
    mutation(child2,rateOfMutation)
    return (child1,child2)

def mutation(chromosome: Chromosome, rateOfMutation: float):
    i = 0
    while i < len(str(chromosome)):
        if(random.randint(1,101) < 100*rateOfMutation):
            chromosome.flipBit(i)
        i+=1

def termination(parentGeneration: Population, childGeneration: Population, targetHammingDistance: int) -> bool:
    hammingDistance = 0
    i = 0
    while(i < len(str(parentGeneration))):
        if(str(parentGeneration)[i] != str(childGeneration)[i]):
            hammingDistance += 1
        i+=1
    return hammingDistance <= targetHammingDistance

def geneticAlgorithm(rateOfMutation: float, targetHammingDistance: int, context: Context) -> Chromosome:
    i = 1
    print("\n\nGeneration #"+str(i))
    generation = Population.initFromPopulation(context.initialPopulation())
    print("Population:\t",generation.asList())
    snapshot = Population.initFromPopulation(generation)
    fitnessDict = {}
    for chromosome in generation.chromosomes():
        fitnessDict[chromosome] = context.fitnessFunction()(chromosome)
    print("Fitnesses:\t",fitnessDict)
    selected = selection(fitnessDict)
    print("Selections:\t",selected)
    #kill the least fit
    generation.chromosomes().remove(selected[2])
    generation.chromosomes().remove(selected[3])

    #mate the most fit
    children = crossover(selected[0], selected[1], rateOfMutation)
    generation.chromosomes().append(children[0])
    generation.chromosomes().append(children[1])

    #print(generation.asList())

    while(not termination(snapshot, generation,targetHammingDistance)):
        i+=1
        print("\n\nGeneration #"+str(i))
        print("Population:\t",generation.asList())
        snapshot = Population.initFromPopulation(generation)
        fitnessDict = {}
        for chromosome in generation.chromosomes():
            fitnessDict[chromosome] = context.fitnessFunction()(chromosome)
        print("Fitnesses:\t",fitnessDict)
        selected = selection(fitnessDict)
        print("Selections:\t",selected)
        #kill the least fit
        generation.chromosomes().remove(selected[2])
        generation.chromosomes().remove(selected[3])

        #mate the most fit
        children = crossover(selected[0], selected[1], rateOfMutation)
        generation.chromosomes().append(children[0])
        generation.chromosomes().append(children[1])

        print(generation.asList())
    return selected[0]

def testGeneration(rateOfMutation: float, targetHammingDistance: int, context: Context) -> Chromosome:
    i = 1
    print("\n\nGeneration #"+str(i))
    generation = Population.initFromPopulation(context.initialPopulation())
    print("Population:\t",generation.asList())
    snapshot = Population.initFromPopulation(generation)
    fitnessDict = {}
    for chromosome in generation.chromosomes():
        fitnessDict[chromosome] = context.fitnessFunction()(chromosome)
    print("Fitnesses:\t",fitnessDict)
    selected = selection(fitnessDict)
    print("Selections:\t",selected)
    #kill the least fit
    generation.chromosomes().remove(selected[2])
    generation.chromosomes().remove(selected[3])

    #mate the most fit
    children = crossover(selected[0], selected[1], rateOfMutation)
    generation.chromosomes().append(children[0])
    generation.chromosomes().append(children[1])
    print(generation.asList())
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
