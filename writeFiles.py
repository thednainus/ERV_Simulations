class writeInfos(object):

    def __init__(self, model, numSimulation, virusMutationRate, hostMutationRate, rateNewCopy):

        self.model = model
        self.numSimulation = numSimulation
        self.virusMutationRate = virusMutationRate
        self.hostMutationRate = hostMutationRate
        self.rateNewCopy = rateNewCopy

    # Function that writes newick format trees to files in the computer
    def writeFile_Newick(self, prefix, newick):

        OutFileName = prefix+self.model+'_VirusMutationRate%s_RateNewCopy%s.tre' %(self.virusMutationRate, self.rateNewCopy)
        OutFile = open(OutFileName,'a')

        OutFile.write('# Simulation Number: %s, Virus Mutation Rate: %s, Rate New Copy: %s\n' %(self.numSimulation, self.virusMutationRate, self.rateNewCopy))
        OutFile.write(newick+'\n')

        OutFile.close()

    # Function that writes some information about the simulation into files in the computer
    def write_info(self, header, timeFinal):

        OutFileName = 'Info_'+self.model+'_VirusMutationRate%s_RateNewCopy%s.txt' %(self.virusMutationRate, self.rateNewCopy)
        OutFile = open(OutFileName,'a')

        if self.numSimulation == 0:
            OutFile.write(header)

        OutFile.write('%s\t%s\t%s\t%s\t%s\n' %(self.model, self.numSimulation, self.virusMutationRate, self.rateNewCopy, timeFinal))
        OutFile.close()