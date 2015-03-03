import random,math

class ERV_Randoms(object):

    def __init__(self, numHostGen, hostMutationRate, rateNewCopy):

        self.numHostGen = numHostGen
        self.hostMutationRate = hostMutationRate
        self.rateNewCopy = rateNewCopy
        random.seed()

    # Function that calculates the waiting time for the Transposon model
    def ERV_Time_Transposon(self, tip_list):

        randNum = random.uniform(0,1)
        waiting_time = -math.log(randNum)/(self.rateNewCopy*(len(tip_list)))

        return waiting_time

    # Function that calculates the waiting time for the Master model
    def ERV_Time_Master(self):

        randNum = random.uniform(0,1)
        waiting_time = -math.log(randNum)/self.rateNewCopy

        return waiting_time

    # Function that randomly choose which copy will generate the new copy
    def choose_seq(self, tip_list):

        seq = random.choice(tip_list)		# It will random choose which copy will generate the new copy

        while seq == None:
            seq = random.choice(tip_list)

        return seq

    # Function that randomly choose which sequence (ERV copy) in the phylogenetic tree will be replaced.
    def choose_SeqToReplace(self, seq1, seq2, tip_list, total_N, max_N):

        prob = float(total_N-1)/(max_N-1) # This probability increases as it increases the number of tips in the phylogenetic tree
        seqToReplace = None
        randN = random.random()

        if randN < prob: # if  random number is less than prob I replace the sequence in the phylogenetic tree
            seqToReplace = random.choice(tip_list) # Randomly chooses a Sequence in the list Tips to be replaced
            # The copy to be selected has to be different from seq1 and seq2. This makes easier to to simulate trees
            # because I will know that seq1 is more closely related to seq2
            while seqToReplace == seq1 or seqToReplace == seq2:
                seqToReplace = random.choice(tip_list)

        return seqToReplace

    # Function that will return the Final Time that the last element was generated
    def final_Time(self, totalTime, waiting_time):

        totalTime -= waiting_time # I subtract the value in waiting_time because in my main code, I first sum to check if totalTime will be bigger than number os host generations
        timeFinal = self.numHostGen - totalTime
        timeFinal = timeFinal * self.hostMutationRate # Final time in substitution per site

        return timeFinal