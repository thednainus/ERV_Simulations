import time
from ete2 import Tree
from ete2.parser import newick
from Tree_utilities import TreeUtil
from Random_Functions import ERV_Randoms
from writeFiles import writeInfos


def main():
    model = 'TransposonImmortal'
    n = 50  # maximum number of copies
    n_Sim = 100  # number of trees to simulate
    numHostGen = 1000000  # Total Number of host generations
    totalN = 1  # It counts the number of elements that is in the list
    totalTips = 1  # It counts that number os elements that are being generated
    totalTime = 0  # It will count total Time

    hostMutationRate = 1.2 * 10**-8
    virusMutationRate = [1.0 * 10**-4, 1.0 * 10**-5, 3.0 * 10**-5, 1.0 * 10**-6, 1.0 * 10**-7, 1.2 * 10**-8]
    virusGenRate = [0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007, 0.0008, 0.0009] # Rate for the chosen element to release a new copy in the genome

    tips = ['Seq0']

    randoms = ERV_Randoms(numHostGen, hostMutationRate, 0)

    for sim in range(n_Sim):

        for rate in virusGenRate:

            for virusRate in virusMutationRate:

                startTime1 = time.time()

                t = Tree()
                ut = Tree()

                randoms.rateNewCopy = rate
                utils = TreeUtil(model, sim, virusRate, rate) # Initializes treeUtil utility python class
                write2Files = writeInfos(model, sim, virusRate, hostMutationRate, rate) # Initializes writeInfos python class
                header = 'Model\tSimulation Number\tERV Mutation Rate\tRate of generating a new copy\tLast Host Generation\n' # Header that will be passed to write_info function in writeInfos python class

                while totalTime <= numHostGen:

                    waiting_time = randoms.ERV_Time_Transposon(tips)
                    totalTime += waiting_time

                    if totalTime > numHostGen:
                        break

                    copy1 = randoms.choose_seq(tips)
                    copy2 = 'Seq' + str(totalTips)  # This is the new copy being generated

                    # It will append the new elements that are being generated, so it can be randomly selected later, when replacement starts to occur.
                    tips.append(copy2)

                    # The master copy will accumulate mutations only according to the host mutation rate
                    # while the new copy being generated will have a branch length that is a composite rate between the virus and host mutation rates.
                    copy1_time = waiting_time * hostMutationRate
                    copy2_time = (virusRate) + (waiting_time * hostMutationRate)


                    if totalN == 1:  # This mean that in the list samples, there is Seq0 and it was appended Seq1 (because I will sum 1 to TotalN, only in the end of this if

                        # To create a tree in which branch length will be in substitutions per site
                        child1 = t.add_child(name = copy1, dist = copy1_time)
                        child2 = t.add_child(name = copy2, dist = copy2_time)

                        # to create the ultrametric trees. Branch lengths in the tree will be a reflection of time
                        u_child1 = ut.add_child(name = copy1, dist = waiting_time)
                        u_child2 = ut.add_child(name = copy2, dist = waiting_time)

                        totalN += 1
                        totalTips += 1

                    elif totalN == 2:

                        t = utils.add_new_seq(t, copy1, copy2, copy1_time, copy2_time)
                        ut = utils.add_new_seq(ut, copy1, copy2, waiting_time, waiting_time)

                        totalN += 1
                        totalTips += 1

                    # When I reach two new copies in my list (TotalN = 3) than, the new copy may be replaced by a new element
                    elif totalN > 2 and totalN < n:

                        seqToReplace = randoms.choose_SeqToReplace(copy1, copy2, tips, totalN, n)

                        # If the SeqToReplace is None than it means only a new branch will be added, and no replacement will happen
                        if seqToReplace == None:

                            t = utils.add_new_seq(t, copy1, copy2, copy1_time, copy2_time)
                            ut = utils.add_new_seq(ut, copy1, copy2, waiting_time, waiting_time)

                            totalN += 1

                        # But if a sequence is selected, than the sequence in SeqToReplace will be removed from the newick format
                        else:
                            tips.remove(seqToReplace)  # It only removes a copy from Tips list when seqToReplace is different of None

                            t = utils.add_new_seq(t, copy1, copy2, copy1_time, copy2_time)
                            t = utils.delete_seq(t, seqToReplace)

                            ut = utils.add_new_seq(ut, copy1, copy2, waiting_time, waiting_time)
                            ut = utils.delete_seq(ut, seqToReplace)

                        totalTips += 1

                    elif totalN >= n:

                        seqToReplace = randoms.choose_SeqToReplace(copy1, copy2, tips, totalN, n)

                        tips.remove(seqToReplace)  # It only removes the sequence from Tips list when seqToReplace is different of None

                        t = utils.add_new_seq(t, copy1, copy2, copy1_time, copy2_time)
                        t = utils.delete_seq(t, seqToReplace)

                        ut = utils.add_new_seq(ut, copy1, copy2, waiting_time, waiting_time)
                        ut = utils.delete_seq(ut, seqToReplace)

                        totalTips += 1

                if totalTime > numHostGen:
                    lastGen = totalTime - waiting_time
                    timeFinal = randoms.final_Time(totalTime, waiting_time)

                    t = utils.update_leaves(t, None, None, timeFinal)
                    ut = utils.update_leaves(ut, None, None, (numHostGen - lastGen)) # numHostGen - lastGen will give the
                                                                                    # time in generations

                newick.set_float_format('%0.16f')
                t_newick = t.write(format=5)
                ut_newick = ut.write(format=5)

                write2Files.writeFile_Newick("subs_per_site_", t_newick)
                write2Files.writeFile_Newick("ultrametric_", ut_newick)
                write2Files.write_info(header, lastGen)

                # Reinitialize variables for the next simulation
                totalN = 1
                totalTips = 1
                totalTime = 0
                tips = ['Seq0']

                print "Elapsed time of first simulation is: %.5f" % (time.time() - startTime1)

if __name__ == "__main__":
    main()