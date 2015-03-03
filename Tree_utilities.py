class TreeUtil(object):

    def __init__(self, model, numSimulation, virusMutationRate, rateNewCopy):

        self.model = model
        self.numSimulation = numSimulation
        self.virusMutationRate = virusMutationRate
        self.rateNewCopy = rateNewCopy

    # Function that updates the branch length for leaves in phylogenetic tree
    def update_leaves(self, tree, seq1, seq2, seq1_time):

        # updates all leaf branches with the exception of branches that were just added to the tree
        if seq1 != None and seq2 != None:
            for leaf in tree.iter_leaves():
                if leaf.name != seq1 and leaf.name != seq2:
                    leaf.dist = leaf.dist + seq1_time

        # if seq1 and seq2 equal "none" than all leaf branches will be updated
        else:
            for leaf in tree.iter_leaves():
                leaf.dist = leaf.dist + seq1_time

        return tree

    # Function that adds a child and a sister to the phylogenetic tree (as explained in the ete2 package)
    def add_new_seq(self, tree, seq1, seq2, seq1_time, seq2_time):

        # this for loop searches for the branch that has seq1 and add to this branch a new child and a new sister
        for leaf in tree.iter_leaves():
            if leaf.name == seq1:
                newChild = leaf.add_child(name = seq1, dist = seq1_time)
                newSister = newChild.add_sister(name = seq2, dist = seq2_time)
                break

        # After adding the new child and sister, the code needs to update the leaves of all other tips in the tree (all accumulates mutations according to the host mutation rate.
        tree = self.update_leaves(tree, seq1, seq2, seq1_time)

        return tree

    # Function that deletes a tip from the phylogenetic tree
    def delete_seq(self, tree, seqToReplace):

        for leaf in tree.iter_leaves():
            if leaf.name == seqToReplace:
                leaf.delete(prevent_nondicotomic=True, preserve_branch_length=True)

        return tree
