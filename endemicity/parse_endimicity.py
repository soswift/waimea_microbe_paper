# CODE CONVERTED FROM JUPYTER NOTEBOOK

from Bio import Phylo
import io
from collections import defaultdict
import numpy as np

# TODO: for finding nearest neighbors, consider computing the distance for each sequence
# to refernce points and cosider points based on on their distances to the closes reference points


default_annot = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

placement_data  = eval("".join(open("waimea.complete_placement.json").readlines()))
f = io.StringIO(placement_data['tree'])
tree = Phylo.read(f, 'newick')

annotations = dict([(x.strip().split("\t")[0], [y.strip() for y in x.rstrip().split("\t")[1].split(";")])
                       for x in open("ref/99_otu_taxonomy.txt").readlines()])





placements_per_branch = {}
placements_per_query = {}
for p in placement_data['placements']:
    if str(p['p'][0][0]) not in placements_per_branch:
        placements_per_branch[str(p['p'][0][0])]=[]
    placements_per_branch[str(p['p'][0][0])].append(p['nm'][0][0])

    placements_per_query[p['nm'][0][0]] = p['p'][0]






# Not used
# def find_common_annotations_from_terminals(clade):
#     terminals = clade.get_terminals()
#     matrix = []
#     for term in terminals:
#         matrix.append(annotations[term.name])
#
#     matrix = np.array(matrix)
#     non_unique = np.where((matrix[1:,:] !=  matrix[:-1, :]).sum(axis=0) > 0)[0]
#     if len(non_unique) > 0:
#         last_matchin_pos = max(non_unique)-1
#         valid_annot = False
#     else:
#         last_matchin_pos = -1
#
#     # print(matrix)
#     while True:
#         # print(last_matchin_pos)
#         if list(set(matrix[:, last_matchin_pos]))[0].split("__")[1] != '':
#             break
#     return matrix[0][last_matchin_pos]


def find_common_annotations_from_children(clade):

    if clade.is_terminal():
        return annotations[clade.name]

    matrix = []
    for c in clade:
        if c.is_terminal():
            matrix.append(annotations[c.name])
        else:
            matrix.append(c.annotation)
        
    matrix = np.array(matrix)
    #print("matrix is ")
    #print(matrix)

    non_unique = np.where((matrix[1:,:] !=  matrix[:-1, :]).sum(axis=0) > 0)[0]
    
    if len(non_unique) > 0:
        last_matchin_pos = max(non_unique)-1
        valid_annot = False
    else:
        last_matchin_pos = 6
        
    # print(matrix)
    while True:
        # print(last_matchin_pos)
        if list(set(matrix[:, last_matchin_pos]))[0].split("__")[1] != '':
            break
        last_matchin_pos -= 1

    last_matchin_pos = max(0, last_matchin_pos)

    annot = matrix[0][:last_matchin_pos+1]
    return list(annot) + default_annot[len(annot):]

comment_to_clade = {}
post_order = tree.get_nonterminals(order="postorder")
i = 0
for clade in post_order:
    i += 1
    clade.annotation = find_common_annotations_from_children(clade)
    comment_to_clade[clade.comment] = clade
    if i % 5_000 == 0 :
        print(i)

for c in tree.get_terminals():
    c.annotation = annotations[c.name]
    comment_to_clade[c.comment] = c


print("Computing assignments per tax level")
#

#==============================================
assignments_per_tax_level = defaultdict(list)
for c in post_order:
    if c.comment in placements_per_branch:
        for item in c.annotation:
             assignments_per_tax_level[item].append(c.comment)

for c in terminals:
    if c.comment in placements_per_branch:
        for item in c.annotation:
             assignments_per_tax_level[item].append(c.comment)


def find_suspicious_annots():
    for c  in post_order:
        child_annots = []
        for x in c:
            try:
                child_annots.append(np.array(x.annotation))
            except:
                child_annots.append(np.array(annotations[x.name]))

        if sum(child_annots[0] == child_annots[1]) != 7:
            print(c.annotation)
            print(list(child_annots[0]))
            print(list(child_annots[1]))
            print("*" * 60)



def compute_distance(query_1, query_2):
    
    branch_1_comment, _, _, graft_loc_1, pendant_length_1 = placements_per_query[query_1]
    clade_branch_1 = comment_to_clade[str(branch_1_comment)]

    branch_2_comment, _, _, graft_loc_2, pendant_length_2 = placements_per_query[query_2]
    clade_branch_2 = comment_to_clade[str(branch_1_comment)]

    
    # if they are both on the same branch, the distance is
    # pendant_length_1 + abs(graft_loc_1 - graft_loc_2) + pendant_length_2
    if branch_1_comment == branch_1_comment:
        return pendant_length_1 +  abs(graft_loc_1 - graft_loc_2) + pendant_length_2
    
    # else if branch of query_1 is prent of branch of query_2
    # ((branch_length_1 - graft_loc_1) + pendant_length_1) + (graft_loc_2 + pendant_length_2) + branch_length of all intermediat branches
    elif clade_branch_1.is_parent_of(clade_branch_2):
        print("I am in 2")
        intermediate_branches_lenght = clade_branch_1.distance(clade_branch_2) - clade_branch_2.branch_length - clade_branch_1.branch_length
        return (clade_branch_1.branch_length - graft_loc_1) + pendant_length_1 + (graft_loc_2 + pendant_length_2) + intermediate_branches_lenght
        
    # else if branch of query_2 is prent of branch of query_1
    # ((branch_length_2 - graft_loc_2) + pendant_length_2) + (graft_loc_1 + pendant_length_1) + branch_length of all intermediat branches
    elif clade_branch_2.is_parent_of(clade_branch_2):
        print("I am in 3")
        intermediate_branches_lenght = clade_branch_2.distance(clade_branch_1) - clade_branch_2.branch_length -	clade_branch_1.branch_length
        return (clade_branch_2.branch_length - graft_loc_2) + pendant_length_2 + (graft_loc_1 + pendant_length_1) + intermediate_branches_lenght

    # the are on different clades.
    # pendant_length_1 + graft_loc_1 + distance_branch_1_common_ancestor +
    #   pendant_length_2 + graft_loc_2 + distance_branch_2_common_ancesto
    else:
        print("I am in 4")
        common_ancestor = tree.common_ancestor([clade_branch_1, clade_branch_2])
        distance_branch_1_common_ancestor = common_ancestor.distance(clade_branch_1) - clade_branch_1.branch_length
        distance_branch_2_common_ancestor = common_ancestor.distance(clade_branch_2) - clade_branch_2.branch_length
        return distance_branch_1_common_ancestor + graft_loc_1 + pendant_length_1 + distance_branch_2_common_ancestor + graft_loc_2 + pendant_length_2


def find_closest_sequences_matching_tax(query_1, show_nb_neighbors=5_000, consider_nb_neighbors=10_000):
    ## Fix so that we add branches to  assignments_per_tax_level only if they have hits
    
    # We find branches with placement in in the same taxonomic levels
    branch_comment = placements_per_query[query_1][0]
    branch_annotation = comment_to_clade[str(branch_comment)].annotation
    distances = {}
    done = False
    for lvl in range(len(branch_annotation)-1,-1,-1):
        # print(f"1 I am in level {lvl}")
        # import ipdb
        # ipdb.set_trace()
        if branch_annotation[lvl].split("__")[1] == '':
            continue
        print(f"Considering {branch_annotation[lvl]}")
        i = 0
        for branch in assignments_per_tax_level[branch_annotation[lvl]]:
            # if placement already in distances that means the branch was already
            # considered at a lower tax level. Move to next branch
            if placements_per_branch[branch][0] in distances:
                continue
            for placement in placements_per_branch[branch]:
                i+=1
                # print(f"\t\t3 I am comuting distance with placement {placement}")
                distances[placement] = compute_distance(query_1, placement)
                consider_nb_neighbors -= 1
                if consider_nb_neighbors == 0 :
                    done = True
                    break
            if done:
                break
        print(f"\tfound {i} placements at this level ")
        if done:
            break
    print(f"consider_nb_neighbors is {consider_nb_neighbors}")
    return sorted(distances.items(), key= lambda x: x[1])
                




