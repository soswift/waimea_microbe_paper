import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import sys



# process each index in parallel
# "python3 process_chunk.py split_data/{}.csv H2 >> outputs/{}.out"
# chosen_indices = ["connectance", "web asymmetry", "nestedness","NODF", "H2"]


pandas2ri.activate()

otus = pd.read_csv(sys.argv[1], index_col= "OTU_ID")

bipartite = importr('bipartite')
results = bipartite.networklevel(otus, sys.argv[2])
# results = bipartite.networklevel(otus, chose_indices)

print("{}\t{}\t{}".format(sys.argv[1], sys.argv[2], results)
