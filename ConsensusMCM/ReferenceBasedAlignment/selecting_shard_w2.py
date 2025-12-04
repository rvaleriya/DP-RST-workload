import numpy as np
import os
import pandas as pd
import ot
np.random.seed(2025)

for task in ['Colon/PC3','Colon/PC10','Lung/PC3','Lung/PC10']:
    means_list_files = sorted(os.listdir('{}/means'.format(task)))

    # Read means
    means =[]
    for file in means_list_files:
        df = pd.read_csv('{}/means/{}'.format(task,file))
        means.append(df.to_numpy()[:,1:])

    # Read proportions
    proportions_list_files = sorted(os.listdir('{}/proportions'.format(task)))
    pis=[]
    for file in proportions_list_files:
        df = pd.read_csv('{}/proportions/{}'.format(task,file))
        pi = df.to_numpy()[:,1]
        pi =pi/np.sum(pi)
        pis.append(pi)


    N = len(means)

    # Construct the posterior similarity matrix
    similarity_matrix = np.zeros((N,N))
    for i in range(N):
        for j in range(i+1,N):
            pii = pis[i]
            meani = means[i]
            pij = pis[j]
            meanj = means[j]

            #  Wasserstein between two mixing measures (not really since it seems that you marginalizes out the mixing measures)
            M= ot.dist(meani,meanj)
            w = ot.emd2(pii,pij,M)
            similarity_matrix[i, j] = w
            similarity_matrix[j, i] = w
    # Average distance
    average_distance = np.mean(similarity_matrix,axis=0)

    # select the shard that have smallest average distance
    index = np.argmin(average_distance)
    print('Select shard {} for {}'.format(index+1,task))









