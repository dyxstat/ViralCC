#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scipy.sparse as scisp
from math import log,exp,sqrt
import logging
import igraph as ig
import leidenalg
from sklearn.metrics import silhouette_score
import os


# package logger
logger = logging.getLogger(__name__)


class ClusterBin:
    def __init__(self, path , viral_info , map_combine , random_seed = 42):
        '''
        perc: threshold of spurious contacts
        viral_info: viral information
        min_combine: integrative graph
        '''

        #information of the viral contigs
        self.path = path
        self.viral_name = []
        for i in range(len(viral_info)):
            temp = viral_info[i]
            self.viral_name.append(temp.name)
        
        
        ###########use leiden algorithm to do clustering#########
        self.map_combine = scisp.coo_matrix(map_combine)
        
        vcount = self.map_combine.shape[0]
        sources = self.map_combine.row
        targets = self.map_combine.col
        index = sources>targets
        sources = sources[index]
        targets = targets[index]
        edgelist = list(zip(sources, targets))
        g = ig.Graph(vcount, edgelist)

        SIL_score = []
        cluster_range = np.arange(2,50)
        for n_clusters in cluster_range:
            part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition , resolution_parameter = n_clusters , n_iterations= -1 , seed = random_seed )
            part = list(part)

            label_pred = np.ones(self.map_combine.shape[0])
            for ci in range(len(part)):
                for contig in part[ci]:
                    label_pred[contig] = ci
                    
            if len(set(label_pred)) == self.map_combine.shape[0]:
                break

            SIL_score.append(silhouette_score(self.map_combine.todense() , np.array(label_pred)))

        optimal = SIL_score.index(max(SIL_score))
        part = leidenalg.find_partition(g , leidenalg.RBConfigurationVertexPartition, resolution_parameter = cluster_range[optimal] , n_iterations= -1 , seed = random_seed)
        part = list(part)
        logger.info('the number of generated viral bins is {}'.format(len(part)))

        self.dist_cluster = {}
        for ci in range(len(part)):
            for id in part[ci]:
                self.dist_cluster[self.viral_name[id]] = 'group'+str(ci)

        self._write_cluster()
        
            

    def _write_cluster(self):
        ########create file for checkm################
        with open(os.path.join(self.path ,'cluster_viral_contig.txt'),'w') as out:
            for key , value in self.dist_cluster.items():
                out.write(str(key)+ '\t' +str(value))
                out.write('\n') 
