import scipy as sp
import scipy.stats as ss
import numpy as np
import networkx as nx
from convergence_ABM import *
from graph_generators import *
import networkx.algorithms.community as nx_comm

# calculate structural statistics
def get_struct_stats(G, node_comm):
    
    results = {}
    
    ## DEGREE METRICS

    # average degree
    deg_seq = list(dict(G.degree()).values())
    results["deg_avg"] = np.mean(deg_seq)
    # standard deviation degree
    results["deg_std"] = np.std(deg_seq)
    # assortativity degree
    results["deg_assort"] = nx.degree_pearson_correlation_coefficient(G)

    # COMMUNITY METRICS

    # average group size
    unique, comm_size = np.unique(node_comm, return_counts=True)
    results["comm_avg"] = len(unique)
    results["comm_avg"] = np.mean(comm_size)
    # standard deviation group size
    results["comm_std"] = np.std(comm_size)
    # to calculate the coefficient of variation divide the standard deviation by the mean
    #comm_stdn = comm_std/comm_avg

    # PARTITION METRICS

    # modularity
    results["modularity"] = get_comm_sets(G, node_comm)
    # % boundary nodes
    # % boundary edges
    pbe, pbn, bound_edges, bound_nodes = get_boundary_data(G,node_comm)
    results["pbe"] = pbe
    results["pbn"] = pbn
    
    # CLUSTERING METRICS

    # transitivity
    results["trans"] = nx.transitivity(G)
    
    # average clustering 
    clust = list(nx.clustering(G).values())
    results["clust_avg"] = np.mean(clust)
    # standard deviation clustering
    results["clust_std"] = np.std(clust)
    
    return results

# calculate structural statistics
def get_struct_stats_min(G):
    
    results = {}
    
    ## DEGREE METRICS

    # average degree
    deg_seq = list(dict(G.degree()).values())
    results["deg_avg"] = np.mean(deg_seq)
    # standard deviation degree
    results["deg_std"] = np.std(deg_seq)
    # assortativity degree
    results["deg_assort"] = nx.degree_pearson_correlation_coefficient(G)
   
    # CLUSTERING METRICS
    G1 = nx.Graph(G)
    # transitivity
    results["trans"] = nx.transitivity(G1)
    
    # average clustering 
    clust = list(nx.clustering(G1).values())
    results["clust_avg"] = np.mean(clust)
    # standard deviation clustering
    results["clust_std"] = np.std(clust)
    
    return results

def get_boundary_data(G,node_comm):
    # get boundary edges from edgelist
    bound_edges = []
    for (e1, e2) in G.edges():
        if node_comm[e1] != node_comm[e2]:
            bound_edges.append((e1, e2))
    # get boundary nodes from boundary edges
    bound_nodes = np.unique(bound_edges)
    # percentage boundary edges 
    pbe = len(bound_edges)/G.number_of_edges()
    # percentage boundary nodes
    pbn = len(bound_nodes)/G.number_of_nodes()
    # return boundary statistics
    return pbe, pbn, bound_edges, bound_nodes

# get community sets
def get_comm_sets(G, node_comm):
    # get group numbers
    unique, comm_size = np.unique(node_comm, return_counts=True)
    # get partition
    communities = []
    for c in unique:
        communities.append(set(np.where(node_comm == c)[0]))
    # return list
    return nx_comm.modularity(G, communities)

# get community assignment
def get_sbm_comms(node_comm):
    # get community size
    unique, comm_size = np.unique(node_comm, return_counts=True)
    comm_sbm = []
    # generate community assignment
    for i, c in enumerate(unique):
        comm_sbm = comm_sbm + [int(c)]*comm_size[i]
    # cast to array
    return np.array(comm_sbm)

# get probability matrix for SBM
def get_sbm_probs(G, node_comm):
    # get adjacency matrix
    A = nx.to_numpy_array(G)
    # get group sizes
    unique, comm_size = np.unique(node_comm, return_counts=True)
    # set group number
    k = len(comm_size)
    # initialise
    sbm_p = np.zeros([k, k])
    # calculate probabilities
    for i in range(k):
        for j in range(k):
            # get nodes belonging to the selected groups
            i_pos = np.where(node_comm==i)[0]
            j_pos = np.where(node_comm==j)[0]
            # count unique edges
            edges = A[np.ix_(i_pos,j_pos)].sum()/2
            if i == j:
                # calculate max edges
                edges_max = comm_size[i]*(comm_size[i] - 1)/2
                # calculate probability
                edges_p = edges/edges_max
                sbm_p[i][i] = edges_p
            else:
                # calculate max edges
                edges_max = comm_size[i]*comm_size[j]/2
                # calculate probability
                edges_p = edges/edges_max
                sbm_p[i][j] = edges_p            
                sbm_p[j][i] = edges_p
    return sbm_p

#simulate several runs updating the graph at each run
def simulate_bulk_update(config_param_c, rounds, get_model, G, node_comm):
    # get error points
    err_len = int(config_param_c["max_time"]/config_param_c["log_interval"])
    # initialise vectors
    errs = np.zeros([rounds,err_len])
    ms = np.zeros(rounds)
    # simulate rounds times
    for i in range(rounds):
        G1 = get_model(G, node_comm)
        config_param_c["graph"] = G1
        sim = EventDrivenNetworkModel(config_param_c)
        err, m = sim.run_simulation()
        errs[i,:] = err
        ms[i] = m
    # average results
    err = np.mean(errs, axis = 0)
    m = np.mean(ms)
    # return results
    return errs, ms, err, m

# simulate several runs with same graph
def simulate_bulk(config_param_c, rounds):
    # get error points
    err_len = int(config_param_c["max_time"]/config_param_c["log_interval"])
    # initialise vectors
    errs = np.zeros([rounds,err_len])
    ms = np.zeros(rounds)
    # simulate rounds times
    for i in range(rounds):
        sim = EventDrivenNetworkModel(config_param_c)
        err, m = sim.run_simulation()
        errs[i,:] = err
        ms[i] = m
    # average results
    err = np.mean(errs, axis = 0)
    m = np.mean(ms)
    # return results
    return errs, ms, err, m

# get null models for a given graph
def get_null_ER(G, node_comm):
    deg_seq = list(dict(G.degree()).values())
    deg_avg = np.mean(deg_seq)
    s = nx.number_of_nodes(G)
    return get_connected_erdos(s, deg_avg/s)    

def get_null_CM(G, node_comm):
    deg_seq = list(dict(G.degree()).values())
    return get_connected_config(deg_seq)

def get_null_SBM(G, node_comm):
    unique, comm_size = np.unique(node_comm, return_counts=True)
    sbm_p = get_sbm_probs(G, node_comm)
    return get_connected_sbm(comm_size, sbm_p)    

def get_null_models(G, node_comm):   
    # (1) Erdos-Renyi
    G_ER = get_null_ER(G, node_comm)   
    # (2) Configuration Model
    G_CM = get_null_CM(G, node_comm)    
    # (3) Stochastic Block Model
    G_SBM = get_null_SBM(G, node_comm) 
    return G_ER, G_CM, G_SBM

def get_deg(G):
    # return degree list
    return list(dict(G.degree()).values())

def get_deg_avg(G):
    # return degree list
    return np.mean(get_deg(G))
