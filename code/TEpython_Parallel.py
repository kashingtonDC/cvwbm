#!/usr/bin/env python
# Edom Moges (ESDL)

# ## This module contains the TE code translation from Matlab
# 
#       (MatLab version written by Laurel L. and modified by Dino B. Translation to Python by Edom M.)
# 
# The following functions are included in this module:
# 
# 1. Mutual information
# 
#     1. mutinfo_new(M, nbins) - Calculates mutual information I(x,y).
#     
#     
# 2. Tranfer entropy
# 
#     1. transen_new(M, lag, nbins) - Calculates transfer information - TE(x,y) x to y. x source M[:,0] and y the sink M[:,1].
#     
#     
# 3. Intermediate functions
# 
#     1. LagData_new - shifts a matrix so that it is rearranged to be ready for TE calculation as in Knutt et al., 2005
#     2. jointentropy_new(M, nbins) - Calculates the joint entropy H(x,y) 
#     3. jointentropy3_new(M, nbins) - Calculates the joint entropy for three variables H(x,y,z)
#     4. shuffle( M ) - shuffles the entries of the matrix M in time while keeping NaNs (blank data values) NaNs. So that, Monte Carlo is possible
#     5. transenshuffle_new(M, lag, nbins) - Calculates the transfer entropy for a shuffled time series that has already been lined up with LagData
#     
#     
# 4. Monte Carlo analysis of mutual information and transfer entropy
# 
#     1. mutinfo_crit_new( M, nbins, alpha, numiter) - Finds critical values of mutual information statistics that needs to be exceeded for statistical significance
#     2. transen_crit_new( M, lag, alpha, numiter, nbins) - Finds the critical value of the transfer entropy statistic that needs to be exceeded for statistical signficance
#     
# 
# 5. All in one code
#    RunNewTE2VarsPar(DataMatrix, LabelCell, shift, SinkNodes=None, SourceNodes=None,
                    # maxLag=100, minSamples=200, numShuffles = 500, sigLevel=0.05, numBins=[11,11,11],ncores=4) - runs all together in #parallel mode.


#get_ipython().run_line_magic('matplotlib', 'inline')
import pandas as pd
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import copy
import os
from joblib import Parallel, delayed
from tqdm import tqdm
np.random.seed(50)


# In[2]:


def checkMakeDir2(dirName):  # 
    result = dirName
    result2 = dirName*2
    return result, result2


# ### Mutual information

# In[8]:

def computeEntropy(M,nbin):
    N, binEdges1d=np.histogram(M,bins=nbin) #Which bin the data column is in
    p2 = N/sum(N)
    
    # Shanon entropy
    p2gt0 = p2[p2>0] # py
    log2p2gt0 = np.log2(p2gt0)
    H = (-sum(p2gt0*log2p2gt0))
    return H
    


def mutinfo_new(M, nbins):
    # Calculates mutual information
    # M is an array with two columns [ source, sink]
    # nbins list of number of bins in 1D, 2D and 3D, with three elements
    ths = 10e-4
    counts1, binEdges1=np.histogram(M[:,0],bins=nbins[1]) # Source Variable
    binEdges1[0] = binEdges1[0]-ths
    binEdges1[len(binEdges1)-1]=binEdges1[len(binEdges1)-1]+ths
    col1cat = np.digitize(M[:,0], binEdges1, right=False)
    
    counts2, binEdges2=np.histogram(M[:,1],bins=nbins[1]) # Sink Variable
    binEdges2[0] = binEdges2[0]-ths
    binEdges2[len(binEdges2)-1]=binEdges2[len(binEdges2)-1]+ths
    col2cat = np.digitize(M[:,1], binEdges2, right=False)  # which bin (ID) is the data located
    #print(col1cat)
    # convert 1D histogram to a 2D histogram
    jointentcat = (col1cat-1)*nbins[1]+col2cat #This classifies the joint entropy bin into a number between 1 and nbins^2. 0 is assigned to rows with misisng data.
    nbins_2 = nbins[1]**2
    N = np.bincount(jointentcat)[1:] # Number of datapoints within each joint entropy bin.
    p = N/sum(N); # Vector of probabilities
    
    # 1D probability/histogram
    N1, binEdges1d1=np.histogram(M[:,0],bins=nbins[0]) # Which bin the first data column is in
    N2, binEdges1d2=np.histogram(M[:,1],bins=nbins[0]) #Which bin the second data column is in
    
    p1 = N1/sum(N1)
    p2 = N2/sum(N2)
    
    
    # Shanon entropy
    pgt0 = p[p>0]    # px,y
    p1gt0 = p1[p1>0] # px
    p2gt0 = p2[p2>0] # py
    
    
    log2p2gt0 = np.log2(p2gt0)
    #Shannon entropy of the sink variable. Used to normalize mutual informaiton in the next line.
    Hy = (-sum(p2gt0*log2p2gt0))
    # Mutual information, in bits. Joint entropy is scaled to the number of bins in a single dimension.
    I = ( (-sum(p1gt0*np.log2(p1gt0)) - sum(p2gt0*log2p2gt0) ) + (sum(pgt0*np.log2(pgt0)))*np.log2(nbins[0])/np.log2(nbins[1]))/Hy 
    # double integral in the last component is done as a 1D. 
        
    #return nbins_2, jointentcat,p , sum(N), I, Hy
    return I
    
    


# ## Intermediate functions

# In[10]:


def LagData_new( M_unlagged, shift ):
    # LagData Shifts two time-series so that a matrix is generated that allows easy computation of Knutt et al 2005 based TE computation
    # M_unlagged is a matrix [X Y..n], where X and Y are column vectors of the
    # variables to be compared. shift is a row vector that says how much each
    # variable in M_unlagged is to be shifted by.
    
    nR,nC = np.shape(M_unlagged)
    maxShift = max(shift)
    minShift = min(shift)
    newlength = nR - maxShift + minShift
    M_lagged = np.nan*np.ones([newlength, nC]) #[source_lagged(1:n-lag), sink_unlagged(lag:n), sink_lagged(1:n-lag)]
    
    #@@@@@@@@@@@@@@######## Dino's verson uses shift of [0, 0, -lag ] for the shuffle case of transfer entropy (transenshuffle_new)
    for ii in range(np.shape(M_lagged)[1]):
        M_lagged[:,ii] = M_unlagged[(shift[ii]-minShift):(np.shape(M_unlagged)[0]-maxShift+shift[ii]), ii]
    
    return M_lagged

# Alternatively
#     lag = np.abs(shift[0])
#     M_lagged[:,0] = M_unlagged[0:(nR-lag), 0]
#     M_lagged[:,1] = M_unlagged[lag:(nR),1]
#     M_lagged[:,2] = M_unlagged[0:(nR-lag),2]
#    return M_lagged


# In[12]:


def jointentropy_new(M, nbins):
    # Calculates the joint entropy H(x,y)
    # M is two dimensional column matrix for which joint entropy is to be computed
    # H is the normalized joint entropy
    # nvalidpoints is the number of rows (samples) used to calculate the joint entropy
    
    ths = 10e-4
    counts1, binEdges1=np.histogram(M[:,0],bins=nbins) # Source Variable [ ]
    binEdges1[0] = binEdges1[0]-ths
    binEdges1[len(binEdges1)-1]=binEdges1[len(binEdges1)-1]+ths
    col1cat = np.digitize(M[:,0], binEdges1, right=False)
    
    counts2, binEdges2=np.histogram(M[:,1],bins=nbins) # Sink Variable
    binEdges2[0] = binEdges2[0]-ths
    binEdges2[len(binEdges2)-1]=binEdges2[len(binEdges2)-1]+ths
    col2cat = np.digitize(M[:,1], binEdges2, right=False)  # which bin (ID) is the data located
    #print(col1cat)
    # convert 1D histogram to a 2D histogram
    jointentcat = (col1cat-1)*nbins+col2cat #This classifies the joint entropy bin into a number between 1 and nbins^2. 0 is assigned to rows with misisng data.
    nbins_2 = nbins**2
    N = np.bincount(jointentcat)[1:] # Number of datapoints within each joint entropy bin.
    p = N/sum(N); # Vector of probabilities
    
    pgt0 = p[p>0]    # p(x,y)
    H = -sum(pgt0*np.log2(pgt0))
    nvalidpoints = sum(N)
    
    return H, nvalidpoints


# In[14]:


def jointentropy3_new(M, nbins):
    # Calculates the joint entropy for three variables H(x,y,z)
    # M is a three-column matrix that contains the input vectors of data.
    # nvalidpoints is the number of rows (samples) used to calculate the joint entropy
    
    ths = 10e-4
    counts1, binEdges1=np.histogram(M[:,0],bins=nbins) # Source Variable [ ]
    binEdges1[0] = binEdges1[0]-ths
    binEdges1[len(binEdges1)-1]=binEdges1[len(binEdges1)-1]+ths
    col1cat = np.digitize(M[:,0], binEdges1, right=False)
    
    counts2, binEdges2=np.histogram(M[:,1],bins=nbins) # Sink Variable
    binEdges2[0] = binEdges2[0]-ths
    binEdges2[len(binEdges2)-1]=binEdges2[len(binEdges2)-1]+ths
    col2cat = np.digitize(M[:,1], binEdges2, right=False)  # which bin (ID) is the data located
    
    counts3, binEdges3=np.histogram(M[:,2],bins=nbins) # Source Variable [ ]
    binEdges3[0] = binEdges3[0]-ths
    binEdges3[len(binEdges3)-1]=binEdges3[len(binEdges3)-1]+ths
    col3cat = np.digitize(M[:,2], binEdges3, right=False)
    
    # This classifies the joint entropy bin into a number between 1 and nbins^2. 0 is assigned to rows with misisng data.
    jointentcat = (col1cat-1)*nbins**2 + (col2cat-1)*nbins + col3cat 
    
    #print(np.asarray((jointentcat,col1cat,col2cat, col3cat)).T)
    
    nbins_3 = nbins**3
    N = np.bincount(jointentcat)[1:] # Number of datapoints within each joint entropy bin.
    sumN = sum(N)
    
    
    p = N/sumN  # Vector of probabilities
    pgt0 = p[p>0]
    H = -sum(pgt0*np.log2(pgt0))
    nvalidpoints = sumN
    
    
    return H, nvalidpoints 


# In[16]:


def shuffle( M ):
    # shuffles the entries of the matrix M in time while keeping NaNs (blank data values) NaNs.
    # M is the matrix where the columns are individual variables and the rows are entries in time
    
    Mss = np.ones(np.shape(M))*np.nan # Initialize
    
    for n in range(np.shape(M)[1]): # Columns are shuffled separately
        notnans = np.argwhere(~np.isnan(M[:,n]))
        R = np.random.rand(np.shape(notnans)[0],1) #np.random.rand(5,1)
        I = np.argsort(R,axis=0)
        #print(notnans[:,0])
        #print(notnans[I,0])
        #print('a',M[notnans[:,0],n])
        Mss[notnans[:,0],n] = M[notnans[I[:],0],n].reshape(np.shape(M)[0],)
        
    return  Mss
    


# ## Transfer entropy




# In[20]:


def transen_new2(M, shift, nbins): # with shift as an input different lags btween source and sink are possible
    # shift [-lag of source, 0, - lag of sink] # lag of sink usually being 1
    # Calculates transfer information
    # M is an array with two columns [ source, sink]
    # nbins list of number of bins in 1D, 2D and 3D, with three elements
    # lag is the time lag of interest. 
    # M4 is the lagged subset of data transfer entropy was run on. 
    
    M4 = LagData_new(np.column_stack((M, M[:,1])), shift) # source, sink, sink is input then 
    #  M4 becomes [source_lagged(1:n-lag), sink_unlagged(lag:n), sink_lagged(1:n-lag)]  => H(Xt-T, Yt, Yt-T)
    
    M4[np.argwhere(np.isnan(np.sum(M4,axis=1))), :] = np.nan # Reset rows with any NaN entry to NaN.
    M4short = M4[np.argwhere(~np.isnan(np.sum(M4,axis=1))),:] # Time series without NaN that will be passed on for shuffling.
    
    M1 = M4[:,(0,2)]  # [source_lagged(1:n-lag), sink_lagged(1:n-lag)]  =>H(Xt-T,Yt-T)
    M2 = M4[:,(1,2)] # [sink_unlagged(lag:n), sink_lagged(1:n-lag)]    =>H(Yt,Yt-T)
    
    #@@@@@@@@@@@@@@######## Dino uses M4[:,1]  to be predicted 
    M3 = M4[:,2]      # [sink_unlagged(lag:n)] to be predicted is used with DINO. BUT, need CORRECTION =>H(Yt) should be corrected to H(Yt-T) M[:,2]
                      # Knutt et al indicates lagged being used H(Yt-T). Thus, M4[:,2]
    # Now calculate the joint and marginal entropy components:
    T1, n_valid_pairs1 = jointentropy_new(M1,nbins[1])
    T2, n_valid_pairs2 = jointentropy_new(M2,nbins[1])
    
    # Entropy for the single predictor
    n3, valueatn = np.histogram(M3, nbins[0]) # results in count [n3] and the corresponding value
    n3gt0 = n3[n3>0]
    sumn3gt0 = sum(n3gt0)
    T3 = -sum((n3gt0/sumn3gt0)*(np.log2(n3gt0/sumn3gt0))) # Nonnormalized Shannon entropy of variable Y
    
    # Three variable entropy
    T4, n_valid_pairs4 = jointentropy3_new(M4,nbins[2])
    
    Tn = T3 # This is the Shannon entropy of Y, used to normalize the value of transfer entropy obtained below.
    
    log2nbins1 = np.log2(nbins[0])
    log2nbins2 = np.log2(nbins[1])
    log2nbins3 = np.log2(nbins[2])
    log2nbins1_2 = log2nbins1/log2nbins2
    log2nbins1_3 = log2nbins1/log2nbins3
    T1 = T1*log2nbins1_2
    T2 = T2*log2nbins1_2
    T4 = T4*log2nbins1_3
    
    T = (T1+T2-T3-T4)/Tn # Knuth formulation of transfer entropy
    
    N = min([n_valid_pairs1, n_valid_pairs2, n_valid_pairs4]) # Number of valid matched pairs used in the calculation
    
    return T, N, M4short
    

def transenshuffle_new2(M, shift, nbins):
    
    # Calculates the transfer entropy for a shuffled time series that has already been lined up with LagData
    
    # Calculates the transfer entropy of X>Y, the amount by which knowledge
    #   of variable X at a time lag reduces the uncertainty in variable Y. M =
    #   [X Y], and lag is the time lag of interest. nbins is the number of bins
    #   used to discretize the probability distributions.
    
    
    Minput = shuffle(M[:,(0,1)])
    T, N,c = transen_new2(Minput, shift, nbins)# use it but not understood why [0 0 -lag] is used instead of [-lag 0 -lag]
   
    
    return T


# ## Critical values of Mutual information and Transfer entropy

def mutualshuffle_new2(M,nbins):
    Mss = shuffle(M)
    MIss = mutinfo_new(Mss,nbins)
    return MIss


def mutinfo_crit_newPar( M, nbins, alpha, numiter,ncores):
    
    MIss = Parallel(n_jobs=ncores)(delayed(mutualshuffle_new2)(M, nbins) for ii in range(numiter))
    MIss = np.sort(MIss)
    MIcrit = MIss[round((1-alpha)*numiter)] # develop a histogram and peak the 95% quantile significance level with alpha = 0.05
    
    return MIcrit


def transen_crit_new2( M, shift, alpha, numiter, nbins,ncores):
    
    # Finds the critical value of the transfer entropy statistic
    # that needs to be exceeded for statistical signficance.
    # M = matrix of unshifted variables, e.g., [X Y] for calculating the X>Y transfer entropy. 
    # lag = time lag. 
    # alpha = significance level. 
    # numiter = number of Monte Carlo shufflings to perform. 
    # nbins = number of bins to use to discretize the probability distributions.
    
    # Serial
#     Tss = np.ones([numiter])*np.nan # Initializing shuffled transfer entropy table
#     #print(Tss)
    
#     for ii in range(numiter):
#         Tss[ii] = transenshuffle_new2(M, shift, nbins) # Calculates TE for each Monte Carlo Shuffling
    # parallel
    Tss = Parallel(n_jobs=ncores)(delayed(transenshuffle_new2)(M, shift, nbins) for ii in range(numiter))
    ## print(Tss)
    
    Tss = np.sort(Tss)
    Tcrit = Tss[round((1-alpha)*numiter)] # develop a histogram and peaks the 1-aplpha (95%) quantile significance level with alpha (= 0.05)
    
    return Tcrit




# ## Parallel TE & I calculater


# number of monteCarlo shuffle - kills the time - going from 100 to 1000 very time consuming. Parallel!!
# maxLag also takes a lot of time. Number of lag considered. 3*365
# number of source variables -- 20
def RunNewTE2VarsPar(DataMatrix, LabelCell, shift, SinkNodes=None, SourceNodes=None,
                     maxLag=100, minSamples=200, numShuffles = 500, sigLevel=0.05, numBins=[11,11,11],ncores=4):
    # computes TE assumes a data matrix with time in first columns and vars on others
         
    # Inputs
    # DataMatrix - data matrix with time in the first column
    # LabelCell - variable name of each data matrix entry
    # shift - shift for sink variable (-1,0,-1) for self optimality
    # Source_nodes - array of column indices for source variables [2]
    # Sink_nodes - array of column of indices for sink variales [3:end]
    # resultsDir - directory for results ./Results/
    # maxLag - maximum lag (3*365) 3 years
    # minSamples - minimum number of valid samples for TE (suggestion 200)
    # numShuffles - number of MonteCarlo shuffle iterations (suggestion 500)
    # sigLevel - significance level (suggested 0.05)
    # numBins - number of bins to use in 1, 2, and 3 dimensions default [11,11,11]
    # ncores - pc cores to use 
    # Outputs
    # Imat - mutual information
    # Icritmat - significance threshold
    # Tfirstmat - first T > Tcrit
    # Tbiggestmat - Tmax for T > Tcrit
    # Tcube_store - all T for all sink, source, lag combinations
    # Tcritcube_store - all Tcrits for all sink, source, lag combinations
    
    if DataMatrix.size == 0:
        return 'no dataMatrix'
    
    if LabelCell.size == 0:
        return 'no variable names'
    
    if SourceNodes is None:
        SourceNodes = np.arange(2,np.shape(DataMatrix)[1])
        
    if SinkNodes is None:
        SinkNodes = np.array([1])
        
    nSources = len(SourceNodes)
    nSinks = len(SinkNodes)
    
    # Start clock
    # print('Beginning 2-variable analysis (serial) ...')
    # Tot = tic
    # print(SourceNodes,SinkNodes)
    # =========================================
    ## Shrink input matrices to include only variables that are used
    # now the order is time, sinks, sources
    
    #@@@@@@@@@@@@@@@@@@@@@
    # from Pd to np.array
    dataMat = np.column_stack((DataMatrix[:,0], DataMatrix[:,SinkNodes], DataMatrix[:,SourceNodes])) # date, sink, sources
    labCell = np.r_[np.array([LabelCell[0]]), np.array(LabelCell[SinkNodes]), np.array(LabelCell[SourceNodes])]
              #np.r_[np.array(LabelCell[0]), np.array(LabelCell[1]), np.array(LabelCell[[2,3,4]])]
    #Or labCell = np.column_stack((LabelCell[:,0], LabelCell[:,SinkNodes], LabelCell[:,SourceNodes]))
     
    
    del DataMatrix # or set it to empty DataMatrix = []
    del LabelCell
    
    # =============================================
    # Initialize output matrices
    # mutual information between sources and sinks
    # the sink is daily mean Q, and all pairwise interactions are evaluated
    
    Imat = np.ones([nSinks,nSources])*np.nan # row value = # sink vars, col values = # source vars;

    # significance threshold

    Icritmat = copy.deepcopy(Imat)

    # first T > Tcrit
    Tfirstmat = copy.deepcopy(Imat)


    # Tmax for T > Tcrit
    Tbiggestmat = copy.deepcopy(Imat)

    # All T for all sink, source, lag combinations
    Tcube_store = np.ones([nSinks,nSources,maxLag])*np.nan

    # All Tcrits for all sink, source, lag combinations
    Tcritcube_store = copy.deepcopy(Tcube_store)
    
    
    # =============================================
    # LOOP OVER ALL PAIRS OF SOURCE AND SINK VARIABLES TO CALCULATE MI and TE
    
    for mySinkIter in range(nSinks):  # loop over Sink nodes (information receivers) [ 0]
        
        mySinkNum = SinkNodes[mySinkIter]
        mySinkInd = 1 + mySinkIter  # exclude time
        
        
        # extract sub-matrices for the ease of computation
        Ivec = Imat[mySinkIter,:]
        Icritvec = Icritmat[mySinkIter,:]
        Tfirstvec = Tfirstmat[mySinkIter,:]
        Tbiggestvec = Tbiggestmat[mySinkIter,:]
        Tmat_store = np.reshape(Tcube_store[mySinkIter,:,:],[nSources,maxLag])
        Tcritmat_store = np.reshape(Tcritcube_store[mySinkIter,:,:], [nSources,maxLag])
        sinkName = labCell[mySinkInd] # Text name of the Sink variable
        MmySink = dataMat[:,mySinkInd] # Select the sink variable to run
        
        print('sink node = ', mySinkIter)
        
        for mySourceIter in tqdm(range(nSources)): # Loop over the source nodes 
            #print(mySourceIter)
            
            mySourceNum = SourceNodes[mySourceIter]
            mySourceInd = 1 + nSinks + mySourceIter
            Mmysource = dataMat[:,mySourceInd] # Select source variables
            sourceName = labCell[mySourceInd]  # Name of the source variable
            print('Source node ', mySourceNum-1, sourceName, ':=>',  'Sink node ', mySinkNum, sinkName)
            #print('Lag ', 'Sink', 'Source')
            
            M = np.column_stack((Mmysource, MmySink)) # Source followed by Sink
            M = M.astype('float')
            #print(M.shape)
            # MUTUAL INFORMATION
            I = mutinfo_new(M,numBins) # computes mutual information
            Ivec[mySourceIter] = I     # save it in a matrix
            Icrit = mutinfo_crit_newPar(M=M, alpha=sigLevel, nbins=numBins,numiter = numShuffles,ncores=ncores)
            Icritvec[mySourceIter] = Icrit
            
            
            # TRANSFER ENTROPY
            T = np.ones([maxLag])*np.nan # intialize the TE vector over the range of lags examined
            Tcrit = copy.deepcopy(T) # Initialize the vector of the critical TE
                        
            #for lag in tqdm(range(maxLag)): #[0 to 364] in a year i.e., no lag day
            for lag in range(maxLag): #[0 to 364] in a year i.e., no lag day
                t, N, Mshort = transen_new2(M=M, shift=[-lag,shift[1],shift[2]], nbins=numBins) # Computes TE for at a given lag of 'lag'
                #print(Mshort, type(Mshort),Mshort.shape)
                Mshort = Mshort.reshape(Mshort.shape[0],Mshort.shape[2])
                if N >= minSamples: # enough length to compute TE
                    T[lag] = t      # save TE computed
                    Tcrit[lag] = transen_crit_new2(M=Mshort, shift=[-lag,shift[1],shift[2]], alpha= sigLevel,nbins=numBins,numiter=numShuffles,ncores=ncores) # TE critical
                #print(lag, mySinkIter, mySourceIter)    
            
            # Save the first and biggest value of T over the significance threshold
            
            TgTcrit = np.argwhere(T >= Tcrit)  # np.argwhere(np.array([5,6,9,18]) > np.array([3,9,2,9]))
            
            if any(TgTcrit):
                Tfirstvec[mySourceIter] = T[TgTcrit[0,0]]
                Tbiggestvec[mySourceIter] = max(T[TgTcrit[:,0]]) # @@@@@ Should be T-Tcrit biggest!!!!!!
            
            #print(Tcrit.shape, T.shape, Tcritcube_store.shape)
            Tmat_store[mySourceIter,:] = T
            Tcritmat_store[mySourceIter,:] = Tcrit
            
           
        # replace column vectors from source iterations into matrices
        Imat[mySinkIter, :] = Ivec
        Icritmat[mySinkIter, :] = Icritvec
        Tfirstmat[mySinkIter,:] = Tfirstvec
        Tbiggestmat[mySinkIter,:] = Tbiggestvec
        Tcube_store[mySinkIter,:,:] = Tmat_store
        Tcritcube_store[mySinkIter,:,:] = Tcritmat_store
              
    # save results (modify to save just relevant variables)
    # save([resultsDir 'TE_analysis_workspace.mat'], '-v7.3');

    # Stop clock
    #print('Finished 2-variable analysis (Parallel)!');    
        
       
    return  Imat, Icritmat, Tfirstmat,  Tbiggestmat, Tcube_store,  Tcritcube_store # | sink | source | lag |

