#### Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from numpy.random import rand
import seaborn as sb
import pandas as pd

from matplotlib import cm
from cycler import cycler


from sklearn.neighbors.kde import KernelDensity

from joblib import Parallel, cpu_count, delayed


import sys


##### Function to compute a realization of the process DEPRECIATED - USE .GO EXECUTIBLE
def transcription_realization(b0,g0,endt,kaps,phis_ids,phis_sc,site_par,gene_par,alpha,t0 = 0,see_real = True):
    '''Create a sample trajectory of a our stochastic regulation network (linear model) using the
    thinning gillespie method. Includes gene activity decay and linear interactions.
    Dependency: numpy, numpy.random.rand() and numpy.random.exponential()

    b0 - initial vector of site on/off
    g0 - initial transcript
    endt - float ending time
    kaps - list of genes for each site that ``bind" there. Doesn't have to be only one gene (which it
                would be for site splitting)
    phis_ids - list of sites (by index) that regulates it for each gene
    phis_sc - plus minus of sites above
    site_par - dataframe of lambdas,lambda hats,mus,nus,alphas. These DO NOT account for ``shared epigenetic"
                    parameters due to ``site splitting". (columns: 'lambda','lambda_hat','mu','nu')
    gene_par - dataframe of gammas, decay parameters (columns 'gamma','decay')

    '''
    t = [t0]
    jump_ts = [t0]
    real_jts = [t0]
    b = [np.array(b0)]
    g = [np.array(g0)]
    numg = len(g0)
    dec = gene_par.loc[:,'decay'].values
    gam = gene_par.loc[:,'gamma'].values
    lamhat = site_par.loc[:,'lambda_hat']
    stead_g_lis = [(gam + [np.dot(b[-1][phis_ids[i]],phis_sc[i]) for i in range(numg)])/dec]
    lamsum = sum(lamhat)
    epign = (site_par.loc[:,'lambda'].values*site_par.loc[:,'mu'].values)/(site_par.loc[:,'mu'].values+np.array(alpha)**site_par.loc[:,'nu'].values)
    while t[-1] < endt:
        tk = t[-1]
        gk = g[-1].copy()
        ###lam0 = Uniform upper bound on sum of all reaction propesenties
        steady_gs = (gam + [np.dot(b[-1][phis_ids[i]],phis_sc[i]) for i in range(numg)])/dec
        stead_g_lis = stead_g_lis + [steady_gs]
        gmx = np.maximum(g[-1],steady_gs)
        kapgs = [sum(li) for li in gmx[kaps]]
        lam0 = np.dot(epign,kapgs) + lamsum
        #### Find next possible jump time as exp(lam0)
        rrand = rand()
        Delt = np.log(1/rrand)/lam0#np.random.exponential(1/lam0)
        jump_ts += [tk+Delt]
        ##### Calculate the gs
        for tm in np.linspace(tk,jump_ts[-1],1):
            nxtg = np.exp(tk-tm)*(gk - steady_gs) + steady_gs
            g = g + [nxtg]
            t = t + [tm]
        ### Calculate all the intensities and cumulitive ones.
        fones = (1-b[-1])*(epign)*[sum(li) for li in g[-1][kaps]]
        ftwos = b[-1]*lamhat
        all_fs = np.concatenate([fones,ftwos])
        the_cu = np.array([sum(all_fs[:k]) for k in range(1,len(all_fs)+1)])
        the_cu2 = np.concatenate([the_cu,[lam0]])/lam0
        ### choose
        urand = rand()
        choice = min(np.argwhere(the_cu2>urand))[0]
        #now do the update
        new_b = b[-1].copy()
        if choice < len(new_b):#binding event
            new_b[choice] = new_b[choice] + 1
            real_jts += [tk+Delt]
        elif choice < 2*len(new_b):#unbinding event
            new_b[choice-len(new_b)] = new_b[choice-len(new_b)] - 1
            real_jts += [tk+Delt]
        b = b + [new_b]
    if see_real:
        return [np.array(b),np.array(g),t,jump_ts, stead_g_lis, real_jts]
    else:
        return [np.array(b),np.array(g)]



######## Long run realization to estimate distribution using ergodic thm, returns stats.DEPRECIATED - USE .GO EXECUTIBLE
def transcription_distribution(realization_params,eqtol = 0.001):
    '''Run a single realization and compute rolling average, rolling standard deviation, and point at wich
    the distribution seems to settle down into equilibrium. Return the g values past this point, the rolling
    statistics, and the realization in full'''
    #run the model realization
    b0,g0,endt,kaps,phis_ids,phis_sc,site_par,gene_par,alphas = realization_params
    b,g,t,jts,sgs,real_jumps = transcription_realization(b0,g0,endt,kaps,phis_ids,phis_sc,site_par,gene_par,alphas)
    arg = np.array(g)

    #get rolling average and standard deviation
    thav = [arg[0,:]]

    tpts = len(t)

    for i in range(1,tpts):
        thav += [np.average(arg[:i,:], axis = 0)]

    vars = [arg[0,:]]
    for i in range(1,tpts):
        vars += [np.std(arg[:i,:], axis = 0)]

    thav = np.array(thav)
    vars = np.array(vars)

    #find where rollimg average and standard deviation seem to stop moving.
    grads_mean  = np.array(np.diff(thav, axis = 0))
    grads_var = np.array(np.diff(vars, axis = 0 ))

    ateqafter = max([np.linspace(0,grads_mean.shape[0],50).astype(int)[np.where([(all(abs(grads_mean[i:,j]) < eqtol)) and all(abs(grads_var[i:,j]) < 2*eqtol) for i in np.linspace(0,grads_mean.shape[0],50).astype(int)])[0][0]] for j in range(grads_mean.shape[1])])

    equil_g = np.array(g)[ateqafter:,:]

    return [equil_g, thav, [thav-vars,thav+vars], t, np.array(b), arg,jts,sgs,real_jumps]


####### Compile long runs into one set of samplesDEPRECIATED - USE .GO EXECUTIBLE
def mc_equilib(realization_params,num_thru = 1):
    '''Run multiple realizations and compile into (1) all transcript values in equilibrium distribution put together,
    (2) the rolling statistics from each realization (list) and (3) each realization (list).
    Runs realizations in parallel. Will run num_thru*cpu_count simulations, doing cpu_count of them
    in parallel'''
    cpcnt = cpu_count() - 2 #count available cpus, leave 2 for other tasks
    all_runs = []
    #run realizations
    for i in range(num_thru):
        run1 = Parallel(n_jobs=cpcnt)(delayed(transcription_distribution)(realization_params) for i in range(cpcnt))
        all_runs += run1
    #combine samples from equilibrium
    full_on_equil_dist = np.concatenate([run[0] for run in all_runs])
    #make list of realizations
    list_of_rolling_stats = [[run[1],run[2],run[3]] for run in all_runs]
    list_of_realizations = [[run[4],run[5],run[3],run[6],run[7],run[8]] for run in all_runs]
    #
    return full_on_equil_dist,list_of_rolling_stats,list_of_realizations

###### Plot a single realization
def plot_realization(bb,gg,tt,tf_list = [],gene_nums = [],stdgs = [],gene_names = [],tf_names = [],lege = True):
    '''plot a realization. Can include the local equilibrium (i.e., the
    attracting value of g between jumps) by inputting stdgs.'''
    if len(tf_list) != 0:
        if len(gene_nums) == 0:
            gene_nums = list(set(range(gg.shape[1]))-set(tf_list))
        if len(gene_names) == 0:
            gene_names = ['Gene ' + str(i) for i in list(set(range(gg.shape[1]))-set(tf_list))]
        if len(tf_names) != len(tf_list):
            tf_names = ['TF: '+str(i) for i in tf_list]
        plt.rc('axes', prop_cycle=(cycler('color', [cm.tab20(i*(1/len(gene_nums))) for i in range(len(gene_nums))])))
        fig,ax = plt.subplots(1,2,constrained_layout=True,figsize = (7,5))
        ax[0].plot(tt,gg[:,tf_list])
        if len(stdgs) != 0:
            ax[0].step(tt,np.array(stdgs)[:,tf_list])
        if lege:
            ax[0].legend(tf_names)
        ax[0].set_title('TF Activity')
        # ax[1,0].step(tt,bb[:,tf_list])
        # ax[1,0].set_title('Bound State')
        # ax[1,0].legend([i + ' Promoter Site' for i in tf_names])
        ax[1].plot(tt,gg[:,gene_nums])
        if len(stdgs) !=0:
            ax[1].step(tt,np.array(stdgs)[:,gene_nums])
        if lege:
            ax[1].legend(gene_names)
        ax[1].set_title('Selected Gene Activity')
        # ax[1,1].step(tt,np.array(bb)[:,gene_nums])
        # ax[1,1].set_title('Bound State')
        # ax[1,1].legend([i + ' Promoter Site' for i in gene_nums])
    else:
        if len(gene_nums) == 0:
            gene_nums = range(gg.shape[1])
        if len(gene_names) == 0:
            gene_names = ['Gene ' + str(i) for i in range(gg.shape[1])]
        plt.rc('axes', prop_cycle=(cycler('color', [cm.tab20(i*(1/len(gene_nums))) for i in range(len(gene_nums))])))
        fig,ax = plt.subplots(1,1,constrained_layout=True,figsize = (7,5))
        ax.plot(tt,gg[:,gene_nums])
        if len(stdgs) !=0:
            ax.step(tt,np.array(stdgs)[:,gene_nums])
        if lege:
            ax.legend(gene_names)
        ax.set_title('Gene Activity')
        # ax[1].step(tt,np.array(bb)[:,gene_nums])
        # ax[1].set_title('Bound State')
        # ax[1].legend([i + ' Promoter Site' for i in gene_names])
    return fig

##### Plot equilibrium distribution from set of samples (2d array - each row a sampled vector, each column a gene)
def plot_dists(eqg,gene_names = []):
    '''plot estimated equilibrium distribution from samples'''

    if eqg.shape[1] % int(np.sqrt(eqg.shape[1])):
        num_y = int(np.sqrt(eqg.shape[1]))
        num_x = eqg.shape[1]//num_y + 1
    else:
        num_y = int(np.sqrt(eqg.shape[1]))
        num_x = num_y
    # plt_sz = round(eqg.shape[1]/12)
    wi = 14
    he = wi*(num_y/num_x)

    fig = plt.figure(figsize=(wi,he))

    # fig,axs = plt.subplots(3*plt_sz,4*plt_sz,constrained_layout=True,figsize=(14,7))
    fig.suptitle('Equilibrium Gene Expression Distribution')

    rectlis = [[i/num_x+0.03 ,0.03 + 0.95*((1-1/num_y)-j/num_y),0.8*(1/num_x),0.8*(1/num_y)] for i in range(num_x) for j in range(num_y) ]

    for i in range(eqg.shape[1]):
        ax = fig.add_axes(rectlis[i])
        sb.distplot(eqg[:,i],ax=ax)
        if len(gene_names):
            ax.set_title(gene_names[i])

    return fig

##### Plot equilibrium distribution from set of samples (dataframe indexed by GeneID)
def plot_dists2(eqg):
    '''plot estimated equilibrium distribution from samples'''

    if len(eqg) % int(np.sqrt(len(eqg))):
        num_y = int(np.sqrt(len(eqg)))
        num_x = len(eqg)//num_y + 1
    else:
        num_y = int(np.sqrt(len(eqg)))
        num_x = num_y
    # plt_sz = round(eqg.shape[1]/12)
    wi = 14
    he = wi*(num_y/num_x)

    fig = plt.figure(figsize=(wi,he))

    # fig,axs = plt.subplots(3*plt_sz,4*plt_sz,constrained_layout=True,figsize=(14,7))
    fig.suptitle('Equilibrium Gene Expression Distribution')

    rectlis = [[i/num_x+0.03 ,0.03 + 0.95*((1-1/num_y)-j/num_y),0.8*(1/num_x),0.8*(1/num_y)] for i in range(num_x) for j in range(num_y) ]

    for i in range(len(eqg)):
        ax = fig.add_axes(rectlis[i])
        sb.distplot(eqg.Samples[i],ax=ax)
        ax.set_title(eqg.index[i])

    return fig

##### Plot equilibrium distribution from set of samples (dataframe with estimated distribution)
def plot_dists3(eqg, rollstats = False):
    '''plot estimated equilibrium distribution from samples'''

    if len(eqg) % int(np.sqrt(len(eqg))):
        num_y = int(np.sqrt(len(eqg)))
        num_x = len(eqg)//num_y + 1
    else:
        num_y = int(np.sqrt(len(eqg)))
        num_x = num_y
    # plt_sz = round(eqg.shape[1]/12)
    wi = 10
    he = wi*(num_y/num_x)*0.8

    fig = plt.figure(figsize=(wi,he))

    # fig,axs = plt.subplots(3*plt_sz,4*plt_sz,constrained_layout=True,figsize=(14,7))
    fig.suptitle('Equilibrium Gene Expression Distribution')

    xcorner = 0.07 #spacing left of leftmost plots
    ycorner = 0.03 #spacing below bottom of bottom plots
    titlespacer = 0.95 #give this proportion to the plots, rest goes to title.
    plotspacer = titlespacer*0.85 #scale plots to fit subplot titles.

    rectlis = [[i/num_x+xcorner ,ycorner + titlespacer*((1-1/num_y)-j/num_y),plotspacer*(1/num_x),plotspacer*(1/num_y)] for i in range(num_x) for j in range(num_y) ]

    for i in range(len(eqg)):
        ax = fig.add_axes(rectlis[i])
        ax.plot(eqg.iloc[i].XDomain,eqg.iloc[i].DistEst)
        ax.set_title(eqg.iloc[i].GeneID)

    if rollstats:
        fig2 = plt.figure(figsize=(wi,he))
        fig2.suptitle('Rolling Statistics')
        for i in range(len(eqg)):
            ax = fig2.add_axes(rectlis[i])
            ax.plot(eqg.iloc[i].TimePts,eqg.iloc[i].RealMeans, label = "Mean")
            ax.plot(eqg.iloc[i].TimePts,eqg.iloc[i].RealVars, label = "Variance")
            ax.set_title(eqg.iloc[i].GeneID)
            ax.legend()

        return fig,fig2
    else:

        return fig

##### Plot a realization (genes only) rolling average and std dev.
def plot_rolling(means,first_stds,t,gene_names = []):
    '''plot rolling statistics from a realization. Only plots, doesn't compute the statistics

    means is n by k, where n is time points and k is genes
    first_stds is same but list of 2

    '''
    if means.shape[1] % int(np.sqrt(means.shape[1])):
        num_y = int(np.sqrt(means.shape[1]))
        num_x = means.shape[1]//num_y + 1
    else:
        num_y = int(np.sqrt(means.shape[1]))
        num_x = num_y

    wi = 14
    he = wi*(num_y/num_x)

    fig = plt.figure(constrained_layout=True,figsize=(wi,he))

    # fig2,axs2 = plt.subplots(3*plt_sz,4*plt_sz,constrained_layout=True,figsize=(14,7))
    fig.suptitle('Rolling Gene Expression Average')

    rectlis = [[i/num_x,(1-1/num_y)-j/num_y,0.83*(1/num_x),0.83*(1/num_y)] for i in range(num_x) for j in range(num_y) ]


    for i in range(means.shape[1]):
        ax = fig.add_axes(rectlis[i])
        ax.fill_between(t,first_stds[0][:,i],first_stds[1][:,i], alpha=0.5,color = 'red')
        ax.plot(t,means[:,i])
        if len(gene_names):
            ax.set_title(gene_names[i])

    return fig




############### Parameter estimation functions.

###compute log-liklihood
def compute_logLE(t1,t2,state,datapoints,alphas,network,params, bw = 1):
    '''

        Using monte carlo method, compute liklihood estimator
        L = sum_{datapoints i}(log(1/Nbw)(sum_{timesteps t}K((y_i^t - x_i)/bw)

        where K is some kernel - let's use a symmetric Gaussian first: K(z) = 1/(2pi)^(d/2) exp(-1/2 |z|^2)

        network is network structure parameters (not being estimated) (binder_inds,phis_ids,phis_sc)
        params is list of 2 DFs of parameters to be estimated
        alphas is list of alpha arrays for each datapoint (so match data!)

        state = [[b1,g1],...,[bn,gn]] are states of system, one for each datapoint.

    '''
    cpcnt = cpu_count() - 2
    binder_inds,phis_ids,phis_sc = network
    bb = [st[0] for st in state]
    gg = [st[1] for st in state]
    ###compute the realizations
    # realizations = [transcription_realization(bb[i],gg[i],t2,binder_inds,phis_ids,phis_sc,params[0],params[1],alphas[i],t0 = t1,see_real = False) for i in range(len(bb))]

    realizations = Parallel(n_jobs = cpcnt)(delayed(transcription_realization)(bb[i],gg[i],t2,binder_inds,phis_ids,phis_sc,params[0],params[1],alphas[i],t0 = t1,see_real = False) for i in range(len(bb)))
    ###then from those the estimate
    estim = sum([np.log((1/(bw*len(realizations[ii])))*sum([np.exp(-0.5*(np.linalg.norm(yy-datapoints[ii])/bw)**2) for yy in realizations[ii][1]])) for ii in range(len(datapoints))])
    return estim,[[yy[0][-1],yy[1][-1]] for yy in realizations]


def compute_avgdist(t1,t2,state,datapoints,alphas,network,params,delta = 1):
    '''

        Using monte carlo method, compute liklihood estimator
        L = sum_{datapoints}((1/N)log(sum_{timesteps}(indicator(realization in ball around data point))))

        network is network structure parameters (not being estimated) (binder_inds,phis_ids,phis_sc)
        params is list of 2 DFs of parameters to be estimated
        alphas is list of alpha arrays for each datapoint (so match data!)

        state = [[b1,g1],...,[bn,gn]] are states of system, one for each datapoint.

    '''
    cpcnt = cpu_count() - 2
    binder_inds,phis_ids,phis_sc = network
    bb = [st[0] for st in state]
    gg = [st[1] for st in state]
    ###compute the realizations
    # realizations = [transcription_realization(bb[i],gg[i],t2,binder_inds,phis_ids,phis_sc,params[0],params[1],alphas[i],t0 = t1,see_real = False) for i in range(len(bb))]

    realizations = Parallel(n_jobs = cpcnt)(delayed(transcription_realization)(bb[i],gg[i],t2,binder_inds,phis_ids,phis_sc,params[0],params[1],alphas[i],t0 = t1,see_real = False) for i in range(len(bb)))
    ###then from those the estimate
    estim = sum([(1/len(realizations[ii]))*sum([(np.linalg.norm(yy-datapoints[ii])) for yy in realizations[ii][1]]) for ii in range(len(datapoints))])/(len(datapoints))
    return estim,[[yy[0][-1],yy[1][-1]] for yy in realizations]


### simulated annealing parameter estimation.
def fit_params_mle(data,network,alphas,site_names,arrow_sites,min_g,heat_cycles = 1,cycle_length = 10,initial_check = 5):
    '''Using a method similar to simulated annealing, try to maximize the liklihood estimate

    alphas should be list of maps (pd series works) with keys in site_names

    '''
    sitep = pd.DataFrame(rand(len(network[0]),2), columns = ['lambda','lambda_hat'])
    epigenp = pd.DataFrame(rand(len(site_names),2), columns = ['mu','nu'], index = site_names)
    genep = pd.DataFrame(rand(len(data[0]),2), columns = ['gamma','decay'])
    genep.loc[:,'gamma'] = genep.loc[:,'gamma'] + min_g
    sitep.loc[:,'mu'] = [epigenp.loc[arrow_sites[i],'mu'] for i in sitep.index]
    sitep.loc[:,'nu'] = [epigenp.loc[arrow_sites[i],'nu'] for i in sitep.index]
    alphas_arr_li = [np.array([alp[asite] for asite in arrow_sites]) for alp in alphas]
    t = 0
    state = [[np.zeros(len(arrow_sites)),xx] for xx in data]
    mle,state = compute_logLE(t,initial_check,state,data,alphas_arr_li,network,[sitep,genep])
    mle_li = [mle]
    t += initial_check
    for i in range(heat_cycles):
        DT = 0.5
        ##### Perturb all of the parameters, save to a new DF
        site_new = pd.DataFrame(columns = ['lambda','lambda_hat','mu','nu'])
        epigen_new = pd.DataFrame(columns = ['mu','nu'], index = site_names)
        gene_new = pd.DataFrame(columns = ['gamma','decay'])
        gene_new.loc[:,'gamma'] = np.maximum(min_g,genep.loc[:,'gamma'] + (np.random.rand(len(genep.loc[:,'gamma']))-0.5))
        gene_new.loc[:,'decay'] = np.maximum(0,genep.loc[:,'decay'] + (np.random.rand(len(genep.loc[:,'decay']))-0.5))
        site_new.loc[:,'lambda'] = np.maximum(0,sitep.loc[:,'lambda']  + (np.random.rand(len(sitep.loc[:,'lambda']))-0.5))
        site_new.loc[:,'lambda_hat'] = np.maximum(0,sitep.loc[:,'lambda_hat']  + (np.random.rand(len(sitep.loc[:,'lambda_hat']))-0.5))
        epigen_new.loc[:,'mu'] = np.maximum(0,epigenp.loc[:,'mu'] + (np.random.rand(len(epigenp.loc[:,'mu']))-0.5))
        epigen_new.loc[:,'nu'] = epigenp.loc[:,'nu'] + (np.random.rand(len(epigenp.loc[:,'nu']))-0.5)
        site_new.loc[:,'mu'] = [epigenp.loc[arrow_sites[i],'mu'] for i in sitep.index]
        site_new.loc[:,'nu'] = [epigenp.loc[arrow_sites[i],'nu'] for i in sitep.index]
        for j in range(cycle_length):
            n_mle,state = compute_logLE(t,t+DT,state,data,alphas_arr_li,network,[sitep,genep])
            t += DT
            if n_mle > mle:
                sitep = site_new.copy()
                genep = gene_new.copy()
                epigenp = epigen.copy()
                mle = n_mle
                mle_li += [mle]
            site_new = pd.DataFrame(columns = ['lambda','lambda_hat','mu','nu'])
            epigen_new = pd.DataFrame(columns = ['mu','nu'], index = site_names)
            gene_new = pd.DataFrame(columns = ['gamma','decay'])
            gene_new.loc[:,'gamma'] = np.maximum(min_g,genep.loc[:,'gamma'] + (np.random.rand(len(genep.loc[:,'gamma']))-0.5))
            gene_new.loc[:,'decay'] = np.maximum(0,genep.loc[:,'decay'] + (np.random.rand(len(genep.loc[:,'decay']))-0.5))
            site_new.loc[:,'lambda'] = np.maximum(0,sitep.loc[:,'lambda']  + (np.random.rand(len(sitep.loc[:,'lambda']))-0.5))
            site_new.loc[:,'lambda_hat'] = np.maximum(0,sitep.loc[:,'lambda_hat']  + (np.random.rand(len(sitep.loc[:,'lambda_hat']))-0.5))
            epigen_new.loc[:,'mu'] = np.maximum(0,epigenp.loc[:,'mu'] + (np.random.rand(len(epigenp.loc[:,'mu']))-0.5))
            epigen_new.loc[:,'nu'] = epigenp.loc[:,'nu'] + (np.random.rand(len(epigenp.loc[:,'nu']))-0.5)
            site_new.loc[:,'mu'] = [epigenp.loc[arrow_sites[i],'mu'] for i in sitep.index]
            site_new.loc[:,'nu'] = [epigenp.loc[arrow_sites[i],'nu'] for i in sitep.index]
            DT = DT*1.5
    return [sitep,epigenp,genep,mle_li]
