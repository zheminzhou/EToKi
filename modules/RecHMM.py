import numpy as np, pandas as pd, sys, os, copy, argparse, re
from numba import jit
from time import time
import functools, datetime
from multiprocessing import Pool
try:
    xrange(2)
except :
    xrange = range
try:
    from configure import uopen
except :
    from .configure import uopen

def _iter_branch_measure(obj, arg) :
    return obj.iter_branch_measure(arg)

def _iter_viterbi(obj, arg) :
    return obj.viterbi(arg)    

@jit(nopython=True, fastmath=True)
def update_distant_transition(transition, emission, dist_transition, dist_transition_adj) :
    interval = dist_transition.shape[0]
    dist_transition[0] = transition

    saturate_id = 0
    ss = 0.
    for i in range(interval-1) :
        saturate_id = i
        t = np.dot(transition*emission[0], dist_transition[i])
        s = np.sum(t)/transition.shape[0]
        dist_transition[i+1] = t/s
        ss = ss + np.log(s)
        dist_transition_adj[i+1] = ss
        if np.sum(np.abs(dist_transition[i+1] - dist_transition[i])) <= 1e-10 :
            dist_transition[i+2:] = dist_transition[i+1]
            dist_transition_adj[i+1:] = dist_transition_adj[i+1] + np.arange(interval-i-1)*np.log(s)
            break
    return dist_transition, dist_transition_adj, saturate_id


class recHMM(object) :
    def __init__(self, mode=1) :
        self.max_iteration = 200
        self.n_base = None
        self.mode = ['legacy', 'hybrid', 'intra', 'both'][mode]
        self.n_a, self.n_b = [[2, 2], [4, 3], [2, 3], [3, 3]][mode]

    def fit(self, mutations, branches=None, sequences=None, missing=[], categories=None, init=None, cool_down=5) :
        self.observations = self.prepare_branches(mutations, sequences, missing)
        self.branches = branches if branches is not None else np.arange(len(self.observations)).astype(str)
        self.categories = { 'noRec':{} }
        for c, assigns in categories.items() :
            self.categories[c] = np.zeros(shape=[len(self.branches)], dtype=int)
            if assigns.get('*', -1) == 0 :
                self.categories[c][:] = np.arange(len(self.branches))
            else :
                for i, n in enumerate(self.branches) :
                    self.categories[c][i] = categories[c].get(n, 0)

        models = self.initiate(self.observations, init=init)
        return self.BaumWelch(models, self.max_iteration, cool_down=cool_down)

    def save(self, fout):
        import json
        model = copy.deepcopy(self.model)
        model['model'] = dict(
            n_a = self.n_a,
            n_b = self.n_b,
            n_base = self.n_base,
        )
        for k in ('v', 'v2', 'R', 'EventFreq', 'theta', 'delta') :
            model[k] = model[k].tolist()
        if 'posterior' in model :
            model['posterior'] = { k:v.tolist() for k,v in model['posterior'].items() }
        if 'categories' in model:
            for k in ('R/theta', 'nu', 'delta') :
                model['categories'][k] = model['categories'][k].tolist()
        json.dump(model, fout)
        return self.model

    def load(self, fin) :
        import json
        model = json.load(fin)
        for k in ('v', 'v2', 'R', 'EventFreq', 'theta', 'delta') :
            model[k] = np.array(model[k])
        if 'posterior' in model :
            model['posterior'] = { k:np.array(v) for k,v in model['posterior'].items() }
        if 'categories' in model:
            for k in ('R/theta', 'nu', 'delta') :
                model['categories'][k] = np.array(model['categories'][k])
        self.__dict__.update(model.pop('model', {}))
        self.model = model
        return self.model


    def initiate(self, observations, init) :
        criteria = np.array(init.split(',')).astype(float)
        intervals = np.sort(np.concatenate([ np.diff(obs[obs.T[3] > 0, 5]) for observation in observations for obs in observation ]))
        criteria = (criteria * intervals.size).astype(int)
        cutoffs = np.unique(intervals[criteria])

        self.models = []
        for cutoff in cutoffs :
            mut_summary, rec = [], []
            for x, observation in enumerate(observations) :
                mut = []
                for obs in observation :
                    substitution = obs[obs.T[3] > 0]
                    if len(substitution) :
                        substitution = np.vstack([obs[-1, :-1].tolist()+[-(obs[-1, -1] - substitution[-1, -1])], substitution])
                        dist = np.diff(substitution.T[-1])
                        inRec = np.concatenate([[False], (dist <= cutoff)[1:], [False]])
                        edges = np.concatenate([[0], np.diff(inRec.astype(int))])
                        mutSites = substitution[~(inRec | np.concatenate([[False], inRec[:-1]])), 3]
                        edgeSites = substitution[np.where(edges != 0), 3]

                        recRegion = []
                        for s, e in np.vstack([np.where(edges > 0), np.where(edges < 0)]).T :
                            r = substitution[np.arange(s, e+1)]
                            recRegion.append([1, r[-1, -1] - r[0, -1], r.shape[0]-1, (np.sum((r[1:-1, 3]-1)/r[1:-1, 3]) + np.sum((r.T[3]-1)/r.T[3]))/2])
                        if len(recRegion) :
                            recRegion = np.array(recRegion)
                            rec.append(recRegion)
                            mut.append([ substitution[-1, -1] - substitution[0, -1] - np.sum(recRegion.T[1]), substitution.shape[0]-1 - np.sum(recRegion.T[2]), 
                                                 np.sum((substitution[1:, 3] - 1)/substitution[1:, 3])-np.sum(recRegion.T[3]), recRegion.shape[0] ])
                        else :
                            mut.append([ substitution[-1, -1] - substitution[0, -1], substitution.shape[0]-1, 
                                                 np.sum((substitution[1:, 3] - 1)/substitution[1:, 3]), 0 ])
                    else :
                        mut.append([ obs[-1, -1] - obs[0, -1], 0, 0, 0 ])
                mut_summary.append(np.sum(mut, 0))
            
            mut_summary = np.array(mut_summary)
            #mut_summary.T[2] -= mut_summary.T[1]
            
            EventFreq = np.sum(mut_summary.T[[1, 3]]/mut_summary.T[0], 0)

            n_br, (bases, muts, homos, recs) = mut_summary.shape[0], np.sum(mut_summary, 0)
            model = dict(theta = np.array([ muts/np.sum(EventFreq)/self.n_base for id in np.unique(self.categories['R/theta']) ]),
                         h     = [ max(0.01, min(0.95, homos/muts)), max(0.01, min(0.95, 3.*homos/muts)) ],
                         probability = -1e300,
                         diff        = 1e300,
                         EventFreq   = EventFreq,
                         id          = len(self.models) + 1,
                         ite         = 0,
                         categories  = copy.deepcopy(self.categories)
                        )

            rec = np.vstack(rec)
            #rec.T[3] -= rec.T[2]

            if rec.shape[0] :
                rec2, rec = np.sum(rec[rec.T[3]/rec.T[2] > 1.5*model['h'][0]], 0), np.sum(rec[rec.T[3]/rec.T[2] <= 1.5*model['h'][0]], 0)
                rec[0], rec2[0] = max(rec[0], 1), max(rec2[0], 1)

                model['h'][1]  = max(model['h'][1], (rec2[3]+0.5)/(rec2[2]+1))
                model['delta'] = np.array([ max(0.00002, min(0.005, (rec[0]+rec2[0])/(rec2[1] + rec[1]))) for id in np.unique(self.categories['delta']) ])
                model['v']     = np.array([ np.max([0.001, np.min([0.7, (rec[2]+rec2[2])/(rec2[1] + rec[1])])]) for id in np.unique(self.categories['nu']) ])
                model['v2']    = np.array([ np.min([0.7, (rec2[2]+1)/(rec2[1]+1)]) for id in np.unique(self.categories['nu']) ])

                p = np.array([rec[0], np.sqrt(rec[0]*rec2[0]), rec2[0]])[:(self.n_a-1)]
                p /= np.sum(p)
                model['R'] = np.array([ recs/np.sum(EventFreq)/self.n_base * p for id in np.unique(self.categories['R/theta']) ])
                tot_event = model['theta'][0] + np.sum(model['R'][0])
                model['theta'] /= tot_event
                model['R'] /= tot_event
                self.screen_out('Initiate', model)
                self.models.append(model)
        return self.models


    def BaumWelch(self, models, max_iteration, cool_down=5) :
        n_model = len(models)
        for ite in xrange(max_iteration) :
            new_models = []
            self.model = models[0]

            for model in models:
                if 'diff' in model and model['diff'] < 0.001 :
                    new_models.append(model)
                else :
                    print('')
                    self.screen_out('Assess', model)
                    #t = time()
                    branch_params = self.update_branch_parameters(model)
                    branch_measures = self.get_branch_measures(branch_params, self.observations)
                    #print(time() - t)
                    prediction = self.estimation(model, branch_measures)
                    prediction['diff'] = -prediction['probability'] if not model['probability'] else prediction['probability'] - model['probability']
                    if prediction['diff'] > 0 :
                        prediction['ite'] = ite+1
                        self.screen_out('Update', prediction)
                        new_models.append(prediction)
                    else :
                        curr_model = copy.deepcopy(model)
                        curr_model['diff'] = prediction['diff']
                        self.screen_out('Freeze', curr_model)
                        new_models.append(curr_model)
                        if ite <= min(cool_down, 50) :
                            prediction['id'] = np.round(prediction['id'] + 0.01, 3)
                            prediction['ite'] = ite+1
                            prediction['diff'] = 1
                            new_models.append(prediction)

            new_models = sorted(new_models, key=lambda x:-x['probability'])
            if ite > 0 and ite % cool_down == 0 :
                while len(new_models) > 1 and n_model - len(new_models) < ite/cool_down :
                    if new_models[-1]['probability'] + new_models[-1]['diff']/2 >= new_models[-2]['probability'] + new_models[-2]['diff']/2 :
                        new_models[-1]['diff'] /= 2
                        new_models[-2]['diff'] /= 2
                        break
                    else :
                        self.screen_out('Delete', new_models[-1])
                        new_models = new_models[:-1]
                self.verify_model(new_models)
            models = new_models
        self.screen_out('Report', models[0])
        return models[0]

    def verify_model(self, models) :
        for model in models :
            if 'low_cov' not in model['categories'] :
                model['categories']['low_cov'] = {}
            for brId, (theta, v) in enumerate(zip(model['posterior']['theta'], model['posterior']['v'])) :
                if v[0] > .05 * theta[0] and v[1]*theta[0] < theta[1]*v[0] and np.sum(model['categories']['nu'] == model['categories']['nu'][brId]) > 2 :
                    print('Model {0}: Branch {1} is too divergent for Rec identification. Move to new category. '.format(model['id'], self.branches[brId]))
                    model['categories']['nu'][brId] = new_id = model['categories']['nu'][brId] + 1
                    if new_id >= model['v'].shape[0] :
                        model['v'] = np.concatenate([model['v'], [0.]])
                        model['v2'] = np.concatenate([model['v2'], [model['v2'][-1]]])
                    if model['v'][new_id] < theta[1]/theta[0]*1.1 :
                        model['v'][new_id] = theta[1]/theta[0]*1.1
                    if model['v2'][new_id] < theta[1]/theta[0]*0.5 :
                        model['v2'][new_id] = theta[1]/theta[0]*0.5
                    model['probability'] = -1e300
                    model['diff'] = 1e300
                elif model['posterior']['R'][brId, 0] * 2. < self.n_base :
                    tId = int(model['categories']['R/theta'][brId])
                    if tId not in model['categories']['low_cov'] :
                        if np.sum(model['categories']['R/theta']== tId) > 1:
                            model['categories']['low_cov'][tId] = int(model['theta'].size)
                            model['categories']['low_cov'][int(model['theta'].size)] = int(model['theta'].size)
                            model['R'] = np.vstack([model['R'], [model['R'][tId]]])
                            model['theta'] = np.concatenate([model['theta'], [model['theta'][tId]]])
                        else :
                            model['categories']['low_cov'][tId] = tId
                    if model['categories']['low_cov'][tId] != tId :
                        print('Model {0}: Branch {1} got too much recombined regions. Move to new category. '.format(model['id'], self.branches[brId]))
                        model['categories']['R/theta'][brId] = model['categories']['low_cov'][tId]
                        model['probability'] = -1e300
                        model['diff'] = 1e300

    def screen_out(self, action, model) :
        if not verbose :
            return
        print('{2}\t{0} model {id}[{ite}] - BIC: {3:.8e} - EventFreq: {1:.3e}; theta: {theta[0]:.3f}; R: {5:.3f}; delta: {4:.3e};  Nu: {v[0]:.3e},{v2[0]:.3e}; h: {h[0]:.3f},{h[1]:.3f}'.format(action, np.sum(model['EventFreq']), str(datetime.datetime.now())[:19],  -2*model['probability'] + self.n_a*self.n_b*np.log(self.n_base*len(self.observations)), 1/model['delta'][0], np.sum(model['R'][0]), **model))
        sys.stdout.flush()

    def estimation(self, model, branch_measures) :
        #h = [1./n+(n-2.)/n*model['h'][0] for n in np.arange(1, self.n_b)]
        h = [ ((1-model['h'][0])*(model['h'][0]**n))**(1/(n+1)) for n in np.arange(self.n_b-1) ]

        probability, br = 0., []

        posterior = dict( theta=np.zeros([len(branch_measures), 2]),
                          h=np.zeros([len(branch_measures), 4]),
                          R=np.zeros([len(branch_measures), 4]),
                          delta=np.zeros([len(branch_measures), 3, 2]),
                          v=np.zeros([len(branch_measures), 2]),
                          v2=np.zeros([len(branch_measures), 2]),
                          probability=np.zeros(len(branch_measures)), )

        for id, measures in enumerate(branch_measures) :
            a, b = measures['a'], measures['b']
            if self.n_a == 2 and self.n_b > 2 :
                posterior['theta'][id, :] = [ np.sum(b[0]) + np.sum(b[1])*h[0]/2, np.sum(b[0, 1:]) + b[1, 1] ]
                posterior['h'][id] = [ np.sum(b[0, 1:]), np.sum(b[0, 2:] * [(n-1)/n for n in np.arange(2, self.n_b)]), np.sum(b[1, 1:]), np.sum(b[1, 2:] * [(n-1)/n for n in np.arange(2, self.n_b)]) ]
                posterior['v'][id] = [np.sum(b[1]), np.sum(b[1, 1:])]
                posterior['v2'][id] = [np.sum(b[1]), b[1, 2]]
            else :
                posterior['theta'][id, :] = [ np.sum(b[0]) + np.sum(b[2:3])*h[0]/2, np.sum(b[0, 1:]) + np.sum(b[2:3, 1]) ]
                posterior['h'][id] = [ np.sum(b[:2, 1:]), np.sum(np.sum(b[:2, 2:], 0) * [(n-1)/n for n in np.arange(2, self.n_b)]), np.sum(b[2:, 1:]), np.sum(np.sum(b[2:, 2:], 0) * [(n-1)/n for n in np.arange(2, self.n_b)]) ]
                posterior['v'][id] = [ np.sum(b[1]) + np.sum(b[3:])*h[0]/2, np.sum(b[1, 1:]) + np.sum(b[3:, 1]) ]
                posterior['v2'][id] = [ np.sum(b[2:]), np.sum(b[2:, 2:]) ]

            posterior['R'][id, 0], posterior['R'][id, 1:] = np.sum(a[0]), a[0, 1:]
            posterior['delta'][id] = np.vstack([np.sum(a[1:], 1), a.T[0, 1:]]).T
            posterior['probability'][id] = measures['probability']

        prediction = copy.deepcopy(model)
        prediction['posterior'] = posterior
        prediction['probability'] = np.sum(posterior['probability'])

        EventFreq = np.vstack([posterior['theta'].T[1]/posterior['theta'].T[0], posterior['R'][:, 1:].T/posterior['R'][:, 0]]).T
        prediction['EventFreq'] = np.sum(EventFreq, 1)
        prediction['h'] = [np.sum(posterior['h'].T[1])/np.sum(posterior['h'].T[0]), np.sum(posterior['h'].T[3])/np.sum(posterior['h'].T[2])]

        for id, theta in enumerate(prediction['theta']) :
            ids = (model['categories']['R/theta'] == id)
            #prediction['theta'][id] = np.sum(posterior['theta'][ids, 1])
            #prediction['R'][id] = np.sum(posterior['R'][ids, 1:], 0)#
            prediction['theta'][id] = np.sum(EventFreq[ids, 0]/prediction['EventFreq'][ids])
            prediction['R'][id] = np.sum(EventFreq[ids, 1:].T/prediction['EventFreq'][ids], 1)
            if np.sum(prediction['R'][id]) < .001 * prediction['theta'][id] :
                prediction['theta'][id] = np.sum(prediction['R'][id])/.001
            tot_event = (prediction['theta'][id]+np.sum(prediction['R'][id]))
            prediction['theta'][id] = prediction['theta'][id]/tot_event
            prediction['R'][id] = prediction['R'][id]/tot_event

        for id, delta in enumerate(prediction['delta']) :
            ids = (model['categories']['delta'] == id)
            prediction['delta'][id] = np.sum(posterior['delta'][ids, :, 1])/np.sum(posterior['delta'][ids, :, 0])
            if prediction['delta'][id] > .01 or prediction['delta'][id] < .00001 :
                prediction['delta'][id] = min(max(model['delta'][id], .00001), .01)

        for id, v in enumerate(prediction['v']) :
            ids = (model['categories']['nu'] == id)
            if np.sum(posterior['v'][ids, 0]) > 0 :
                prediction['v'][id] = np.sum(posterior['v'][ids, 1])/np.sum(posterior['v'][ids, 0])
            if np.sum(posterior['v2'][ids, 0]) > 0 :
                prediction['v2'][id] = np.sum(posterior['v2'][ids, 1])/np.sum(posterior['v2'][ids, 0])

            if prediction['v'][id] < 0.001 or prediction['v'][id] > 0.7:
                prediction['v'][id] = min(max( model['v'][id], 0.001 ), 0.7)
            if prediction['v2'][id] < 0.001 or prediction['v2'][id] > 0.7:
                prediction['v2'][id] = min(max( model['v2'][id], 0.001 ), 0.7)
            #max_br = np.max(prediction['EventFreq'][ids])
            #if self.n_b > 2 and prediction['v2'][id] < max_br * 0.5 + 0.001 :
                #prediction['v2'][id] = max( max_br * 0.5 + 0.001, model['v2'][id] )
        return prediction

    def update_branch_parameters(self, model, lower_limit=False) :
        branch_params = []

        for brId, d in enumerate(model['EventFreq']) :
            rId, dId, vId = model['categories']['R/theta'][brId], model['categories']['delta'][brId], model['categories']['nu'][brId]
            if lower_limit and d < 1./self.n_base:
                d = 1./self.n_base
            m, r = (d * model['theta'][rId], d * model['R'][rId])

            a = np.zeros(shape=[self.n_a, self.n_a])
            a[0, 1:] = r
            a[1:, 0] = model['delta'][dId]
            if brId in model['categories'].get('noRec', {}) :
                m += np.sum(r)
                a[0, 1:] = 1e-300
                a[1:, 0] = 1-1e-6

            np.fill_diagonal(a, 1-np.sum(a, 1))
            b = np.zeros(shape=[self.n_a, self.n_b])
            #h = [1./n+(n-2.)/n*model['h'][0] for n in np.arange(1, self.n_b)]
            h = [ ((1-model['h'][0])*(model['h'][0]**n))**(1/(n+1)) for n in np.arange(self.n_b-1) ]

            v, v2 = model['v'][vId], model['v2'][vId]

            b[0]  = [ 1-m ] + [m * hh for hh in h]
            extra = [ 1-v ] + [v * hh for hh in h]
            intra = [ 1-m*h[0]*0.5-v2,  m*h[0]*0.5 ] + [ v2 ] * (self.n_b - 2) # additional SNPs only get in after intra-population recombinations
            mixed = [ 1-v*h[0]*0.5-v2,  v*h[0]*0.5 ] + [ v2 ] * (self.n_b - 2)
            b[1] = intra if self.n_a == 2 and self.n_b > 2 else extra

            if self.n_a > 2 :
                b[2] = intra
                if self.n_a > 3 :
                    b[3] = mixed
            b[b.T[0] < 0.01, 0] = 0.01
            b[2:, 1] = 1 - b[2:, 0] - b[2:, 2]

            branch_params.append(dict(
                pi = np.array([1.] + [0. for i in range(1, self.n_a)]),
                a = a,
                b = b,
            ))
        return branch_params

    def iter_branch_measure(self, data) :
        obs, param, gammaOnly = data
        interval = np.max([np.max(o.T[4]) for o in obs] + [50])
        dist_transition = np.zeros(shape=[interval, param['a'].shape[0], param['a'].shape[1]] )
        dist_transition_adj = np.zeros(shape=[interval] )
        
        a2, a2x, saturate_id = update_distant_transition( param['a'], param['b'].T, dist_transition, dist_transition_adj )
        new_params = []
        for o in obs :
            alpha_Pr, alpha, beta = self.forward_backward(o, pi=param['pi'], a2s=[a2, a2x], b=param['b'])
            new_param = self.estimate_params(param['a'], param['b'], o, alpha, beta, a2, saturate_id, gammaOnly)
            new_param['probability'] = alpha_Pr
            new_params.append(new_param)
        new_param = {'a':[], 'b':[], 'probability':[], 'gamma':[]}
        for k in new_param :
            new_param[k] = np.sum(np.array([ p.get(k) for p in new_params if k in p ]), 0) if k != 'gamma' \
                else [p.get(k) for p in new_params if k in p]
        
        return new_param

    def get_branch_measures(self, params, observations, gammaOnly=False) :
        #branch_measures = list(map(functools.partial(_iter_branch_measure, self), zip(observations, params, [gammaOnly for p in params])))
        branch_measures = pool.map(functools.partial(_iter_branch_measure, self), zip(observations, params, [gammaOnly for p in params]))
        return branch_measures

    def estimate_params(self, transition, emission, obs, alpha, beta, tr2, saturate_id, gammaOnly=False) : # mode = accurate
        n_a, n_b = self.n_a, self.n_b
        gamma = alpha*beta
        gamma = (gamma.T/np.sum(gamma, 1)).T
    
        na = np.dot(alpha[0], tr2[saturate_id])*emission.T[0]
        nb = np.dot(beta[0], tr2[saturate_id].T)
        ng = na*nb/np.sum(na*nb)
        ne = np.array([na for i in xrange(n_a)]).T * np.array([nb*emission.T[0] for i in xrange(n_a)])*transition
        ne /= np.sum(ne)
    
        a2 = np.zeros(shape=[n_a, n_a])
        b2 = np.zeros(shape=[n_a, n_b])
    
        for o, g in zip(obs, gamma) :
            b2[:, o[3]] += g
        
        for o, s, e in zip(obs[1:], alpha[:-1], beta[1:]) :
            d = o[4] - 1
            if d > 2*saturate_id :
                a2 += (d - 2*saturate_id)*ne
                b2[:, 0] += (d - 2*saturate_id)*ng
                d = 2 * saturate_id
    
            if d > saturate_id :
                a, b = np.zeros(shape=[2, d, n_a])
                a[:], b[:] = na, nb
                a[:saturate_id] = np.dot(s, tr2[:saturate_id])*emission.T[0]
                b[-saturate_id:] = np.dot(e*emission.T[o[3]], tr2[:saturate_id].transpose((0, 2, 1)))[::-1]
            else :
                a = np.dot(s, tr2[:d])*emission.T[0]
                b = np.dot(e*emission.T[o[3]], tr2[:d].transpose((0, 2, 1)))[::-1]
    
            g = a*b
            g = g.T/np.sum(g, 1)
    
            b2[:, 0] += np.sum(g, 1)
            if not gammaOnly :
                s1 = np.zeros(shape=[d+1, n_a, 1])
                s1[0, :, 0] = s
                s1[1:, :, 0] = a
            
                s2 = np.zeros(shape=[d+1, 1, n_a])
                s2[:-1, 0, :] = (b*emission.T[0])
                s2[-1, 0, :] = e*emission.T[o[3]]
    
                t = np.matmul(s1, s2) * transition.reshape([1] + list(transition.shape))
                a2 += np.sum(t.T/np.sum(t, axis=(1,2)), 2).T
    
        a2[0] += gamma[0]
        a2.T[0] += gamma[-1]
        if gammaOnly :
            return dict(b=b2, gamma=gamma)
        else :
            return dict(a=a2, b=b2)

    def forward_backward(self, obs, pi, a2s, b) :
        bv = b.T
        a2, a2x = a2s
        alpha, alpha_Pr = np.zeros(shape=[obs.shape[0], self.n_a]), 0.
        alpha[0] = np.dot(pi, a2[0]) * bv[obs[0, 3]]
        alpha_Pr = np.sum(alpha[0])
        alpha[0], alpha_Pr = alpha[0]/alpha_Pr, np.log(alpha_Pr)
        for id, o in enumerate(obs[1:]) :
            r = np.dot(alpha[id], a2[o[4]-1]) * bv[o[3]]
            s = np.sum(r)
            alpha[id+1] = r/s
            alpha_Pr += np.log(s) + a2x[o[4]-1]

        beta, beta_Pr = np.ones(shape=[obs.shape[0], self.n_a]), 0.
        beta[-1] = np.dot(pi, a2[0].T)
        for i, o in enumerate(obs[:0:-1]) :
            id = obs.shape[0]-1-i
            r = np.dot(beta[id] * bv[o[3]], a2[o[4]-1].T)
            s = np.sum(r)
            beta[id-1] = r/s
            beta_Pr += np.log(s) + a2x[o[4]-1]
        return alpha_Pr, alpha, beta

    def get_brLens(self, branches, n_base) :
        return np.array([ np.sum(branch.T[1] > 0)/float(n_base) for branch in branches ])

    def prepare_branches(self, mutations, sequences, missing, interval=None) :
        if not interval :
            interval = np.sum([s[1] for s in sequences])
        
        
        branches = np.unique(mutations.T[0])
        mutations = np.hstack([mutations, np.zeros([mutations.shape[0], 2], dtype=int)])
        blocks = []
        
        for seqId, (seqName, seqLen) in enumerate(sequences) :
            region = [[seqId, 1, seqLen]]
            if missing.size :
                for ms in missing[(missing.T[0] == seqId) & (missing.T[2] - missing.T[1] + 1 >= np.min([5000, seqLen]))] :
                    region[-1][2] = ms[1] - 1
                    region.append([seqId, ms[2]+1, seqLen])
            region = [ (r[0], r[1], r[2], r[2]-r[1]+1, 0) for r in region if r[2]>=r[1] ]
            blocks.extend(region)
        blocks = np.array(blocks, dtype=int)
        blocks.T[4] = np.arange(blocks.shape[0])
        anchors = [np.vstack([ np.vstack([branches, np.repeat(block[0], branches.size), np.repeat(block[1]-1, branches.size), np.zeros([3, branches.size])]).T, 
                               np.vstack([branches, np.repeat(block[0], branches.size), np.repeat(block[2]+1, branches.size), np.zeros([3, branches.size])]).T ] ) for block in blocks]
        mutations = np.vstack([mutations] + anchors).astype(int)
        mutations = mutations[np.lexsort(mutations.T[::-1])]
        self.n_base = 0
        for seqId, (seqName, seqLen) in enumerate(sequences) :
            regions = blocks[blocks.T[0] == seqId]
            if regions.shape[0] == 0 : continue
            s = np.zeros([2, seqLen+2], dtype=int)
            s[1] = np.arange(seqLen+2)
            if missing.size :
                for ms in missing[(missing.T[0] == seqId)] :
                    s[1, ms[1]:ms[2]+1] = -2
                    s[1, ms[2]+1:] -= ms[2]-ms[1]+1
            for r in regions :
                if r[4] == 55 :
                    print ('')
                s[0, r[1]-1:r[2]+2] = r[4]
                if s[1, r[1]] > 1 :
                    s[1, r[1]-1] = s[1, r[1]] - 1
                else :
                    ss = s[1, r[1]:r[2]+1]
                    s[1, r[2]+1] = np.min(ss[ss>0]) - 1
                s[1, r[1]-1:r[2]+1] -= s[1, r[1]-1]
                if s[1, r[2]] > 0 :
                    s[1, r[2]+1] = s[1, r[2]] + 1
                else :
                    s[1, r[2]+1] = np.max(s[1, r[1]:r[2]+1]) + 1
                
                r[3] = np.sum(s[1, r[1]:r[2]+1]>0)
            self.n_base += int(np.sum(s[1] >= 0))
            mutations[mutations.T[1] == seqId, -2:] = s.T[mutations[mutations.T[1] == seqId, 2]]
            mut = np.copy(mutations[(mutations.T[1] == seqId)]).astype(float)
            mut.T[2] = 0
            mut[:-1, 2] = mut[1:, 5] - mut[:-1, 5]
            anchors = mut[:-1][(mut[:-1, 2] > interval) & (mut[:-1, 3]+mut[1:, 3] > 0)]
            anchors.T[3] = np.ceil(anchors.T[2]/interval)-1
            anchors[anchors.T[3] > 3, 3] = 3
            anchors.T[1] = anchors.T[2] / (anchors.T[3]+1)
            while anchors.size :
                anchors.T[5] += anchors.T[1]
                mutAnchors = np.vstack([anchors.T[0], np.repeat(seqId, anchors.shape[0]), np.repeat(-1, anchors.shape[0]), np.zeros(anchors.shape[0]), anchors.T[4], np.round(anchors.T[5])]).T.astype(int)
                mutations = np.vstack([mutations, mutAnchors])
                anchors.T[3] -= 1
                anchors = anchors[anchors.T[3] > 0]
        mutations = mutations[np.lexsort(mutations.T[[5, 4, 0]])]
        mutations[mutations.T[2] == 0, 2] = 1
        mutations[mutations.T[2] > blocks[mutations.T[4], 2], 2] -= 1
        if self.n_b > 2 :
            self.n_b = int(max(3, np.max(mutations.T[3]) + 1))
        else :
            mutations[mutations.T[3] > 1, 3] = 1
        
        def prepare_obs(obs, blocks, interval=None) :
            res = [ obs[obs.T[4] == blkId] for blkId, block in enumerate(blocks) ]
            for r in res :
                r[1:, 4] = np.diff(r.T[5])
            return res
        return [prepare_obs(mutations[mutations.T[0] == brId], blocks, interval) for brId in np.unique(mutations.T[0])]

    def predict(self, mutations, branches, sequences, missing, marginal, prefix='RecHMM', tree=None) :
        assert self.model, 'No model'
        self.observations = self.prepare_branches(mutations, sequences, missing)
        self.branches = branches if branches is not None else np.arange(len(self.observations)).astype(str)
        self.sequences = sequences if sequences is not None else [[str(id), 0] for id, _ in enumerate(self.observations[0])]

        stats = self.margin_predict(marginal) if marginal > 0. and marginal <= 1. else self.map_predict()
        
        with open(prefix+'.recombination.region', 'w') as rec_out: #, uopen(prefix+'.mutations.status.gz', 'w') as mut_out :
            rec_out.write('#Branch\tname\tmutationRate\trecombinationRate\tMutationCoverage\n')
            rec_out.write('#\tImportation\tseqName\tstart\tend\ttype\tscore\n')
            for name in self.branches :
                stat = stats[name]
                rec_out.write('Branch\t{0}\tM={1:.5e}\tR={2:.5e}\tB={3:.3f}\n'.format(name, stat['M'], stat['R'], stat['weight_p'][0]))
                for r in stat['sketches'] :
                    rec_out.write('\tImportation\t{0}\t{1}\t{2}\t{3}\t{4}\t{5:.3f}\n'.format(name, self.sequences[r[0]][0], r[1], r[2], ['External', 'Internal', 'Mixed   '][r[3]-1], r[4]))
        if tree :
            from ete3 import Tree
            tre = Tree(tree, format=1)
            for node in tre.traverse() :
                if node.name in stats :
                    node.dist = stats[node.name]['M']
            tre.write(format=1, outfile=prefix + '.mutational.tre')
            print('Mutational tree is written in {0}'.format(prefix+'.mutational.tre'))

        print('Imported regions are reported in {0}'.format(prefix+'.recombination.region'))

    def map_predict(self) :
        self.screen_out('Predict recombination sketches using', self.model)
        branch_params = self.update_branch_parameters(self.model, lower_limit=True)
        status = pool.map(functools.partial(_iter_viterbi, self), zip(self.observations, branch_params))
        #status = list(map(functools.partial(_iter_viterbi, self), zip(self.observations, branch_params)))
        res = {}
        for name, dm, dr, stat in zip(self.branches, self.model['posterior']['theta'], self.model['posterior']['R'], status) :
            rec_len = np.sum([ e1-s1+1 for c, s, e, t, s1, e1, p in stat['sketches'] ])
            res[name] = dict(sketches=[ k[:4]+k[6:] for k in stat['sketches']],
                             weight_p=np.array([self.n_base-rec_len, rec_len], dtype=float),
                             M=dm[1]/dm[0],
                             R=np.sum(dr[1:])/dr[0])
        return res

    def margin_predict(self, marginal=0.9) :
        branch_params = self.update_branch_parameters(self.model)
        status = self.get_branch_measures(branch_params, self.observations, gammaOnly=True)
        res = {}
        for name, dm, dr, observation, stat in zip(self.branches, self.model['posterior']['theta'], self.model['posterior']['R'], self.observations, status) :
            path = []
            for obs, gamma in zip(observation, stat['gamma']) :
                for id, (o, s) in enumerate(zip(obs, gamma)) :
                    p = np.argmax(s)
                    if p > 0 and s[0] < 0.5 :
                        if len(path) == 0 or path[-1][3] != p or path[-1][5] != obs[id-1][5] :
                            if o[2] >= 0 :
                                path.append([o[1], o[2], o[2], p, o[5], o[5], 1-s[0]])
                        else :
                            path[-1][5] = o[5]
                            if 1-s[0] > path[-1][6] :
                                path[-1][6] = 1-s[0]
                            if o[2] >= 0 :
                                path[-1][2] = o[2]
            res[name] = dict(sketches=[p[:4]+p[6:] for p in path if p[2]-p[1] > 0 and p[6] >= marginal],
                             weight_p=np.sum(stat['b'], 1),
                             M=dm[1]/dm[0],
                             R=np.sum(dr[1:])/dr[0])
        return res


    def viterbi(self, data) :
        observation,  params = data
        pi, a, b = params['pi'], params['a'], params['b']
        bv = b.T
        regions = []
        for obs in observation :
            rsite = dict(obs[:, np.array([5,2])])
            seqName, n_base = obs[0, 1], obs[-1, -1] + 1
            path = np.zeros(shape=[n_base, self.n_a], dtype=int)
            alpha = np.zeros(shape=[n_base, self.n_a])
            alpha[0] = np.log( np.dot(pi, a) * bv[obs[0, 3]])
    
            i = 0
            p = np.zeros([self.n_a, self.n_a])
            a[a==0], b[b==0] = 1e-300, 1e-300
            
            pa, pb = np.log(a), np.log(b)
            ids = np.arange(self.n_a)
            for br, sn, s0, o, d, s in obs[1:] :
                for dd in np.arange(d-1) :
                    i += 1
                    p.T[:] = alpha[i-1]
                    p = p + pa + pb.T[0]
                    path[i] = np.argmax(p, 0)
                    alpha[i] = p[path[i], ids]
                    if np.max(path[i]) == 0 :
                        j = i + d - 1 - dd
                        alpha[i:j] = alpha[i]
                        alpha[i:j] += ((pa[0, 0] + pb[0, 0]) * np.arange(d-1-dd)).reshape(d-1-dd, 1)
                        i = j - 1
                        break
                i += 1
                p.T[:] = alpha[i-1]
                p = p + pa + pb.T[o]
                path[i] = np.argmax(p, 0)
                alpha[i] = p[path[i], ids]
            alpha[i] += np.log( np.dot(pi, a.T) )
            max_path = np.argmax(alpha[i])
            inrec = np.zeros(path.shape[0])
            for id in np.arange(path.shape[0]-2, 0, -1) :
                max_path = path[id+1, max_path]
                if max_path > 0 :
                    inrec[id] = 1
                    if len(regions) == 0 or regions[-1][4] != id + 1 :
                        regions.append([seqName, -1, -1, max_path, id, id, 1.])
                    else :
                        regions[-1][4] = id
                    if id in rsite :
                        if regions[-1][2] == -1 :
                            regions[-1][2] = rsite[id]
                        regions[-1][1] = rsite[id]
        regions = [r for r in regions if r[2] >= 0]
        return dict(sketches=sorted(regions), gamma=1.-inrec[ obs.T[5] ])

    def report(self, bootstrap, prefix='') :
        if 'posterior' in self.model :
            posterior = self.model['posterior']
        else :
            posterior = []
        n_br = posterior['theta'].shape[0]
        bs = np.random.randint(n_br, size=[bootstrap, n_br])

        h0 = 1-self.model['h'][0]/2

        reports = {}

        reports['probability'] = [self.model['probability']]
        reports['BIC'] = [-2*reports['probability'][0] + self.n_a*self.n_b*np.log(self.n_base*n_br) ]

        EventFreq = np.vstack([posterior['theta'].T[1]/posterior['theta'].T[0], posterior['R'][:, 1:].T/posterior['R'][:, 0]]).T
        reports['EventFreq'] = [np.sum(EventFreq), np.sum(EventFreq[bs], (1,2))]

        #reports['theta'] = [np.sum(posterior['theta'].T[1]), np.sum(posterior['theta'].T[1][bs], 1)]
        #reports['R'] = [np.sum(posterior['R'].T[1:]), np.sum(np.sum(posterior['R'].T[1:], 0)[bs])]
        reports['theta'] = [np.sum(EventFreq[:, 0]/self.model['EventFreq']), np.sum(EventFreq.T[0][bs]/self.model['EventFreq'][bs], 1)]
        reports['R'] = [np.sum(np.sum(EventFreq[:, 1:], 1)/self.model['EventFreq']), np.sum(np.sum(EventFreq.T[1:], 0)[bs]/self.model['EventFreq'][bs], 1)]
        tot = [reports['theta'][0]+reports['R'][0], reports['theta'][1]+reports['R'][1]]
        reports['theta'] = [reports['theta'][0]/tot[0], reports['theta'][1]/tot[1]]
        reports['R'] = [reports['R'][0]/tot[0], reports['R'][1]/tot[1]]

        reports['delta'] = [np.sum(posterior['delta'][:, :, 0])/np.sum(posterior['delta'][:, :, 1]),
                            np.sum(posterior['delta'][:, :, 0][bs], (1,2))/np.sum(posterior['delta'][:, :, 1][bs], (1, 2))]

        reports['nu'] = [np.sum(posterior['v'][:, 1])/np.sum(posterior['v'][:, 0]),
                            np.sum(posterior['v'][:, 1][bs], 1)/np.sum(posterior['v'][:, 0][bs], 1)]

        reports['nu(in)'] = [np.sum(posterior['v2'][:, 1])/np.sum(posterior['v2'][:, 0]),
                            np.sum(posterior['v2'][:, 1][bs], 1)/np.sum(posterior['v2'][:, 0][bs], 1)]

        reports['R/theta'] = [reports['R'][0]/reports['theta'][0], reports['R'][1]/reports['theta'][1]]

        reports['r/m'] = [np.sum(posterior['v'].T[1]+posterior['v2'].T[1])/np.sum(posterior['theta'].T[1]), np.sum(posterior['v'].T[1][bs]+posterior['v2'].T[1][bs], 1)/np.sum(posterior['theta'].T[1][bs], 1)]
        with open(prefix + '.best.model.report', 'w') as fout :
            fout.write( 'Prefix    \tParameter \tValue     \tSTD       \tCI 95% (Low - High)\n' )
            sys.stdout.write( 'Prefix    \tParameter \tValue     \tSTD       \tCI 95% (Low - High)\n' )
            for key in ('R/theta', 'r/m', 'delta', 'nu', 'nu(in)', 'EventFreq', 'theta', 'R') :
                if key == 'delta' :
                    fout.write( '{0}\t{1}\t{2:.4f}\t{3:.4f}\t{4:.4f} - {5:.4f}\n'.format(prefix.ljust(10), key.ljust(10), reports[key][0], np.std(reports[key][1]), *np.sort(reports[key][1])[[int(bootstrap*0.025), int(bootstrap*0.975)]].tolist()) )
                    sys.stdout.write( '{0}\t{1}\t{2:.4f}\t{3:.4f}\t{4:.4f} - {5:.4f}\n'.format(prefix.ljust(10), key.ljust(10), reports[key][0], np.std(reports[key][1]), *np.sort(reports[key][1])[[int(bootstrap*0.025), int(bootstrap*0.975)]].tolist()) )
                else :
                    fout.write( '{0}\t{1}\t{2:.6f}\t{3:.6f}\t{4:.6f} - {5:.6f}\n'.format(prefix.ljust(10), key.ljust(10), reports[key][0], np.std(reports[key][1]), *np.sort(reports[key][1])[[int(bootstrap*0.025), int(bootstrap*0.975)]].tolist()) )
                    sys.stdout.write( '{0}\t{1}\t{2:.6f}\t{3:.6f}\t{4:.6f} - {5:.6f}\n'.format(prefix.ljust(10), key.ljust(10), reports[key][0], np.std(reports[key][1]), *np.sort(reports[key][1])[[int(bootstrap*0.025), int(bootstrap*0.975)]].tolist()) )
            fout.write('{0}\tBIC       \t{1}\n'.format(prefix.ljust(10), reports['BIC'][0]))
            sys.stdout.write('{0}\tBIC       \t{1}\n'.format(prefix.ljust(10), reports['BIC'][0]))
        print('Global parameters are summarized in {0}'.format(prefix + '.best.model.report'))
        return


def parse_arg(a) :
    parser = argparse.ArgumentParser(description='Parameters for RecHMM. ', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--data', '-d', help='A list of mutations generated by EnPhyl', required=True)
    parser.add_argument('--model', '-m', help='Read a saved best model.', default='')
    parser.add_argument('--task', '-t', help='task to run. \n0: One rec category from external sources.\n1: Three rec categories considering internal, external and mixed sources [default].', default=1, type=int)
    parser.add_argument('--init', '-i', help='Initiate models with guesses of recombinant proportions. \nDefault: 0.05,0.5,0.95', default='0.05,0.5,0.95')
    parser.add_argument('--prefix', '-p', help='Prefix for all the outputs ', default='RecHMM')
    parser.add_argument('--cool_down', '-c', help='Delete the worst model every N iteration. Default:5', type=int, default=5)
    parser.add_argument('--n_proc', '-n', help='Number of processes. Default: 5. ', type=int, default=5)
    parser.add_argument('--bootstrap', '-b', help='Number of Randomizations for confidence intervals. \nDefault: 1000. ', type=int, default=1000)
    parser.add_argument('--report', '-r', help='Only report the model and do not calculate external sketches. ', default=False, action="store_true")
    parser.add_argument('--marginal', '-M', help='Find recombinant regions using marginal likelihood rather than [DEFAULT] maximum likelihood method. \n[DEFAULT] 0 to use Viterbi algorithm to find most likely path.\n Otherwise (0, 1) use forward-backward algorithm, and report regions with >= M posterior likelihoods as recombinant sketches.', default=0, type=float)
    parser.add_argument('--tree', '-T', help='[INPUT, OPTIONAL] A labelled tree. Only used to generate corresponding mutational tree.', default=None)
    parser.add_argument('--clean', '-v', help='Do not show intermediate results during the iterations.', default=False, action='store_true')
    parser.add_argument('--local_r', '-lr', help='Specify a comma-delimited list of branches that share a different R/theta (Frequency of rec) ratio than the global consensus. Can be specified multiple times. \nUse "*" to assign different value for each branch.', default=[], action="append")
    parser.add_argument('--local_nu', '-ln', help='Specify a comma-delimited list of branches that share a different Nu (SNP density in rec) than the global consensus. Can be specified multiple times. \nUse "*" to assign different value for each branch.', default=[], action="append")
    parser.add_argument('--local_delta', '-ld', help='Specify a comma-delimited list of branches that share a different Delta (Length of rec sketches) than the global consensus. Can be specified multiple times. \nUse "*" to assign different value for each branch.', default=[], action="append")

    args = parser.parse_args(a)
    args.categories = { 'R/theta':{},
                        'nu':{},
                        'delta':{} }
    for variable, categories in zip(['R/theta', 'nu', 'delta'], [args.local_r, args.local_nu, args.local_delta]) :
        for id, category in enumerate(categories) :
            branches = category.split(',')
            args.categories[variable].update({ br:id+1 for br in branches })
        if '*' in args.categories[variable] :
            args.categories[variable] = {'*': 0}
    return args

def RecHMM(args) :
    args = parse_arg(args)
    global pool, verbose
    pool = Pool(args.n_proc)
    verbose = not args.clean

    model = recHMM(mode=args.task)
    
    if not args.report or not args.model :
        sequences, missing = [], []
        with uopen(args.data) as fin :
            for line in fin :
                if line.startswith('##') :
                    if line.startswith('## Sequence_length:') :
                        part = line[2:].strip().split()
                        sequences.append([part[1], int(part[2])])
                    elif line.startswith('## Missing_region:') :
                        part = line[2:].strip().split()
                        missing.append([part[1], int(part[2]), int(part[3])])
                else :
                    break
            data = pd.read_csv(fin, sep='\t', dtype=str, header=None).values
        branches, mutations = {}, []
        seqLens = {seqName:[seqId, seqLen] for seqId, (seqName, seqLen) in enumerate(sequences)}
        for d in data :
            if re.findall(r'^[ACGTacgt]->[ACGTacgt]$', d[4]) :
                if d[1] not in seqLens :
                    seqLens[d[1]] = [len(seqLens), int(d[2])]
                if seqLens[d[1]][1] < int(d[2]) :
                    seqLens[d[1]][1] = int(d[2])
                if d[0] not in branches :
                    branches[d[0]] = len(branches)
                brId, seqId = branches[d[0]], seqLens[d[1]][0]
                mutations.append([brId, seqId, int(d[2]), int(d[3])])
        missing = np.array([ [seqLens.get(m[0], [-1])[0], m[1], m[2]] for m in missing ])
        sequences = [ [n, i[1]] for n, i in sorted(seqLens.items(), key=lambda x:x[1][0])]
        branches = np.array([ br for br, id in sorted(branches.items(), key=lambda x:x[1]) ])
        mutations = np.array(mutations)
        reorder = np.argsort(-np.bincount(mutations.T[0]))
        branches = branches[reorder]
        reorder = np.array([i1 for i1, i2 in sorted(enumerate(reorder), key=lambda x:x[1])])
        mutations.T[0] = reorder[mutations.T[0]]
    if args.model :
        model.load(open(args.model, 'r'))
    else :
        model.fit(mutations, branches=branches, sequences=sequences, missing=missing, categories=args.categories, init=args.init, cool_down=args.cool_down)
        model.save(open(args.prefix + '.best.model.json', 'w'))
        print('Best HMM model is saved in {0}'.format(args.prefix + '.best.model.json'))
    model.report(args.bootstrap, args.prefix)

    if not args.report :
        model.predict(mutations, branches=branches, sequences=sequences, missing=missing, marginal=args.marginal, prefix=args.prefix, tree=args.tree)

pool, verbose = None, True
if __name__ == '__main__' :
    RecHMM(sys.argv[1:])

