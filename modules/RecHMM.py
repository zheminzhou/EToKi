import numpy as np, pandas as pd, sys, os, copy, argparse
import functools, time, datetime
from multiprocessing import Pool

def _iter_branch_measure(obj, arg) :
    return obj.iter_branch_measure(arg)

def _iter_viterbi(obj, arg) :
    return obj.viterbi(arg)

class recHMM(object) :
    def __init__(self, mode=1) :
        self.max_iteration = 200
        self.n_base = None
        self.mode = ['legacy', 'hybrid', 'intra', 'both'][mode]
        self.n_a, self.n_b = [[2, 2], [4, 3], [2, 3], [3, 3]][mode]

    def fit(self, branches, names=None, categories=None, init=None, cool_down=5, global_rec=False) :
        self.prepare_branches(branches)
        self.names = np.array(names) if names else np.arange(len(self.observations)).astype(str)
        self.categories = { 'noRec':{} }
        for c, assigns in categories.iteritems() :
            self.categories[c] = np.zeros(shape=[len(self.names)], dtype=int)
            if assigns.get('*', -1) == 0 :
                self.categories[c][:] = np.arange(len(self.names))
            else :
                for i, n in enumerate(self.names) :
                    self.categories[c][i] = categories[c].get(n, 0)

        models = self.initiate(self.observations, init=init)
        return self.BaumWelch(models, self.max_iteration, cool_down=cool_down, global_rec=global_rec)

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
            model['posterior'] = { k:v.tolist() for k,v in model['posterior'].iteritems() }
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
            model['posterior'] = { k:np.array(v) for k,v in model['posterior'].iteritems() }
        if 'categories' in model:
            for k in ('R/theta', 'nu', 'delta') :
                model['categories'][k] = np.array(model['categories'][k])
        self.__dict__.update(model.pop('model', {}))
        self.model = model
        return self.model


    def initiate(self, observations, init) :
        criteria = np.array(init.split(',')).astype(float)
        intervals = np.sort(np.concatenate([ np.diff(observation[observation.T[2] > 0, 0]) for observation in observations ]))
        criteria = (criteria * intervals.size).astype(int)
        cutoffs = np.unique(intervals[criteria])

        self.models = []
        for cutoff in cutoffs :
            mut_summary, rec = [], []
            for observation in observations :
                substitution = observation[observation.T[2] > 0]
                if len(substitution) :
                    d, p, po = substitution[0, 0] - 1 + self.n_base - substitution[-1, 0], substitution[0, 0], int(substitution[0, 2] > 1)
                    blocks = [[ False, d, 1, po + int(substitution[-1, 2] > 1) ]]
                    for s, d, o, _, _ in substitution[1:] :
                        o, d = int(o>1), s-p
                        inRec = ( d <= cutoff )
                        if blocks[-1][0] != inRec :
                            blocks.append([ inRec, d, 1, o + po ])
                        else :
                            blocks[-1][1:4] = blocks[-1][1] + d, blocks[-1][2] + 1, blocks[-1][3] + o + po
                        p, po = s, o
                    blocks = np.array(blocks, dtype=float)
                    mut_region, rec_region = blocks[blocks.T[0] == 0], blocks[blocks.T[0] == 1]
                    mut_summary.append([ np.sum(mut_region.T[1]), np.sum(mut_region.T[2]), np.sum(mut_region.T[3])/2, np.min([mut_region.shape[0], rec_region.shape[0]]).astype(float) ])
                    rec.extend(rec_region)
            mut_summary = np.array(mut_summary)

            EventFreq = np.sum(mut_summary.T[[1, 3]]/mut_summary.T[0], 0)

            n_br, (bases, muts, homos, recs) = mut_summary.shape[0], np.sum(mut_summary, 0)
            model = dict(theta = np.array([ muts/np.sum(EventFreq)/self.n_base for id in np.unique(self.categories['R/theta']) ]),
                         h     = [ max(0.01, min(0.95, 1./2.*homos/muts)), max(0.01, min(0.95, 5./2.*homos/muts)) ],
                         probability = -1e300,
                         diff        = 1e300,
                         EventFreq   = EventFreq,
                         id          = len(self.models) + 1,
                         ite         = 0,
                         categories  = copy.deepcopy(self.categories)
                        )

            rec = np.array(rec)

            if rec.shape[0] :
                rec2, rec = np.sum(rec[rec.T[3]/rec.T[2] > 2.*2.5*model['h'][0]], 0), np.sum(rec[rec.T[3]/rec.T[2] <= 2.*2.5*model['h'][0]], 0)
                rec[0], rec2[0] = max(rec[0], 1), max(rec2[0], 1)

                model['h'][1]  = max(model['h'][1], 1./2.*(rec2[3]+1)/(rec2[2]+1))
                model['delta'] = np.array([ max(0.00002, min(0.005, (rec[0]+rec2[0])/(rec2[1] + rec[1]))) for id in np.unique(self.categories['delta']) ])
                model['v']     = np.array([ np.max([0.001, np.min([0.7, (rec[2]+rec2[2])/(rec2[1] + rec[1])])]) for id in np.unique(self.categories['nu']) ])
                model['v2']    = np.array([ np.min([0.7, (rec2[2]+1)/(rec2[1]+1)]) for id in np.unique(self.categories['nu']) ])

                p = np.array([rec[0], np.sqrt(rec[0]*rec2[0]), rec2[0]])[:(self.n_a-1)]
                p = p/np.sum(p)
                model['R'] = np.array([ recs/np.sum(EventFreq)/self.n_base * p for id in np.unique(self.categories['R/theta']) ])
                tot_event = model['theta'][0] + np.sum(model['R'][0])
                model['theta'] = model['theta']/tot_event
                model['R'] = model['R']/tot_event
                self.screen_out('Initiate', model)
                self.models.append(model)
        return self.models


    def BaumWelch(self, models, max_iteration, cool_down=5, global_rec=False) :
        n_model = len(models)
        for ite in np.arange(max_iteration) :
            new_models = []
            self.model = models[0]

            for model in models:
                if 'diff' in model and model['diff'] < 0.001 :
                    new_models.append(model)
                else :
                    print ''
                    self.screen_out('Assess', model)
                    branch_params = self.update_branch_parameters(model, global_rec=global_rec)
                    branch_measures = self.get_branch_measures(branch_params, self.observations)

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
                            prediction['id'] += 0.01
                            prediction['ite'] = ite+1
                            prediction['diff'] = 1
                            new_models.append(prediction)

            new_models = sorted(new_models, key=lambda x:-x['probability'])
            if ite > 0 and ite % cool_down == 0 :
                while len(new_models) > 1 and n_model - len(new_models) < ite/cool_down :
                    if new_models[-1]['probability'] + new_models[-1]['diff']/2 >= new_models[-2]['probability'] + new_models[-2]['diff']/2 :
                        new_models[-1]['diff'] = new_models[-1]['diff']/2
                        new_models[-2]['diff'] = new_models[-2]['diff']/2
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
                if v[0] > .0001 * theta[0] and v[1]*theta[0] < theta[1]*v[0] :
                    print 'Model {0}: Branch {1} is too divergent for Rec identification. Move to new category. '.format(model['id'], self.names[brId])
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
                    tId = model['categories']['R/theta'][brId]
                    if tId not in model['categories']['low_cov'] :
                        if np.sum(model['categories']['R/theta']== tId) > 1:
                            model['categories']['low_cov'][tId] = model['theta'].size
                            model['categories']['low_cov'][model['theta'].size] = model['theta'].size
                            model['R'] = np.vstack([model['R'], [model['R'][tId]]])
                            model['theta'] = np.concatenate([model['theta'], [model['theta'][tId]]])
                        else :
                            model['categories']['low_cov'][tId] = tId
                    if model['categories']['low_cov'][tId] != tId :
                        print 'Model {0}: Branch {1} got too much recombined regions. Move to new category. '.format(model['id'], self.names[brId])
                        model['categories']['R/theta'][brId] = model['categories']['low_cov'][tId]
                        model['probability'] = -1e300
                        model['diff'] = 1e300

    def screen_out(self, action, model) :
        if not verbose :
            return
        print '{2}\t{0} model {id}[{ite}] - BIC: {3:.8e} - EventFreq: {1:.3e}; theta: {theta[0]:.3f}; R: {5:.3f}; delta: {4:.3e};  Nu: {v[0]:.3e},{v2[0]:.3e}; h: {h[0]:.3f},{h[1]:.3f}'.format(action, np.sum(model['EventFreq']), str(datetime.datetime.now())[:19],  -2*model['probability'] + self.n_a*self.n_b*np.log(self.n_base*len(self.observations)), 1/model['delta'][0], np.sum(model['R'][0]), **model)
        sys.stdout.flush()

    def estimation(self, model, branch_measures) :
        h = [1-model['h'][0]/2, np.sqrt(model['h'][0]*(2-model['h'][0]))/2.]

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
                posterior['theta'][id, :] = [ np.sum(b[0]) + np.sum(b[1])*h[0], np.sum(b[0, 1:]) + b[1, 1] ]
                posterior['h'][id] = [ b[0, 1]+b[0, 2], b[0, 2], b[1, 1]+b[1, 2], b[1, 2] ]
                posterior['v'][id] = [np.sum(b[1]), np.sum(b[1, 1:])]
                posterior['v2'][id] = [np.sum(b[1]), b[1, 2]]
            else :
                posterior['theta'][id, :] = [ np.sum(b[0]) + np.sum(b[2:3])*h[0], np.sum(b[0, 1:]) + np.sum(b[2:3, 1]) ]
                posterior['h'][id] = [ np.sum(b[:2, 1:3]), np.sum(b[:2, 2]), np.sum(b[2:, 1:3]), np.sum(b[2:, 2]) ]
                posterior['v'][id] = [np.sum(b[1]) + np.sum(b[3:])*h[0], np.sum(b[1, 1:]) + np.sum(b[3:, 1])]
                posterior['v2'][id] = [np.sum(b[2:]), np.sum(b[2:, 2])]

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
            prediction['v'][id] = np.sum(posterior['v'][ids, 1])/np.sum(posterior['v'][ids, 0])
            prediction['v2'][id] = np.sum(posterior['v2'][ids, 1])/np.sum(posterior['v2'][ids, 0])

            if prediction['v'][id] < 0.001 or prediction['v'][id] > 0.7:
                prediction['v'][id] = min(max( model['v'][id], 0.001 ), 0.7)
                max_br = np.max(prediction['EventFreq'][ids])
                if self.n_b > 2 and prediction['v2'][id] < max_br*h[1] :
                    prediction['v2'][id] = max( max_br*h[1], model['v2'][id] )
        return prediction

    def update_branch_parameters(self, model, lower_limit=False, global_rec=False) :
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
            b = np.zeros(shape=[self.n_a, 3])
            h = [1-model['h'][0]/2, np.sqrt(model['h'][0]*(2-model['h'][0]))/2.]

            v, v2 = model['v'][vId], model['v2'][vId]

            b[0]  = [ 1-m,          m*h[0],  m*h[1] ]
            extra = [ 1-v,          v*h[0],  v*h[1] ]
            intra = [ 1-m*h[0]-v2,  m*h[0],  v2     ]
            mixed = [ 1-v*h[0]-v2,  v*h[0],  v2     ]
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
        a2, a2x, saturate_id = self.update_distant_transition(param['a'], param['b'].T, np.max(obs.T[1]))
        alpha_Pr, alpha, beta = self.forward_backward(obs, pi=param['pi'], a2s=[a2, a2x], b=param['b'])
        new_param = self.estimate_params(param['a'], param['b'], obs, alpha, beta, a2, saturate_id, gammaOnly)
        new_param['probability'] = alpha_Pr
        return new_param

    def get_branch_measures(self, params, observations, gammaOnly=False) :
        branch_measures = pool.map(functools.partial(_iter_branch_measure, self), zip(observations, params, [gammaOnly for p in params]))
        return branch_measures

    def estimate_params(self, transition, emission, obs, alpha, beta, tr2, saturate_id, gammaOnly=False) : # mode = accurate
        gamma = alpha*beta
        gamma = (gamma.T/np.sum(gamma, 1)).T

        na = np.dot(alpha[0], tr2[saturate_id])*emission.T[0]
        nb = np.dot(beta[0], tr2[saturate_id].T)
        ng = na*nb/np.sum(na*nb)
        ne = np.array([na for i in np.arange(self.n_a)]).T * np.array([nb*emission.T[0] for i in np.arange(self.n_a)])*transition
        ne = ne/np.sum(ne)

        a2 = np.zeros(shape=[self.n_a, self.n_a])
        b2 = np.zeros(shape=[self.n_a, 3])

        for o, g in zip(obs, gamma) :
            b2[:, o[2]] += g

        for o, s, e in zip(obs[1:], alpha[:-1], beta[1:]) :
            d = o[1] - 1
            if d > 2*saturate_id :
                a2 += (d - 2*saturate_id)*ne
                b2[:, 0] += (d - 2*saturate_id)*ng
                d = 2 * saturate_id

            if d > saturate_id :
                a = np.zeros(shape=[d, self.n_a])
                b = np.zeros(shape=[d, self.n_a])
                a[:], b[:] = na, nb
                a[:saturate_id] = np.dot(s, tr2[:saturate_id])*emission.T[0]
                b[-saturate_id:] = np.dot(e*emission.T[o[2]], tr2[:saturate_id].transpose((0, 2, 1)))[::-1]
            else :
                a = np.dot(s, tr2[:d])*emission.T[0]
                b = np.dot(e*emission.T[o[2]], tr2[:d].transpose((0, 2, 1)))[::-1]

            g = a*b
            g = g.T/np.sum(g, 1)

            b2[:, 0] += np.sum(g, 1)
            if not gammaOnly :
                s1 = np.zeros(shape=[d+1, self.n_a, self.n_a])
                s1[0].T[:] = s
                s1[1:].T[:] = a.T

                s2 = np.zeros(shape=[d+1, self.n_a, self.n_a])
                s2[:-1].T[:] = (b*emission.T[0]).T
                s2 = s2.transpose((0, 2, 1))
                s2[-1, :] = e*emission.T[o[2]]

                t = (s1 * s2) * transition.reshape([1] + list(transition.shape))
                a2 += np.sum(t.T/np.sum(t, axis=(1,2)), 2).T
        a2[0] += gamma[0]
        a2.T[0] += gamma[-1]
        if gammaOnly :
            return dict(b=b2, gamma=gamma)
        else :
            return dict(a=a2, b=b2)


    def update_distant_transition(self, transition, emission, interval) :
        interval = np.max([interval, 50])
        dist_transition = np.zeros(shape=[interval, transition.shape[0], transition.shape[1]] )
        dist_transition[0] = transition
        dist_transition_adj = np.zeros(shape=[interval] )
        ss = 0.
        for id in np.arange(interval-1) :
            t = np.dot(transition*emission[0], dist_transition[id])
            s = np.sum(t)/transition.shape[0]
            dist_transition[id+1] = t/s
            ss = ss + np.log(s)
            dist_transition_adj[id+1] = ss
            if np.sum(np.abs(dist_transition[id+1] - dist_transition[id])) <= 1e-12 :
                dist_transition[id+2:] = dist_transition[id+1]
                dist_transition_adj[id+1:] = dist_transition_adj[id+1] + np.arange(interval-id-1)*np.log(s)
                break

        return dist_transition, dist_transition_adj, id

    def forward_backward(self, obs, pi, a2s, b) :
        bv = b.T
        a2, a2x = a2s
        alpha, alpha_Pr = np.zeros(shape=[obs.shape[0], self.n_a]), 0.
        alpha[0] = np.dot(pi, a2[0]) * bv[obs[0, 2]]
        alpha_Pr = np.sum(alpha[0])
        alpha[0], alpha_Pr = alpha[0]/alpha_Pr, np.log(alpha_Pr)
        for id, o in enumerate(obs[1:]) :
            r = np.dot(alpha[id], a2[o[1]-1]) * bv[o[2]]
            s = np.sum(r)
            alpha[id+1] = r/s
            alpha_Pr += np.log(s) + a2x[o[1]-1]

        beta, beta_Pr = np.ones(shape=[obs.shape[0], self.n_a]), 0.
        beta[-1] = np.dot(pi, a2[0].T)
        for id in np.arange(obs.shape[0]-1, 0, -1) :
            o = obs[id]
            r = np.dot(beta[id] * bv[o[2]], a2[o[1]-1].T)
            s = np.sum(r)
            beta[id-1] = r/s
            beta_Pr += np.log(s) + a2x[o[1]-1]
        return alpha_Pr, alpha, beta

    def get_brLens(self, branches, n_base) :
        return np.array([ np.sum(branch.T[1] > 0)/float(n_base) for branch in branches ])

    def prepare_branches(self, branches, interval=None) :
        branches = [ np.array(branch) for branch in branches ]
        for branch in branches :
            branch[branch.T[2] > 1, 2] = 2

        blocks = np.unique([ np.unique(branch.T[0]) for branch in branches ])
        n_bases = np.concatenate([ [[block, np.max(branch[branch.T[0] == block, 1])] for block in blocks] for branch in branches ])
        n_bases = n_bases[np.argsort(-n_bases.T[1])]
        n_bases = n_bases[np.unique(n_bases.T[0], return_index=True)[1]]
        n_bases = n_bases[np.argsort(n_bases.T[0])]
        p = 0
        for nb in n_bases :
            nb[0] = p
            p += nb[1] + 1
        self.n_base = np.sum(n_bases.T[1]) + 1 + len(n_bases)
        if not interval :
            interval = self.n_base
        def prepare_obs(obs, n_bases, interval=None) :
            observation = [[0, 0, 0, 0, 0]]
            obs = np.vstack([obs, [len(n_bases)-1, n_bases[-1, 1]+1, 0]])
            p = 0
            for cc, ss, o in obs :
                s = n_bases[cc, 0] + ss
                d = s - p
                if d == 0 : continue
                cuts = min(4, int(d/interval)+1)
                itv = float(d)/cuts
                px = p
                for c in range(1, cuts+1) :
                    x = int(c*itv)+p if c < cuts else s
                    observation.append([x, x-px, 0 if c < cuts else o, cc, ss])
                    px = x
                p = s
            return np.array(observation)
        self.observations = [prepare_obs(b, n_bases, interval) for b in branches]

    def predict(self, branches=None, names=None, global_rec=False, marginal=0.8, prefix='RecHMM') :
        assert self.model, 'No model'
        import gzip
        stats = self.margin_predict(branches, names, global_rec, marginal) if marginal > 0. and marginal < 1. else self.map_predict(branches, names, global_rec)
        with open(prefix+'.recombination.region', 'wb') as rec_out, gzip.open(prefix+'.mutations.status.gz', 'wb') as mut_out :
            mut_out.write('#Branch\t#Site\t#Type\t#!W[Mutational]\t#!W[Presence]\n')
            for (name, branch) in zip(names, branches) :
                stat = stats[name]
                for site in branch :
                    mut_out.write( '{0}\t{1}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\n'.format(name, site[0], site[1], site[2], stat['stat'].get((site[0], site[1]), 1.), np.sum(stat['weight_p']).astype(float)/stat['weight_p'][0]) )
                rec_out.write('Branch\t{0}\tM={1:.5e}\tR={2:.5e}\tB={3:.3f}\n'.format(name, stat['M'], stat['R'], stat['weight_p'][0]))
                for r in stat['sketches'] :
                    rec_out.write('\tRecomb\t{0}\t{1}\t{2}\t{3}\t{4:.3f}\n'.format(name, r[0], r[1], ['External', 'Internal', 'Mixed   '][r[2]-1], r[3]))

        print 'Imported regions are reported in {0}'.format(prefix+'.recombination.region')
        print 'Mutation status are reported in {0}'.format(prefix+'.mutations.status.gz')

    def map_predict(self, branches=None, names=None, global_rec=False) :
        self.prepare_branches(branches)
        self.screen_out('Predict recombination sketches using', self.model)
        branch_params = self.update_branch_parameters(self.model, lower_limit=True, global_rec=global_rec)
        status = pool.map(functools.partial(_iter_viterbi, self), zip(self.observations, branch_params))
        res = {}
        for name, dm, dr, obs, stat in zip(names, self.model['posterior']['theta'], self.model['posterior']['R'], self.observations, status) :
            rec_len = np.sum([ e-s+1 for s, e, t, p in stat['sketches'] ])
            res[name] = dict(stat=dict(zip(tuple(zip(obs.T[3], obs.T[4])), stat['gamma'])),
                             sketches=stat['sketches'],
                             weight_p=np.array([self.n_base-rec_len, rec_len], dtype=float),
                             M=dm[1]/dm[0],
                             R=np.sum(dr[1:])/dr[0])
        return res

    def margin_predict(self, branches, names=None, global_rec=False, marginal=0.9) :
        self.prepare_branches(branches, 100)
        self.names = np.array(names) if names else np.arange(len(self.observations)).astype(str)
        branch_params = self.update_branch_parameters(self.model, global_rec=global_rec)
        status = self.get_branch_measures(branch_params, self.observations, gammaOnly=True)
        res = {}
        for name, dm, dr, obs, stat in zip(names, self.model['posterior']['theta'], self.model['posterior']['R'], self.observations, status) :
            path = []
            for id, (o, s) in enumerate(zip(obs, stat['gamma'])) :
                p = np.argmax(s)
                if p > 0 and s[0] < 0.5 :
                    if len(path) == 0 or path[-1][2] != p or path[-1][4] != obs[id-1][0] :
                        if o[2] > 0 :
                            path.append([o[0], o[0], p, o[0], o[0], 1-s[0]])
                    else :
                        path[-1][4] = o[0]
                        if 1-s[0] > path[-1][5] :
                            path[-1][5] = 1-s[0]
                        if o[2] > 0 :
                            path[-1][1] = o[0]
            res[name] = dict(stat=dict(zip(tuple(zip(obs.T[3], obs.T[4])), stat['gamma'].T[0])),
                             sketches=[p[:3]+p[5:] for p in path if p[1]-p[0] > 0 and p[5] >= marginal],
                             weight_p=np.sum(stat['b'], 1),
                             M=dm[1]/dm[0],
                             R=np.sum(dr[1:])/dr[0])
        return res


    def viterbi(self, data) :
        obs, params = data
        pi, a, b = params['pi'], params['a'], params['b']
        bv = b.T
        path = np.zeros(shape=[self.n_base, self.n_a], dtype=int)
        alpha = np.zeros(shape=[self.n_base, self.n_a])
        alpha[0] = np.log( np.dot(pi, a) * bv[obs[0, 2]])

        i = 0
        p = np.zeros([self.n_a, self.n_a])
        a[a==0] = 1e-300
        b[b==0] = 1e-300
        pa, pb = np.log(a), np.log(b)
        ids = np.arange(self.n_a)
        for s, d, o, _, _ in obs[1:] :
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
        regions = []
        inrec = np.zeros(path.shape[0])
        for id in np.arange(path.shape[0]-2, 0, -1) :
            max_path = path[id+1, max_path]
            if max_path > 0 :
                inrec[id] = 1
                if len(regions) == 0 or regions[-1][0] != id + 1 :
                    regions.append([id, id, max_path, 1.])
                else :
                    regions[-1][0] = id
        return dict(sketches=sorted(regions), gamma=1.-inrec[ obs.T[0] ])

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
        with open(prefix + '.best.model.report', 'wb') as fout :
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
        print 'Global parameters are summarized in {0}'.format(prefix + '.best.model.report')
        return


def parse_arg(a) :
    parser = argparse.ArgumentParser(description='Parameters for RecHMM. ', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--data', '-d', help='A list of mutations generated by EnPhyl', required=True)
    parser.add_argument('--model', '-m', help='Read a saved model instead of EM process', default='')
    parser.add_argument('--task', '-t', help='task to run. \n0: One rec category from external sources.\n1: Three rec categories considering internal, external and mixed sources [default].', default=1, type=int)
    parser.add_argument('--init', '-i', help='Initiate models with guesses of recombinant proportions. \nDefault: 0.01,0.05,0.3,0.7,0.95,0.99', default='0.01,0.05,0.3,0.7,0.95,0.99')
    parser.add_argument('--prefix', '-p', help='Prefix for all the outputs ', default='RecHMM')
    parser.add_argument('--cool_down', '-c', help='Delete the worst model every N iteration. Default:5', type=int, default=5)
    parser.add_argument('--n_proc', '-n', help='Number of processes. Default: 5. ', type=int, default=5)
    parser.add_argument('--bootstrap', '-b', help='Number of Randomizations for confidence intervals. \nDefault: 1000. ', type=int, default=1000)
    parser.add_argument('--report', '-r', help='Only report the model and do not calculate external sketches. ', default=False, action="store_true")
    parser.add_argument('--marginal', '-M', help='Calculate recombination sketches. \nSet to 0 to use Viterbi algorithm to find most likely path.\n Otherwise (0, 1) use forward-backward algorithm, and report regions with >= M posterior likelihoods as recombinant sketches [Default: 0.8].', default=0.8, type=float)
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
        data = pd.read_csv(args.data, sep='\t', dtype=str, header=0).as_matrix()
        names, branches, blocks = {}, [], {}
        for d in data :
            if d[1] not in blocks :
                blocks[d[1]] = len(blocks)
            if d[0] not in names :
                names[d[0]] = len(names)
                branches.append([])
            id = names[d[0]]
            branches[id].append([ blocks[d[1]], int(d[2]), int(d[3]) ])

        names = [ n for n, id in sorted(names.iteritems(), key=lambda x:x[1]) ]
        ids = [id for id, b in sorted(enumerate(branches), key=lambda x: len(x[1]), reverse=True)]
        branches = [ branches[id] for id in ids ]
        names = [ names[id] for id in ids ]

    if args.model :
        model.load(open(args.model, 'rb'))
    else :
        model.fit(branches, names=names, categories=args.categories, init=args.init, cool_down=args.cool_down)
        model.save(open(args.prefix + '.best.model.json', 'wb'))
        print 'Best HMM model is saved in {0}'.format(args.prefix + '.best.model.json')
    model.report(args.bootstrap, args.prefix)

    if not args.report :
        model.predict(branches, names=names, marginal=args.marginal, prefix=args.prefix)

pool, verbose = None, True
if __name__ == '__main__' :
    RecHMM(sys.argv[1:])

