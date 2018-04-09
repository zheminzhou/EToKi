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

    def fit(self, branches, names=None, brLens=None, init=None, cool_down=5) :
        self.prepare_branches(branches)
        self.brLens = np.array(brLens) if brLens else None
        self.names = np.array(names) if names else np.arange(len(self.observations)).astype(str)

        models = self.initiate(self.observations, init=init)
        return self.BaumWelch(models, self.max_iteration, cool_down=cool_down)

    def save(self, fout):
        import json
        model = copy.deepcopy(self.model)
        model['p'] = model['p'].tolist()
        model['pi'] = model['pi'].tolist()
        model['EventFreq'] = model['EventFreq'].tolist()
        if 'posterior' in model :
            for b in model['posterior'] :
                b['pi'] = b['pi'].tolist()
                b['a'] = b['a'].tolist()
                b['b'] = b['b'].tolist()
        json.dump(model, fout)
        return self.model
    
    def load(self, fin) :
        import json
        self.model = json.load(fin)
        self.model['p'] = np.array(self.model['p'])
        self.model['pi'] = np.array(self.model['pi'])
        self.model['EventFreq'] = np.array(self.model['EventFreq'])
        if 'posterior' in self.model :
            self.n_base = int(0.5+np.sum(self.model['posterior'][0]['b']))
            for b in self.model['posterior'] :
                b['pi'] = np.array(b['pi'])
                b['a'] = np.array(b['a'])
                b['b'] = np.array(b['b'])
        self.n_a = self.model['p'].size + 1
        if self.n_a == 2 and self.model['h'][0] == self.model['h'][1] :
            self.n_b = 2
        else :
            self.n_b = 3
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
                    for s, d, o in substitution[1:] :
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

            EventFreq = np.array([3./4.*(1-np.exp(-self.brLens)), 0.]) if self.brLens else (mut_summary.T[[1, 3]]/mut_summary.T[0]).T

            n_br, (bases, muts, homos, recs) = mut_summary.shape[0], np.sum(mut_summary, 0)
            model = dict(theta=muts/np.sum(EventFreq)/self.n_base,
                         R=recs/np.sum(EventFreq)/self.n_base,
                         h=[max(0.01, min(0.95, 1./2.*homos/muts)), max(0.01, min(0.95, 5./2.*homos/muts))], 
                         probability=-1e30,
                         diff=1e30,
                         EventFreq=EventFreq,
                         pi = np.zeros(shape=[mut_summary.shape[0], self.n_a]),
                         id = len(self.models) + 1,
                         ite = 0,
                        )
            model['pi'][:, 0] = 1

            rec = np.array(rec)

            if rec.shape[0] :
                rec2, rec = np.sum(rec[rec.T[3]/rec.T[2] > 2.*2.5*model['h'][0]], 0), np.sum(rec[rec.T[3]/rec.T[2] <= 2.*2.5*model['h'][0]], 0)
                rec[0], rec2[0] = max(rec[0], 1), max(rec2[0], 1)
                
                model['h'][1] = max(model['h'][1], 1./2.*(rec2[3]+1)/(rec2[2]+1))
                model['delta'], model['v'], model['v2'] = max(0.00002, min(0.005, (rec[0]+rec2[0])/(rec2[1] + rec[1]))), np.max([0.001, np.min([0.5, (rec[2]+rec2[2])/(rec2[1] + rec[1])])]), np.min([0.5, (rec2[2]+1)/(rec2[1]+1)])

                model['p'] = np.array([rec[0], np.sqrt(rec[0]*rec2[0]), rec2[0]])[:(self.n_a-1)]
                model['p'] = model['p']/np.sum(model['p'])
            n_total = model['theta'] + model['R']
            model['theta'], model['R'] = model['theta']/n_total, model['R']/n_total
            if self.brLens :
                EventFreq[:] = EventFreq.T[0] * np.array([model['theta'], model['R']])
            self.screen_out('Initiate', model)
            self.models.append(model)
        return self.models


    def BaumWelch(self, models, max_iteration, cool_down=5) :
        n_model = len(models)
        for id, model in enumerate(models) :
            if 'id' not in model :
                model['id']
            if 'ite' not in model :
                model['ite'] = id + 1
        for ite in np.arange(max_iteration) :
            new_models = []
            self.model = models[0]

            for model in models:
                if 'diff' in model and model['diff'] < 0.001 :
                    new_models.append(model)
                else :
                    print ''
                    self.screen_out('Assess', model)
                    branch_params = self.update_branch_parameters(model)
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
            models = new_models
        self.screen_out('Report', models[0])
        return models[0]

    def screen_out(self, action, model) :
        print '{2}\t{0} model {id}[{ite}] - BIC: {3:.8e} - EventFreq: {1:.3e}; theta: {theta:.3f}; R: {R:.3f}; delta: {4:.3e};  Nu: {v:.3e},{v2:.3e}; h: {h[0]:.3f},{h[1]:.3f}; REC sources: {p}'.format(action, np.sum(model['EventFreq']), str(datetime.datetime.now())[:19],  -2*model['probability'] + self.n_a*self.n_b*np.log(self.n_base*model['pi'].shape[0]), 1/model['delta'], **model)
        sys.stdout.flush()

    def estimation(self, model, branch_measures) :
        h = [1-model['h'][0]/2, np.sqrt(model['h'][0]*(2-model['h'][0]))/2.]
        
        probability, pi, br = 0., [], []
        a = np.zeros(shape=[self.n_a, self.n_a])
        b = np.zeros(shape=[self.n_a, 3])
        
        for id, measures in enumerate(branch_measures) :
            a += measures['a']
            b += measures['b']
            pi.append(measures['pi'])
            if self.n_a == 2 and self.n_b > 2 :
                br.append([ ( np.sum(measures['b'][0, 1:]) + np.sum(measures['b'][1:2, 1]) )/( np.sum(measures['b'][0]) + np.sum(measures['b'][1:2])*h[0] ), np.sum(measures['a'][0, 1:])/np.sum(measures['a'][0]) ])
            else :
                br.append([ ( np.sum(measures['b'][0, 1:]) + np.sum(measures['b'][2:3, 1]) )/( np.sum(measures['b'][0]) + np.sum(measures['b'][2:3])*h[0] ), np.sum(measures['a'][0, 1:])/np.sum(measures['a'][0]) ])
            if br[-1][0] < .1/self.n_base :
                br[-1][0] = .1/self.n_base
            if br[-1][0]*0.001 > br[-1][1] :
                br[-1][1] = br[-1][0]*0.001
            if np.sum(br[-1]) < 1./self.n_base :
                br[-1] = [br[-1][0]/self.n_base/np.sum(br), br[-1][1]/self.n_base/np.sum(br)]
            probability += measures['probability']

        prediction = copy.deepcopy(model)
        prediction['posterior'] = branch_measures
        prediction['probability'] = probability

        prediction['pi'] = np.array(pi)
        prediction['EventFreq'] = np.array(br)
        prediction['theta'] = np.sum(prediction['EventFreq'].T[0]/np.sum(prediction['EventFreq'], 1))
        prediction['R'] = np.sum(prediction['EventFreq'].T[1]/np.sum(prediction['EventFreq'], 1))
        
        n_br, sum_br = prediction['EventFreq'].shape[0], np.sum(prediction['EventFreq'])
        prediction['p'] = a[0, 1:]/np.sum(a[0, 1:])
        prediction['delta'] = np.sum(a[1:, 0])/np.sum(a[1:])
        
        if self.n_a == 2 and self.n_b == 2 :
            prediction['v'] = np.sum(b[1, 1:])/np.sum(b[1])
            prediction['v2'] = np.sum(b[1, 2:])/np.sum(b[1])
            prediction['h'] = [np.sum(b[:2, 2])/np.sum(b[:2, 1:]), np.sum(b[:2, 2])/np.sum(b[:2, 1:])]
        elif self.n_a == 2 and self.n_b > 2 :
            prediction['v'] = np.sum(b[1, 1:])/np.sum(b[1])
            prediction['v2'] = np.sum(b[1, 2:])/np.sum(b[1])
            prediction['h'] = [np.sum(b[:1, 2])/np.sum(b[:1, 1:]), np.sum(b[1:, 2])/np.sum(b[1:, 1:])]
            
        elif self.n_a == 3 :
            prediction['v'] = np.sum(b[1, 1:])/np.sum(b[1])
            prediction['v2'] = np.sum(b[2:, 2:])/np.sum(b[2:])
            prediction['h'] = [np.sum(b[:2, 2])/np.sum(b[:2, 1:]), np.sum(b[2:, 2])/np.sum(b[2:, 1:])]
        else :
            prediction['v'] = ( np.sum(b[1, 1:]) + np.sum(b[3, 1]) )/( np.sum(b[1]) + np.sum(b[3])*h[0] )
            prediction['v2'] = np.sum(b[2:, 2:])/np.sum(b[2:])
            prediction['h'] = [np.sum(b[:2, 2])/np.sum(b[:2, 1:]), np.sum(b[2:, 2])/np.sum(b[2:, 1:])]
        for i, px in enumerate(prediction['p']) :
            if px < 1e-300 :
                prediction['p'][i] = 1e-300
        n_total = prediction['theta'] + prediction['R']
        prediction['theta'], prediction['R'] = prediction['theta']/n_total, prediction['R']/n_total
        max_br = np.max(np.array(br).T[0])
        if prediction['delta'] > 0.01 or prediction['delta'] < 0.00001 :
            prediction['delta'] = min(max(model['delta'], 0.00001), 0.01)
        if prediction['v'] < 0.0005:
            prediction['v'] = max( 0.0005, model['v'] )
        if self.n_b > 2 and prediction['v2'] < max_br/2. :
            prediction['v2'] = max( max_br/2., model['v2'] )
        return prediction

    def update_branch_parameters(self, model, lower_limit=False) :
        branch_params = []
        
        for d, pi in zip(model['EventFreq'], model['pi']) :
            if lower_limit and np.sum(d) < 1./self.n_base:
                d[:] = d/self.n_base/np.sum(d)
                
            a = np.zeros(shape=[self.n_a, self.n_a])
            a[0, 1:] = model['p']*d[1]
            a[1:, 0] = model['delta']
            np.fill_diagonal(a, 1-np.sum(a, 1))
            b = np.zeros(shape=[self.n_a, 3])
            h = [1-model['h'][0]/2, np.sqrt(model['h'][0]*(2-model['h'][0]))/2.]
            
            mut = d[0]
            
            b[0]  = [ 1-mut,                         mut*h[0],        mut*h[1]        ]
            extra = [ 1-model['v'],                  model['v']*h[0], model['v']*h[1] ]
            intra = [ 1-mut*h[0] - model['v2'],      mut*h[0],        model['v2']     ]
            mixed = [ 1-model['v']*h[0]-model['v2'], model['v']*h[0], model['v2']     ]
            if self.n_a == 2 and self.n_b > 2 :
                b[1] = intra
            else :
                b[1] = extra
            if self.n_a > 2 :
                b[2] = intra
                if self.n_a > 3 :
                    b[3] = mixed

            branch_params.append(dict(
                pi = pi,
                a = a,
                b = b,
            ))
        return branch_params

    def iter_branch_measure(self, data) :
        obs, param = data
        a2, a2x, saturate_id = self.update_distant_transition(param['a'], param['b'].T, np.max(obs.T[1]))
        alpha_Pr, alpha, beta = self.forward_backward(obs, pi=param['pi'], a2s=[a2, a2x], b=param['b'])
        new_param = self.estimate_params(param['a'], param['b'], obs, alpha, beta, a2, saturate_id)
        new_param['probability'] = alpha_Pr
        return new_param

    def get_branch_measures(self, params, observations) :
        branch_measures = pool.map(functools.partial(_iter_branch_measure, self), zip(observations, params))
        return branch_measures

    def estimate_params(self, transition, emission, obs, alpha, beta, tr2, saturate_id, mode='accurate') : # mode = accurate
        gamma = alpha*beta
        gamma = (gamma.T/np.sum(gamma, 1)).T

        na = np.dot(alpha[0], tr2[saturate_id])*emission.T[0]
        nb = np.dot(beta[0], tr2[saturate_id].T)
        ng = na*nb/np.sum(na*nb)
        ne = np.array([na for i in np.arange(self.n_a)]).T * np.array([nb*emission.T[0] for i in np.arange(self.n_a)])*transition
        ne = ne/np.sum(ne)

        a2 = np.zeros(shape=[self.n_a, self.n_a])
        b2 = np.zeros(shape=[self.n_a, 3])

        pi2 = np.zeros(self.n_a)
        pi2[0] = 1.
        
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
        return dict(a=a2, b=b2, pi=pi2)


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

    def prepare_branches(self, branches) :
        branches = [ np.array(branch) for branch in branches ]
        for branch in branches :
            branch[branch.T[1] > 1, 1] = 2
        if self.n_base is None :
            self.n_base = np.max([ np.max(branch.T[0]) for branch in branches ]) + 100
        def prepare_obs(obs, n_base) :
            observation = [[1, 0, 0]] if obs[0][0] != 1 else []
            if obs.T[0, -1] < n_base :
                obs = np.vstack([obs, [n_base, 0]])
            p = 1
            for s, o in obs :
                d = s - p
                observation.append([s, d, o])
                p = s
            return np.array(observation)
        self.observations = [prepare_obs(b, self.n_base) for b in branches]

    def predict(self, branches=None, names=None, fout=sys.stdout) :
        assert self.model, 'No model'

        self.prepare_branches(branches)
        self.screen_out('Predict recombination sketches using', self.model)
        branch_params = self.update_branch_parameters(self.model, lower_limit=True)
        regions = pool.map(functools.partial(_iter_viterbi, self), zip(self.observations, branch_params))
        res = {}
        for name, (dm, dr), obs, (region, snp_status) in zip(names, self.model['EventFreq'], self.observations, regions) :
            res[name] = dict(zip(obs.T[0], snp_status))
            fout.write('Branch\t{0}\tM={1:.5e}\tR={2:.5e}\n'.format( name, dm, dr ))
            for r in region :
                fout.write('\tRecomb\t{0}\t{1}\t{2}\t{3}\n'.format(name, r[0], r[1], ['External', 'Internal', 'Mixed'][r[2]-1]))
            fout.flush()
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
        for s, d, o in obs[1:] :
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
        regions = [] if max_path == 0 else [[i+1, i+1, max_path]]
        inrec = np.zeros(path.shape[0])
        for id in np.arange(path.shape[0]-2, -1, -1) :
            max_path = path[id+1, max_path]
            if max_path > 0 :
                inrec[id] = 1
                if len(regions) == 0 or regions[-1][0] != id + 2 :
                    regions.append([id+1, id+1, max_path])
                else :
                    regions[-1][0] = id+1
        return sorted(regions), inrec[obs.T[0]-1]
    
    def get_parameter(self, bootstrap, prefix='') :
        if 'posterior' in self.model :
            posterior = self.model['posterior']
        else :
            posterior = []
        n_br = len(posterior)
        bs = np.random.randint(n_br, size=[bootstrap, n_br])
        
        aa = np.array([d['a'] for d in posterior])[bs]
        b = np.array([d['b'] for d in posterior])
        bb = b[bs]
        a_sum = np.sum([d['a'] for d in posterior], 0)
        b_sum = np.sum(b, 0)
        
        pp = np.array([d['probability'] for d in posterior])[bs]

        h0 = 1-self.model['h'][0]/2

        reports = {}
        
        reports['probability'] = [self.model['probability'], 
                                  np.sum(pp, 1)]
        reports['BIC'] = [-2*reports['probability'][0] + self.n_a*self.n_b*np.log(self.n_base*n_br), -2*reports['probability'][1] + self.n_a*self.n_b*np.log(self.n_base*n_br) ]
        reports['p'] = [self.model['p'], 
                        np.sum(aa, 1)[:, 0, 1:].T/np.sum(np.sum(aa, 1)[:, 0, 1:], 1)]

        if self.n_a == 2 and self.n_b > 2 :
            reports['EventFreq'] = [np.sum(self.model['EventFreq']), 
                                    np.array([(np.sum(bb[:, :, 0, 1:], 2)+bb[:,:,1, 1])/(np.sum(bb[:, :, 0, :], 2)+np.sum(bb[:, :, 1, :], 2)*h0), np.sum(aa[:, :, 0, 1:], 2)/np.sum(aa[:, :, 0, :], 2)])]#.transpose((1,2,0))]
        else :
            reports['EventFreq'] = [np.sum(self.model['EventFreq']), 
                                    np.array([(np.sum(bb[:, :, 0, 1:], 2)+bb[:,:,2, 1])/(np.sum(bb[:, :, 0, :], 2)+np.sum(bb[:, :, 2, :], 2)*h0), np.sum(aa[:, :, 0, 1:], 2)/np.sum(aa[:, :, 0, :], 2)])]#.transpose((1,2,0))]
        x = np.sum(reports['EventFreq'][1]/np.sum(reports['EventFreq'][1], 0), 2)
        reports['R'] = [np.sum(self.model['EventFreq'].T[1]/np.sum(self.model['EventFreq'], 1)), x[1]]
        reports['theta'] = [np.sum(self.model['EventFreq'].T[0]/np.sum(self.model['EventFreq'], 1)), x[0]]
        reports['EventFreq'][1] = np.sum(reports['EventFreq'][1], (0,2))        
        
        reports['delta'] = [1/self.model['delta'], 
                            np.sum(np.sum(aa, 1)[:, 1:, 0], 1)/np.sum(np.sum(aa, 1)[:, 1:, :], (1,2))]
        reports['delta'][1][reports['delta'][1] < 1e-5] = 1e-5
        reports['delta'][1] = 1./reports['delta'][1]

        if self.n_a == 2 and self.n_b == 2 :
            nu    = ( np.sum(np.sum(bb, 1)[:, 1, 1:], 1) )/( np.sum(np.sum(bb, 1)[:, 1, :], 1) )
            nu2   = ( np.sum(np.sum(bb, 1)[:, 1, 2:], 1) )/( np.sum(np.sum(bb, 1)[:, 1, :], 1) )
            h     = [ np.sum(b[:, :2, 2], 1)/np.sum(b[:, :2, 1:], (1, 2)), 
                      np.sum(b[:, :2, 2], 1)/np.sum(b[:, :2, 1:], (1, 2)) ]
        elif self.n_a == 2 and self.n_b > 2 :
            nu    = ( np.sum(np.sum(bb, 1)[:, 1, 1:], 1) )/( np.sum(np.sum(bb, 1)[:, 1, :], 1) )
            nu2   = ( np.sum(np.sum(bb, 1)[:, 1, 2:], 1) )/( np.sum(np.sum(bb, 1)[:, 1, :], 1) )
            h     = [ np.sum(b[:, :1, 2], 1)/np.sum(b[:, :1, 1:], (1, 2)), 
                      np.sum(b[:, 1:, 2], 1)/np.sum(b[:, 1:, 1:], (1, 2)) ]
        elif self.n_a == 3 :
            nu    = ( np.sum(np.sum(bb, 1)[:, 1, 1:], 1) )/( np.sum(np.sum(bb, 1)[:, 1, :], 1) )
            nu2   = ( np.sum(np.sum(bb, 1)[:, 2, 2:], 1) )/( np.sum(np.sum(bb, 1)[:, 2, :], 1) )
            h     = [ np.sum(b[:, :2, 2], 1)/np.sum(b[:, :2, 1:], (1, 2)), 
                      np.sum(b[:, 2:, 2], 1)/np.sum(b[:, 2:, 1:], (1, 2)) ]
        else :
            nu    = ( np.sum(np.sum(bb, 1)[:, 1, 1:], 1) + np.sum(np.sum(bb, 1)[:, 3, 1:2], 1) )/( np.sum(np.sum(bb, 1)[:, 1, :], 1) + np.sum(np.sum(bb, 1)[:, 3, :], 1)*h0 )
            nu2   = ( np.sum(np.sum(bb, 1)[:, 2:, 2:], (1,2)) )/( np.sum(np.sum(bb, 1)[:, 2:, :], (1,2)) )
            h     = [ np.sum(b[:, :2, 2], 1)/np.sum(b[:, :2, 1:], (1, 2)), 
                      np.sum(b[:, 2:, 2], 1)/np.sum(b[:, 2:, 1:], (1, 2)) ]
        n_total = [reports['theta'][0] + reports['R'][0], reports['theta'][1] + reports['R'][1]]
        reports['R']     = [reports['R'][0]/n_total[0],     reports['R'][1]/n_total[1]]
        reports['theta'] = [reports['theta'][0]/n_total[0], reports['theta'][1]/n_total[1]]
        if reports['R'][0] < 0.001 :
            reports['R'][0] = 0.001
            reports['theta'][0] = 0.999
        
        reports['nu'] = [(np.sum(b_sum[1, 1:]) + np.sum(b_sum[3, 1]))/(np.sum(b_sum[1, :])+np.sum(b_sum[3, :])*h0), nu]
        reports['nu(in)'] = [np.sum(b_sum[2:, 2:])/np.sum(b_sum[2:, :]), nu2]
        reports['h'] = [self.model['h'], h]

        reports['R/theta'] = [reports['R'][0]/reports['theta'][0], reports['R'][1]/reports['theta'][1]]
        
        s0 = np.sum(np.sum(b, 0)[:, 1:], 1)
        s1 = np.sum(np.sum(bb, 1)[:, :, 1:], 2)
        reports['r/m'] = [np.sum(s0[1:])/s0[0], np.sum(s1[:, 1:], 1)/s1[:, 0]]
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
    args = parser.parse_args(a)
    return args

def RecHMM(args) :
    args = parse_arg(args)
    global pool
    pool = Pool(args.n_proc)
    
    model = recHMM(mode=args.task)
    if not args.report or not args.model :
        data = pd.read_csv(args.data, sep='\t', dtype=str, header=0).as_matrix()
        names, branches = {}, []
        for d in data :
            if d[0] not in names :
                names[d[0]] = len(names)
                branches.append([])
            id = names[d[0]]
            branches[id].append([int(d[2]), int(d[3])])
        
        names = [ n for n, id in sorted(names.iteritems(), key=lambda x:x[1]) ]
    
        ids = [id for id, b in sorted(enumerate(branches), key=lambda x: len(x[1]), reverse=True)]
        branches = [branches[id] for id in ids]
        names = [names[id] for id in ids]

    if args.model :
        model.load(open(args.model, 'rb'))
    else :
        model.fit(branches, names=names, init=args.init, cool_down=args.cool_down)
        model.save(open(args.prefix + '.best.model.json', 'wb'))
        print 'Best HMM model is saved in {0}'.format(args.prefix + '.best.model.json')
    model.get_parameter(args.bootstrap, args.prefix)
    
    if not args.report :
        region_out = args.prefix + '.recombination.region'
        snp_status = model.predict(branches, names=names, fout=open(region_out, 'w'))
        print 'Imported regions are reported in {0}'.format(region_out)
        import gzip
        with gzip.open(args.prefix+'.mutations.status.gz', 'wb') as fout :
            for d in data :
                fout.write('{0}\t{1}\n'.format('\t'.join(d), snp_status[d[0]][int(d[2])]))
        print 'Mutation status are reported in {0}'.format(args.prefix+'.mutations.status.gz')

pool = None
if __name__ == '__main__' :
    RecHMM(sys.argv[1:])