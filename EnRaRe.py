import numpy as np, pandas as pd, sys, os, copy, argparse
import functools, time, datetime
from multiprocessing import Pool

def _iter_branch_measure(obj, arg) :
    return obj.iter_branch_measure(arg)


class RecHMM(object) :
    def __init__(self, mode=1) :
        self.max_iteration = 100
        self.n_base = None
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
        json.dump(model, fout)
    def load(self, fin) :
        import json
        self.model = json.load(fin)
        self.model['p'] = np.array(self.model['p'])
        self.model['pi'] = np.array(self.model['pi'])
        self.model['EventFreq'] = np.array(self.model['EventFreq'])

    def initiate(self, observations, init='0.02,0.2,0.4,0.6,0.8') :
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
                    d, p, po = substitution[0, 0] - 1 + self.n_base - substitution[-1, 0], substitution[0, 0], int(substitution[0, 2] > 0)
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
                         h=[0.5*homos/muts, 2.5*homos/muts],
                         probability=-1e30,
                         diff=1e30,
                         EventFreq=EventFreq,
                         pi = np.zeros(shape=[mut_summary.shape[0], self.n_a]),
                         id = len(self.models) + 1,
                         ite = 0,
                        )
            model['pi'][:] = 1./self.n_a

            rec = np.array(rec)

            if rec.shape[0] :
                rec2, rec = np.sum(rec[rec.T[3]/rec.T[2] > 5.*model['h'][0]], 0), np.sum(rec[rec.T[3]/rec.T[2] <= 5.*model['h'][0]], 0)
                rec[0], rec2[0] = max(rec[0], 1), max(rec2[0], 1)

                model['h'][1] = max(model['h'][1], 0.5*rec2[3]/rec2[2])
                model['delta'], model['v'] = np.min([0.5, (rec[0]+rec2[0])/(rec2[1] + rec[1])]), np.min([0.5, (rec[2]+rec2[2])/(rec2[1] + rec[1])])

                model['p'] = np.array([rec[0], np.sqrt(rec[0]*rec2[0]), rec2[0]])[:(self.n_a-1)]
                model['p'] = model['p']/np.sum(model['p'])
            self.screen_out('Initiate', model)
            self.models.append(model)
        return self.models


    def BaumWelch(self, models, max_iteration, cool_down=5) :
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
                        model['diff'] = 0. #prediction['diff']
                        self.screen_out('Freeze', model)
                        new_models.append(model)
            new_models = sorted(new_models, key=lambda x:-x['probability'])
            if ite > 0 and ite % cool_down == 0 and len(new_models) > 1 :
                self.screen_out('Delete', new_models[-1])
                new_models = new_models[:-1]
            models = new_models
        self.screen_out('Report', models[0])
        return models[0]

    def screen_out(self, action, model) :
        print '{2}\t{0} model {id}[{ite}] - BIC: {3:.8e} [Pr: {probability:.5e}] - EventFreq: {1:.3e}; theta: {theta:.3f}; R: {R:.3f}; Delta: {delta:.3e}; v: {v:.3e}; h: {h[0]:.3f},{h[1]:.3f}; REC sources: {p}'.format(action, np.sum(model['EventFreq']), str(datetime.datetime.now())[:19],  -2*model['probability'] + self.n_a*self.n_b*np.log(self.n_base*model['pi'].shape[0]), **model)
        sys.stdout.flush()

    def estimation(self, model, branch_measures) :
        probability, pi, br = 0., [], []
        a = np.zeros(shape=[self.n_a, self.n_a])
        b = np.zeros(shape=[self.n_a, 3])
        for measures in branch_measures :
            a += measures['a']
            b += measures['b']
            pi.append(measures['pi'])
            br.append([ np.sum(measures['b'][0, 1:])/np.sum(measures['b'][0]), np.sum(measures['a'][0, 1:])/np.sum(measures['a'][0]) ])
            probability += measures['probability']

        prediction = copy.deepcopy(model)
        prediction['probability'] = probability

        prediction['pi'] = np.array(pi)
        prediction['EventFreq'] = np.array(br)

        prediction['R'] = np.sum(a[0, 1:])/np.sum(a[0])*len(branch_measures)/np.sum(br)
        prediction['theta'] = (np.sum(b[0, 1:])+np.sum(b[2:3, 1:]))/(np.sum(b[0])/len(branch_measures) + np.sum(b[2:3]))/np.sum(br)
        prediction['v'] = (np.sum(b[1, 1:]) + np.sum(b[3:, 1:]))/(np.sum(b[1]) + np.sum(b[3:]))

        prediction['p'] = a[0, 1:]/np.sum(a[0, 1:])
        if self.n_a == 2 :
            if self.n_b > 2 :
                prediction['h'] = [np.sum(b[:1, 2])/np.sum(b[:1, 1:]), np.sum(b[1:, 2])/np.sum(b[1:, 1:])]
            else :
                prediction['h'] = [np.sum(b[:2, 2])/np.sum(b[:2, 1:]), np.sum(b[:2, 2])/np.sum(b[:2, 1:])]
        else :
            prediction['h'] = [np.sum(b[:2, 2])/np.sum(b[:2, 1:]), np.sum(b[2:, 2])/np.sum(b[2:, 1:])]

        prediction['delta'] = np.sum(a[1:, 0])/np.sum(a[1:])

        return prediction

    def update_branch_parameters(self, model) :
        branch_params = []
        div_sum = np.sum(model['EventFreq'])
        for d, pi in zip(model['EventFreq'], model['pi']) :
            d = np.sum(d)
            toRec = model['p']*model['R']*d
            a = np.zeros(shape=[self.n_a, self.n_a])
            a[0, 1:] = model['p']*model['R']*d
            a[1:, 0] = model['delta']
            np.fill_diagonal(a, 1-np.sum(a, 1))
            b = np.zeros(shape=[self.n_a, 3])
            b[0] = [1-model['theta']*d, model['theta']*d*(1-model['h'][0]), model['theta']*d*model['h'][0] ]
            if self.n_a == 2 and self.n_b > 2 :
                b[1] = [1-model['v'], model['v']*(1-model['h'][1]), model['v']*model['h'][1] ]
            else :
                b[1] = [1-model['v'], model['v']*(1-model['h'][0]), model['v']*model['h'][0] ]
            if self.n_a > 2 :
                b[2] = [1-model['theta']*div_sum, model['theta']*div_sum*(1-model['h'][1]), model['theta']*div_sum*model['h'][1] ]
                if self.n_a > 3 :
                    b[3] = [1-model['v'], model['v']*(1-model['h'][1]),model['v']*model['h'][1] ]

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
        ot = time.time()

        branch_measures = pool.map(functools.partial(_iter_branch_measure, self), zip(observations, params))
        #print '\t\tTook {0} seconds'.format(time.time() - ot)
        return branch_measures

    def get_interval_drop(self, transition, emission0) :
        p0 = transition[0, 0] * emission0[0]
        return np.array([ (t * e)/p0 for t, e in zip(np.diagonal(transition), emission0) ])

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

        pi2 = (gamma[0] + gamma[-1])/2
        for o, g in zip(obs, gamma) :
            b2[:, o[2]] += g

        for o, s, e in zip(obs[1:], alpha[:-1], beta[1:]) :
            d = o[1]
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

            s1 = np.zeros(shape=[d, self.n_a, self.n_a])
            s1.T[:] = a.T

            s2 = np.zeros(shape=[d, self.n_a, self.n_a])
            s2[:-1].T[:] = (b[1:]*emission.T[0]).T
            s2 = s2.transpose((0, 2, 1))
            s2[-1, :] = e*emission.T[o[2]]

            t = (s1 * s2) * transition.reshape([1] + list(transition.shape))
            a2 += np.sum(t.T/np.sum(t, axis=(1,2)), 2).T

        return dict(a=a2, b=b2, pi=pi2)


    def update_distant_transition(self, transition, emission, interval) :
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
        alpha[0] = pi * bv[obs[0, 2]]
        alpha_Pr = np.sum(alpha[0])
        alpha[0], alpha_Pr = alpha[0]/alpha_Pr, np.log(alpha_Pr)
        for id, o in enumerate(obs[1:]) :
            r = np.dot(alpha[id], a2[o[1]-1]) * bv[o[2]]
            s = np.sum(r)
            alpha[id+1] = r/s
            alpha_Pr += np.log(s) + a2x[o[1]-1]

        beta, beta_Pr = np.ones(shape=[obs.shape[0], self.n_a]), 0.
        beta[-1] = pi
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
        #if self.n_b == 2 :
            #for branch in branches :
                #branch[branch.T[1] > 0, 1] = 1
        #elif self.n_b == 3 :
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
                d, p = s - p, s
                observation.append([s, d, o])
            return np.array(observation)
        self.observations = [prepare_obs(b, self.n_base) for b in branches]

    def predict(self, branches=None, names=None, fout=sys.stdout) :
        assert self.model, 'No model'

        self.prepare_branches(branches)
        self.screen_out('Predict recombination sketches using', self.model)
        branch_params = self.update_branch_parameters(self.model)
        for name, obs, param, (dm, dr) in zip(names, self.observations, branch_params, self.model['EventFreq']) :
            fout.write('Branch\t{0}\tM={1:.5e}\tR={2:.5e}\n'.format( name, -np.log(1-dm*4./3.), -np.log(1-dr) ))
            regions = self.viterbi(obs, pi=param['pi'], a=param['a'], b=param['b'])
            for r in regions :
                fout.write('\tRecomb\t{0}\t{1}\t{2}\t{3}\n'.format(name, *r))


    def viterbi(self, obs, pi, a, b) :
        bv = b.T
        path = np.zeros(shape=[self.n_base, self.n_a], dtype=int)
        alpha = np.zeros(shape=[self.n_base, self.n_a])
        alpha[0] = np.log(pi * bv[obs[0, 2]])

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
        alpha[i] += np.log(pi)
        max_path = np.argmax(alpha[i])
        regions = [] if max_path == 0 else [[i, i, max_path]]
        for id in np.arange(path.shape[0]-2, -1, -1) :
            max_path = path[id+1, max_path]
            if max_path > 0 :
                if len(regions) == 0 or regions[-1][0] != id + 2 :
                    regions.append([id+1, id+1, max_path])
                else :
                    regions[-1][0] = id+1
        return sorted(regions)

def parse_arg() :
    parser = argparse.ArgumentParser(description='Parameters for RecHMM. ', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--data', '-d', help='A list of mutations genereated by EnPhyl', required=True)
    parser.add_argument('--model', '-m', help='Read a saved model instead of EM process', default='')
    parser.add_argument('--task', '-t', help='task to run. \n0: One rec category from external sources.\n1: Three rec categories considering internal, external and mixed sources [default].\n2: One rec category from internal sources.\n3: Two rec categories from internal or external sources.', default=1, type=int)
    parser.add_argument('--init', '-i', help='Initiate models with guesses of recombinant proportions. \nDefault: 0.02,0.2,0.4,0.6,0.8', default='0.02,0.2,0.4,0.6,0.8')
    parser.add_argument('--prefix', '-p', help='Prefix for all the outputs ', default='RecHMM')
    parser.add_argument('--cool_down', '-c', help='Delete the worst model every N iteration. Default:5', type=int, default=5)
    parser.add_argument('--n_proc', '-n', help='Number of processes. Default: 5. ', type=int, default=5)
    args = parser.parse_args()
    return args



if __name__ == '__main__' :
    args = parse_arg()

    data = pd.read_csv(args.data, sep='\t', dtype=str, header=0).as_matrix()
    names, branches = {}, []
    for d in data :
        if d[0] not in names :
            names[d[0]] = len(names)
            branches.append([])
        id = names[d[0]]
        branches[id].append([int(d[2]), int(d[3])])
    model = RecHMM(mode=args.task)
    names = [ n for n, id in sorted(names.iteritems(), key=lambda x:x[1]) ]

    ids = [id for id, b in sorted(enumerate(branches), key=lambda x: len(x[1]), reverse=True)] #[:20]
    branches = [branches[id] for id in ids]
    names = [names[id] for id in ids]

    pool = Pool(args.n_proc)
    if args.model :
        model.load(open(args.model, 'rb'))
    else :
        model.fit(branches, names=names, init=args.init, cool_down=args.cool_down)
        model.save(open(args.prefix + '.best.model.json', 'wb'))

    model.predict(branches, names=names, fout=open(args.prefix + '.recombination.region', 'w'))

