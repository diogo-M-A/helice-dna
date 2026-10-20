"""Projeção de mistura com frequências fixas e restrição ao simplex.

NumPy é opcional no app geral e necessário para painéis grandes/reamostragem.
Nenhum componente é interpretado como identidade étnica ou nacionalidade.
"""
import numpy as np


def project_simplex(v):
    values=np.sort(v)[::-1]
    cumulative=np.cumsum(values)-1
    valid=np.flatnonzero(values-cumulative/np.arange(1,len(v)+1)>0)
    theta=cumulative[valid[-1]]/(valid[-1]+1)
    return np.maximum(v-theta,0)


def fit_mixture(frequencies,dosages,weights=None,initial=None,tolerance=1e-8,max_iter=1000):
    """Maximiza a log-verossimilhança binomial, côncava em q.

Newton em conjunto ativo + busca de passo; critério de parada é o gap KKT.
Retorna resíduos para não confundir um limite de iterações com convergência.
"""
    f=np.asarray(frequencies,dtype=np.float64)
    g=np.asarray(dosages,dtype=np.float64)
    if f.ndim!=2 or g.shape!=(len(f),) or not np.isfinite(f).all() or not np.isfinite(g).all():
        raise ValueError('Matrizes de projeção inválidas.')
    if len(f)==0 or f.shape[1]<2 or (f<0).any() or (f>1).any() or (g<0).any() or (g>2).any():
        raise ValueError('Frequências ou dosagens fora do intervalo permitido.')
    w=np.ones(len(g)) if weights is None else np.asarray(weights,dtype=np.float64)
    if w.shape!=g.shape or not np.isfinite(w).all() or (w<0).any() or w.sum()<=0:
        raise ValueError('Pesos de reamostragem inválidos.')
    w=w/w.sum()
    k=f.shape[1]
    q=np.full(k,1/k) if initial is None else project_simplex(np.asarray(initial,dtype=float))

    def objective(x):
        p=np.clip(f@x,1e-10,1-1e-10)
        return float(w@(g*np.log(p)+(2-g)*np.log1p(-p)))

    value=objective(q)
    converged=False
    for iteration in range(max_iter):
        p=np.clip(f@q,1e-10,1-1e-10)
        gradient=f.T@(w*(g/p-(2-g)/(1-p)))
        gap=float(gradient.max()-gradient@q)
        if gap<tolerance:
            converged=True
            break
        curvature=w*(g/p**2+(2-g)/(1-p)**2)
        hessian=f.T@(curvature[:,None]*f)
        active=set(np.flatnonzero(q>1e-12).tolist())
        active.add(int(np.argmax(gradient)))
        direction=None
        for _ in range(k):
            indices=sorted(active)
            matrix=np.zeros((len(indices)+1,len(indices)+1))
            matrix[:-1,:-1]=hessian[np.ix_(indices,indices)]+np.eye(len(indices))*1e-10
            matrix[-1,:-1]=1
            matrix[:-1,-1]=1
            rhs=np.r_[gradient[indices],0.0]
            try: solution=np.linalg.solve(matrix,rhs)[:-1]
            except np.linalg.LinAlgError: break
            negative_zero=[indices[i] for i,d in enumerate(solution) if d<0 and q[indices[i]]<=1e-12]
            if negative_zero:
                active.difference_update(negative_zero)
                continue
            direction=np.zeros(k)
            direction[indices]=solution
            break
        if direction is None or float(gradient@direction)<=1e-14:
            scale=max(float(np.linalg.eigvalsh(hessian)[-1]),1)
            direction=project_simplex(q+gradient/scale)-q
        negative=direction<0
        alpha=min(1.0,float(np.min(-q[negative]/direction[negative]))) if negative.any() else 1.0
        improvement=float(gradient@direction)
        accepted=False
        for _ in range(60):
            candidate=np.maximum(q+alpha*direction,0)
            candidate/=candidate.sum()
            candidate_value=objective(candidate)
            if candidate_value>=value+1e-4*alpha*improvement-1e-14:
                q=candidate
                value=candidate_value
                accepted=True
                break
            alpha*=.5
        if not accepted:
            break
    p=np.clip(f@q,1e-10,1-1e-10)
    gradient=f.T@(w*(g/p-(2-g)/(1-p)))
    gap=float(gradient.max()-gradient@q)
    converged=gap<tolerance
    return {'q':q,'converged':converged,'iterations':iteration+1,'kkt_gap':gap,'log_likelihood_per_marker':value}


def bootstrap_blocks(frequencies,dosages,blocks,replicates=200,seed=20261019,initial=None,progress=lambda *args:None):
    """Reamostra blocos genômicos, mantendo o painel de frequências fixo."""
    f=np.asarray(frequencies,dtype=float)
    g=np.asarray(dosages,dtype=float)
    names=sorted(set(blocks))
    mapping={name:i for i,name in enumerate(names)}
    block_index=np.asarray([mapping[b] for b in blocks])
    rng=np.random.default_rng(seed)
    values=[]
    failures=0
    for i in range(replicates):
        multiplicities=np.bincount(rng.integers(0,len(names),len(names)),minlength=len(names))
        fit=fit_mixture(f,g,weights=multiplicities[block_index],initial=initial,tolerance=1e-7)
        if not fit['converged']: failures+=1
        values.append(fit['q'])
        if (i+1)%20==0:progress('Reamostrando blocos genômicos',i+1)
    values=np.asarray(values)
    return {'replicates':replicates,'seed':seed,'blocks':len(names),'nonconverged':failures,
            'lower':np.percentile(values,2.5,axis=0),'upper':np.percentile(values,97.5,axis=0),'samples':values}
