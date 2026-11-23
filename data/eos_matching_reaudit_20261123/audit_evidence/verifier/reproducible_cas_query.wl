Clear[t,s,a0,a1,a2,a3,e0,mu,v];
p=a0+a1 t+a2 t^2+a3 t^3;
aa=p+D[p,t]/s+D[p,{t,2}]/s^2+D[p,{t,3}]/s^3;
cc=e0+(aa/.t->0);
ee=cc Exp[s t]-aa;
gg=ee+p-D[p,t]/s;
bridge=FullSimplify[{(ee/.t->0)-e0,D[ee,t]-s (ee+p),D[gg,{t,2}]-(cc s^2 Exp[s t]-12 a3/s),D[gg,{t,3}]-cc s^3 Exp[s t]},Assumptions->s>0];
a=1201/800;delta=1/1000000;
nn=20 mu-a^3+delta-3 a mu^2+mu^3;
pp=1+Integrate[nn/.mu->v,{v,1,mu}];
en=mu nn-pp;
margin=Expand[mu D[nn,mu]-nn];
<|"bridge_exact_residuals"->bridge,
"counterexample_first_law_residual"->FullSimplify[D[en,mu]-mu D[nn,mu]],
"counterexample_margin"->margin,
"sampled_margin_min_401"->N[Min[Table[margin/.mu->1+k/400,{k,0,400}]],25],
"exact_margin_at_a"->(margin/.mu->a),
"cs2_at_a"->N[(nn/(mu D[nn,mu]))/.mu->a,25],
"n_at_1"->N[nn/.mu->1,25],
"epsilon_at_1"->N[en/.mu->1,25],
"minimum_dn_dmu"->N[D[nn,mu]/.mu->a,25]|>
