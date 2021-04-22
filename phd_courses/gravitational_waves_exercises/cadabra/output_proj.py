#!/usr/local/bin/cadabra2
{r,t,\phi,\theta}::Coordinate;
{\mu,\nu,\rho,\sigma,\lambda,\kappa,\chi,\gamma}::Indices(values={t,r,\phi,\theta}, position=fixed);
\partial{#}::PartialDerivative;
g_{\mu\nu}::Metric.
g^{\mu\nu}::InverseMetric.
ss:= { g_{t t} = -(1-2 M/r),
       g_{r r} = 1/(1-2 M/r), 
       g_{\theta\theta} = r**2, 
       g_{\phi\phi}=r**2 \sin(\theta)**2
     }.

complete(ss, $g^{\mu\nu}$);
ch:= \Gamma^{\mu}_{\nu\rho} = 1/2 g^{\mu\sigma} ( 
                                   \partial_{\rho}{g_{\nu\sigma}} 
                                  +\partial_{\nu}{g_{\rho\sigma}}
                                  -\partial_{\sigma}{g_{\nu\rho}} ):
                          
evaluate(ch, ss, rhsonly=True);
rm:= R^{\rho}_{\sigma\mu\nu} = \partial_{\mu}{\Gamma^{\rho}_{\nu\sigma}}
                                  -\partial_{\nu}{\Gamma^{\rho}_{\mu\sigma}}
                                  +\Gamma^{\rho}_{\mu\lambda} \Gamma^{\lambda}_{\nu\sigma}
                                  -\Gamma^{\rho}_{\nu\lambda} \Gamma^{\lambda}_{\mu\sigma};
substitute(rm, ch)
evaluate(rm, ss, rhsonly=True);
acc := a^{\mu} = R^{\mu}_{\nu \rho \sigma} u^{\nu} u^{\rho} \xi^{\sigma};

velocity := {
	u^{t} = 1 / \sqrt{1 - 2 M / r},
	u^{r} = 0,
	u^{\phi} = 0,
	u^{\theta} = 0,
};

displacement := {
	\xi^{t} = 0,
	\xi^{r} = h \sqrt{1 - 2 M / r},
	\xi^{\phi} = 0,
	\xi^{\theta} = 0,
};	

evaluate (acc, ss+velocity+displacement, rhsonly=True);
acc_modulus := a = \sqrt{ a^{\mu} a^{\nu} g_{\mu \nu}};

substitute (acc_modulus, acc);
substitute (acc_modulus, rm);

evaluate(acc_modulus, ss+velocity+displacement, rhsonly=True);