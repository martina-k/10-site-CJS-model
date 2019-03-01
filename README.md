# 10-site-model
Metapopulation model for Stockholm archipelago data, with 10 islands, and, in some models, an Elsewhere site (all other sites)  

Age is used as state, with with 5 stages: 1-year,2-year,3-year,4-year, adults.  
Adults are in turn separated into recently ringed and previously ringed.  
Birds are only ringed as chicks and adults (recently ringed).

Models with numbers work well (i.e. compile, run and produce results)

Models xny and x do not work well - they give 'invalid parent error' during most attempts to run them

Model 11ny produce reasonable estimates for: 
p at islands 1-8
phi for 1-year-olds, and for 2-4-years-olds (estimated jointly)

Problems with 11ny:
p for islands 9-10 is somewhat difficult to interpret conceptually
phi for adults is low (~ 84%, should be >90%), but could be due to heterogeneity
phi at island 9 is low, whereas it is not in model x runs

Therefore it would be great to make model xny work, to combine the advantages of model 11ny and model x.
So far I have received 'invalid parent error' every time.