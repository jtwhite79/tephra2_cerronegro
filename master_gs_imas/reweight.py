import pst_handler as phand

pst = phand.pst("base.pst")
pst.proportional_weights(0.25,0.25)
obs = pst.observation_data
obs.index = obs.obgnme
obs = obs.sort(["obgnme","obsval"])
pst.observation_data = obs
pst.zero_order_tikhonov()
#pst.parameter_data.parlbnd = 1.0e-20
#pst.parameter_data.parubnd = 1.0e+20
#grp = pst.parameter_data.groupby("partrans").groups
#pst.parameter_data.loc[grp["none"],"parlbnd"] = -1.0E20
pst.write("reg_base.pst",True)

pst.proportional_weights(0.4,0.25)
obs = pst.observation_data
obs.index = obs.obgnme
imas_grp = obs.groupby(lambda x : "imas" in x.lower()).groups
obs.weight.loc[imas_grp[False]] = 0.0
obs = obs.sort(["obgnme","obsval"])
pst.observation_data = obs
pst.write("reg_base_imas.pst",True)
