# D22SC: D22SG restarted with cold He-rich layers

D22SG (CMFGEN's shock-power width, DVLOC 50 km/s) leaves the He-rich layers
(cells 33-48, 2000-2240 km/s, solar Fe) at 11-13 kK with Fe 99% Fe III and
Mg 97% Mg III, as D22SF did, although its deposition per gram there is 3x
lower. Their temperature is therefore set by the radiation of the deposition
layer, not by local heating. Dessart et al. (2022) find Fe II emission from
exactly these layers.

This run starts from D22SG's final state with T_e = 7000 K in cells 32-50
(Tools/nltex_state_set_te.py). It asks whether a cool, Fe II- and
Mg II-cooled solution exists that a hot warm start never reaches.
