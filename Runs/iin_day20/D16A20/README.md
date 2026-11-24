# D16A20: twin of model A of Dessart et al. (2016) at 20 d

Model A of Dessart, Hillier, Audit, Livne & Waldman (2016, MNRAS 458, 2094):
a ~10 Msun, 1e51 erg ejecta ramming into a 0.4 Msun CSM (100 km/s), H-rich,
the analogue of SN 1998S.  Their CMFGEN spectra post-process HERACLES
snapshots with the temperature held fixed (their Sect. 2), at 2.0, 11.6,
20.0, 34.7, 57.9 and 81.0 d after the onset of interaction (Figs 5, 6).

- Structure at 20.0 d: T, V, rho and the radiative luminosity of the 8
  energy groups against radius from their Fig. 3 (right column), read from
  the vector paths of the PDF (`modelA_20d_raw.npz`; Tools/dessart/
  pdf_vector_paths.py, axes calibrated on the tick marks). Checks: total
  mass 10.44 Msun (ejecta ~10 + CSM 0.4), L_bol = sum of the groups at the
  outer edge 2.38e43 erg/s = 6.2e9 Lsun (their Fig. 2 at 20 d: ~10^9.8).
- `input.str` by `Tools/dessart/build_modelA.py`: 120 cells refined where
  rho and V change fast; the non-monotonic velocity (homologous ejecta to
  5600 km/s, the CDS at 1.55e15 cm, the CSM at 110-220 km/s); temp = their T,
  held fixed (`in_tesolve = f`); lum = L(r), the sum of the groups
  (`lumg`, as heat, `in_nt_lumg = 0`); shock cells off (`in_shock_dvmin`
  1e12: CMFGEN had only the hydro temperatures there).
- Composition: uniform, H 0.633 and He 0.36564 as their HERACLES runs,
  with solar metals (Asplund et al. 2009; not given in the paper) scaled
  into the rest. CMFGEN also had Sc II and Ti II (no SuperLite atoms).
- The comparison: their Halpha profiles of model A (Fig. 6, vector:
  `modelA_halpha_cmfgen.npz`, raw plot coordinates; profiles normalized and
  shifted, compared by shape) and their optical spectrum (Fig. 5, a raster
  image).
