# C14H: the CMFGEN twin of C14

The same model computed by SuperLite and CMFGEN, for a three-way comparison
with SN 2006jc at +20 d.  Made with

    Tools/stella_twin_homologous.py Runs/iin_day20/C14 Runs/iin_day20/C14H

## How it differs from C14

- **Homologous at 30.37 d** (the median r/v of C14).  C14's dense shell is
  not homologous: its velocity falls and rises again with radius, and
  CMFGEN needs a monotonic field.  Each C14 cell's mass, element masses and
  power are spread over the velocity range it spans; radii become r = v t.
  The mass per velocity interval, which sets the line profiles, is kept.
  The shell's density peak sits at 8070-8330 km/s.
- **All power as heat**: 7.69e42 erg/s (STELLA's deposition and the shock
  luminosity), in SuperLite with `in_nt_lumg = 0` and no shock cells, in
  CMFGEN as a prescribed deposition with the same profile (`EDEP_PROFILE`).
  C14 sends the deposition to fast electrons and radiates 3.0e42 erg/s of
  it as CIE shock spectra.
- **H, Ar and Ni left out** (no model atoms in CMFGEN's set), the rest
  renormalized.  The ions are IBNH's, the same in both codes.

## CMFGEN

`VADAT` is IBNH's corrected one (`Runs/iin_day20/IBNH/cmfgen/README.md`)
with this model's radius, velocity, mass, age and power.  `SN_HYDRO_DATA`
and `EDEP_PROFILE` come from the script.  CMFGEN's shock power has a
Gaussian profile only, so the deposition table needs a local patch of
`new_main/subs/get_edep_shock_power.f`: when a file `EDEP_PROFILE`
exists (comment lines starting with `!`, the number of points, then
v [km/s] and power per volume), the deposition follows it, linearly
interpolated in v, instead of the Gaussian, still normalized to
`PRESCRIBED_SHOCK_POWER`.  Without the file CMFGEN is unchanged.
