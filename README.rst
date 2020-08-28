Code for a Spitzer MIR Extinction Paper
=======================================

Routines for Spitzer MIR dust extinction curve paper.
This paper is Gordon, Misselt, et al. (in prep).

In Development!
---------------

Active development.
Data still in changing.
Use at your own risk.

Contributors
------------
Karl Gordon

License
-------

This code is licensed under a 3-clause BSD style license (see the
``LICENSE`` file).

Extinction Curves
-----------------

External packages used include `dust_extintion`, `measure_extinction`, and
`extstar_data`.  The basic information on each star is given in `extstar_data`
including (TBD) the IRS/IUE/STIS spectra.  `dust_extinction` provides the
fitting curves shapes (FM90, AxAvToExv).  `measure_extinction` gives the routines
to read the data and calculate extinction curves.

Curves calculated with `run_all_ext`.

(UV/Opt/)NIR/MIR portion fit with P92 shape using `fit_plasymdrude_all_ext`.

UV portion fit with FM90 shape using `fit_fm90_all_ext`.

Average extinction curve made using
`Figs/plot_mir_mext.py --alav --ave data/all_ext_18feb20_diffuse.dat`.

Figures
-------

1. MIR spectra comparison stars: Figs/plot_mir_specra_standards.py

2. MIR spectra reddened stars: Figs/plot_mir_specra_reddened.py

3. MIR spectra windy stars: Figs/plot_mir_specra_reddened_windy.py

4. Example fit: python utils/fit_mir_ext_powerlaw.py fits/hd112272_hd204172_ext.fits --notitle --nsteps=10000 --burnfrac=0.4

5. MIR E(lambda-V) curves: Figs/plot_mir_mext.py data/all_ext_18feb20.dat

6. UV+MIR extinction: Figs/plot_uv_mir_mext.py data/all_ext_18feb20_pldrude.dat --models

7. Sample properties: Figs/plot_sampprob.py data/all_ext_18feb20_pldrude.dat

8. UV+MIR Average diffuse extinction curve: Figs/plot_ave_ext.py

9. MIR ext literature comparison: Figs/plot_ext_litcomp.py

10. MIR ext dust grain model comparision: Figs/plot_ext_modcomp.py

11. Silicate (10 um) versus various: Figs/plot_silicate.py data/all_ext_18feb20.dat

Tables
------

1. Spiter Data: idl/make_ir_phot_table.pro

2. Opt/NIR Data: idl/make_ir_ophot_table.pro

3. Sightline properties: by hand

4. MIR extinction parameters: by hand

5. MIR general extinction parameters: Tables/gen_mir_ext_params_table.py data/all_ext_18feb20_pldrude.dat

6. MIR silicate extinction parameters: Tables/gen_mir_ext_params_table.py data/all_ext_18feb20_pldrude.dat --sil

7. UV extinction parameters: Tables/gen_uv_ext_params_table.py data/all_ext_18feb20_uv.dat
