

Most of these programs were modified from the WPIT module, available at https://github.com/stourgai/WPIT

First, let's begin by explaing the script "example_Figures_8l_8o_70deg_psi.py".
The script "example_Figures_8l_8o_70deg_psi.py" was modified from the WPIT test example (available at https://github.com/stourgai/WPIT/blob/main/WPIT_tests/Latitudinal%20dependence%20of%20nonlinear%20interaction%20between%20electromagnetic%20ion%20cyclotron%20wave%20and%20terrestrial%20ring%20current%20ions.ipynb).
Here are the modifications we suggest:
(1) It should be noted that the perpendicular wave vector kx is not zero for oblique EMIC waves. Therefore, for each step of the Runge-Kutta algorithm, we should potentially remove the line “kx_run=0”, which appears after the line: “Bxw_run, Byw_run, Bzw_run, Exw_run, Eyw_run, Ezw_run=wave.wave_amplitudes_bell(mu_run,P_run,D_run,S_run,Byw0_s,theta0)”  
(2) In each step of the Runge-Kutta method, there is a line of code following “ #calculate the wpi parameters”.
We should replace “B0” with “B_run”. Please refer to Eq. (18) of Su et al. (2014), where β is defined as:
β=-(k⊥ p⊥)/(qB_D )
In this equation, B_D represents the total magnetic field strength that varies with latitude (as shown in Eq. (4) of Su et al., 2014). 
Therefore, B_D should change as a simulated proton moves along the magnetic field line. 
However, the variable “B0” in the WPIT program refers to the magnetic field strength at the proton’s initial location. 
Thus, we should replace “B0” with “B_run”. For parallel EMIC waves, where k⊥ is zero, β will always be zero regardless of the “B0” used. 
In other words, using “B0” does not affect the calculation in this case.
 However, when the wave normal angle of EMIC waves is non-zero, k⊥ will be greater than zero. Under this condition, using “B0” will generate an incorrect β value, which affects the calculation of particle trajectories.

(3) There may be another issue regarding the fourth step of the Runge-Kutta method. We believe that the term “0.5*h*” should be changed to “h*”.

(4) The nonlinear behavior is highly sensitive to even a small change in the initial latitude of a proton. In our study, each proton was launched from its magnetic mirror point. 
To determine this point, we used Eq. (33) from Summers et al. (2007) to obtain the latitude of the magnetic mirror point. 
For a proton with an initial equatorial pitch angle αeq0=50°, the latitude of the magnetic mirror point is calculated to be “20.1854488137971” degree.

In addition to the modifications mentioned above, other changes may be required:
(5) There may be an issue in Line 6 of the function "dpperdt.py" (https://github.com/stourgai/WPIT/blob/main/WPIT/WPI_mod/EMIC_ion_mod/dpperdt.py). 
The term "ppar_arg-pwL_arg" should be changed to "ppar_arg+pwL_arg". Please refer to Eq. (15) in the work of Su et al. (2014).

(6) There may be an issue in Line 11 of the function "dalphadt.py (https://github.com/stourgai/WPIT/blob/main/WPIT/WPI_mod/EMIC_ion_mod/dalphadt.py). 
We believe that the term "fac2b-fac2c" should be revised to "fac2b+fac2c".  Please refer to Eq. (28) in Tourgaidis & Sarris (2022) for further details.

After applying these modifications, we could try to run the script "example_Figures_8l_8o_70deg_psi.py" to reproduce the case of 70° wave normal angle (WNA). 


Second, we explain the scripts attached to this page.
(1) "plot_dispersion_resLat.py"  calculates the dispersion relationship of EMIC waves
(2) "plot_s_affected_by_Bwy.py" calculates the dimensionless parameter |S| affected by the EMIC wave amplitude.
(3) "plot_trajectory_EMIC_Bwy.py" calculates the trajectories of three selected protons, corresponding to EMIC wave amplitudes of 0.1 nT, 1.0 nT, and 3.0 nT.
(4) "plot_s_affected_by_freq.py" computes the parameter |S| for 0.7Ω_He_eq, 0.8Ω_He_eq, and 0.9Ω_He_eq.
(5) "plot_trajectory_EMIC_freq.py" calculates the proton trajectories for 0.7Ω_He_eq, 0.8Ω_He_eq, and 0.9Ω_He_eq.
(6) "plot_s_affected_by_wna.py" plots the parameter |S| for various WNA values of 10°, 40°, and 70°
(7) "plot_Bessel_function_wna_10deg_40deg_70deg.py" analyzes the discontinuities in the parameter |S|
(8) "plot_trajectory_EMIC_wna.py" calculates the proton trajectories for WNA values of 10°, 40°, and 70°
(9) "plot_diffusion_coef.py" plots pitch angle advection and diffusion coefficients obtained from  "cal_50keV_diffusion_coef_678.py".
(10) "plot_trajectory_EMIC_freq_aeq404142.py" plots the trajectories of three selected protons with αeq0 of 40°, 41°, and 42°. The EMIC amplitude is 3 nT, the frequency is 0.7Ω_He_eq, and the wave normal angle is 0°.

To successfully run scirpts (1)-(10), the scripts "refr_index_full_left.py" and "numlib.py" are also required.
We would like to acknowledge Tourgaidis & Sarris (2022) for sharing the WPIT module on GitHub website (https://github.com/stourgai/WPIT).
For any questions, please feel free to contact me via this email (zhousujob@gmail.com)

References:
Tourgaidis, S., & Sarris, T. (2022). Wave-particle interactions toolset: A python-based toolset to model wave-particle interactions in the magnetosphere. Frontiers in Astronomy and Space Sciences, 9, 1005598. https://doi.org/10.3389/fspas.2022.1005598
Su, Z., Zhu, H., Xiao, F., Zheng, H., Zhang, M., Liu, Y. C.-M., et al. (2014). Latitudinal dependence of nonlinear interaction between electromagnetic ion cyclotron wave and terrestrial ring current ions. Physics Of Plasmas, 21(5). https://doi.org/10.1063/1.4880036
Summers, D., Ni, B., & Meredith, N. P. (2007). Timescales for radiation belt electron acceleration and loss due to resonant wave-particle interactions: 2. Evaluation for VLF chorus, ELF hiss, and electromagnetic ion cyclotron waves. Journal of Geophysical Research: Space Physics, 112(A4). https://doi.org/10.1029/2006JA011993






