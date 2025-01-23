
These scripts are designed to study the nonlinear effects of ring current protons interacting with electromagnetic cyclotron (EMIC) waves.

Most of these programs were modified from the WPIT module, available at https://github.com/stourgai/WPIT

Script "example_Figures_8l_8o_70deg_psi.py" calculates the trajectories of 50 keV protons interacting with the oblique EMIC waves (wave normal angle is 70°).
This script was modified from the WPIT test example (available at https://github.com/stourgai/WPIT/blob/main/WPIT_tests/Latitudinal%20dependence%20of%20nonlinear%20interaction%20between%20electromagnetic%20ion%20cyclotron%20wave%20and%20terrestrial%20ring%20current%20ions.ipynb).

For detailed information regarding the modifications we made, please contact Dr. Su Zhou (zhousujob@gmail.com)

Below, we provide an explanation of the scripts attached to this page:

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







