# NLstep
NLstep is an academic Matlab tool for time step integration of mechanical systems with the special ability to treat sharply or even non-regularized frictional-unilateral contact. Typically, an FE model is used as point of departure, and component mode synthesis is then applied, not only to reduce the model order, but also to obtain a massless boundary, which is useful for mitigating numerical oscillations [1]. The tool permits to simulate the  system's response starting from a given set of initial conditions under possible forcing, imposed contact gap and velocity at the contacts, each of which may have an arbitrary explicit time dependence. Initially closed, preloaded contacts can also be handled. 

A central motivation behind NLstep is to overcome the limitations of Harmonic Balance (as implemented in NLvib [2]). In particular, NLstep permits to study transient dynamics as well as non-periodi steady state dynamics (e.g. deterministic chaos). NLstep is competitive alternative to Harmonic Balance also in the periodic case, especially for sharply or even non-regularized contact, where the transitions (stick/slip; opening/closing) slow down the  harmonic convergence. 

In the following, an overview is given over the tool's key functionality, included examles, and requirements. 

# Key functionality 
Two functions are contained in `SRC` folder: 

`sim_contact3D.m`: This function treats 3D frictional-unilateral contact in a non-regularized way (set-valued Coulomb-Signorini conditions). The model of the mechanical system must have a massless boundary, which is ensured by requiring that it has been obtained by an appopriate component mode synthesis method implemented in the class 'CMS_ROM'. A leapfrog integration scheme is used, where the dynamic internal force balance is stepped in an explicit way, while the static force balance at the boundary is treated in an implicit way. An event caturing strategy is pursued where the Coulomb-Signorini conditions are expressed via a projective non-smooth equation system, which is iteratively solved using projected Jacobi over-relaxation. For more information, please see the header of that function.

`sim_regularNL.m`: This function is designed regular nonlinear behavior. It is slightly more general than the above function in the sense that it is not restricted to contact, but allows to treat quite arbitrary local nonlinear elements, and does not require a massless-boundary model. The local nonlinear elements are defined as in NLvib [2], but in contrast to the public NLvib version, 3D penalty contact is implemented in this function. The time stepping relies on the Newmark's quadrature rules in conjunction with equilibrium averaging. Popular particular forms include the HHT and the Generalized-alpha method. The algebraic equation system governing the coordinate vector at the next time level is solved using the Newton method with a Cholesky factorization of the analytical Jacobian. 

# Included examples 
The included examples illustrate the usage of NLstep. Typically, this entails the sub-tasks: 

        1. Set up an initial FE model.
        
        2. Apply component mode synthesis to obtain a reduced-order model.
        
        3. Define damping, forcing and nonlinear elements.
        
        4. Call the integrator.
        
        5. Carry out basic post-processing.

The following content is found in the `EXAMPLES` folder:

`01_bouncingBar`: common bouncing bar problem; similar results are shown in [1]; features a 1D unilateral contact variant of `sim_contact3D.m`
        
`02_beamFrictionFrequencySweep`: cantilvered Euler-Bernoulli beam subjected to dry friction harmonically driven near resonance with the fundamental bending mode; comparison against Harmonic Balance results; features a 1D frictional contact variant of `sim_contact3D.m`
       
`03_bladeCasingRubbing`: 3D FE model of a rotating compressor blade subjected to frictional impacts with a rigid oval casing; similar results are shown in [1]
       
`04_obliqueContactStickSlipLiftoff`: 3D FE model of a cantilever subjected to oblique contact at its free end; the contact is initially preloaded but undergoes stick-slip-liftoff transitions under dynamic loading; illustrates how node-based quadrature of the contact stress can be implemented; illustrates consistency and benefit of eliminating presumed sticking contacts from the projective equation system in `sim_contact3D.m`
       
`05_squeakNoiseTestRig`: cantilevered beam with a curved plate bolted to its free end, pressed against a flat plate, mimicking squeak noise generation [3]; illustrates simulation with penalty regularized contact

# Requirements
The tool relies on NLvib [2], primarily for the definition of mechanical systems including nonlinear (contact) elements. In fact, the two classes, `FEmodel` and `CMS_ROM`, have been added to NLvib which were only used within NLstep by the time of the initial release (yet these seemed useful for NLvib also). These import / handle data exported from a conventional FE tool, and apply component mode synthesis as required in NLstep, respectively. For more information, please see the header of the respective class. 

The tool should work well with a wide range of Matlab releases, without the need for any toolbox. There is currently one exception: If you want to run NLvib's solve_and_continue.m, as it is done in the example `02_beamFrictionFrequencySweep` to compute the Harmonic Balance approximation, you need the optimization toolbox.

NLstep was developed by Malte Krack and Johann Gross. If you use NLstep, please refer to our article [1]. We always appreciate any kind of feedback you may have. If you encounter any problems, which you cannot solve, please do not hesitate to contact the authors of this code (malte.krack@ila.uni-stuttgart.de; johann.gross@ila.uni-stuttgart.de).

# References 

[1] Monjaraz-Tec, C.; Gross, J.; Krack, M.: "A massless boundary component mode synthesis method for elastodynamic contact problems." Computers & Structures 260 (2022): 106698. https://doi.org/10.1016/j.compstruc.2021.106698 

[2] https://www.ila.uni-stuttgart.de/nlvib/; https://github.com/maltekrack/NLvib/ 

[3] Utzig, L.; Weisheit, K.; Sepahvand, K.; Marburg, S.: "Innovative  squeak noise prediction: An approach using the harmonic balance method and a variable normal contact force".Journal of Sound and %   Vibration 501 (2021): 116077. https://doi.org/10.1016/j.jsv.2021.116077 
