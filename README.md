
# ACTAM

Monte-Carlo tool for the ACTAM Time Projection Chamber (TPC).

TPC will be used for studies of the cluster emissin form nuclei and for other other purposes.

**Note:** axial symmetry is assumed, therefor scattering occurs in YZ-plane. 

## To-do list

 1. Anode structure optimization for recoil ID
 2. Electronic noise [done]
 3. Beam noise [done]
 4. More realistic signal shape [done]

## Parameters in CONFIG.txt

 - **Line1:** Pressure of Argon (`float` value)
 - **Line2:** Code of recoil (beam) particle (`int` value)
 - **Line3:** Limits for kinetic energy of the recoil (two `float` values)
 - **Line4:** Limits of angle wrt z-axis **in degrees** for the recoil (two `float` values)
 - **Line5:** XYZ-position of the recoil wrt center of the TPC volume (three `float` values)

### Code for recoil (beam) particles

 * 0 - proton (H-1)
 * 1 - deuteron (H-2)
 * 2 - triton (H-3)
 * 3 - He-3
 * 4 - alpha (He-4)
 * 11 - electron
 * 22 - gamma 

## How-to simulate N events

 1. Simulate electron beam events (factor 1000 to N) with the code. Rename `out.data` to `out.data.beam`
 2. Simulate noise with `anode_noise_gen.py` from  https://github.com/aleksha/electronic-noise (factor 3 to N)
 3. Simalate recoil with same gas pressure as at the step 1.
 4. Build events (in progerss, currently `analyze2.C`).
 
