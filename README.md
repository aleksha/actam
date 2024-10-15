
# ACTAM

Monte-Carlo tool for the ACTAM Time Projection Chamber (TPC).

TPC will be used for studies of the cluster emissin form nuclei and for other other purposes.

**Note:** axial symmetry is assumed, therefor scattering occurs in YZ-plane. 

## To-do list

 1. Anode structure optimization for recoil ID
 2. Electronic noise
 3. Beam noise
 4. More realistic signal shape

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

