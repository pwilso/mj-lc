Calculating the Coupling Factor for Semiconducting Devices
===================================

The physical processes being modelled can be split into calculations involved photons (optical) and charge carriers (electrical). Radiative recombination in a semiconducting device occurs when and electron/hole pair recombine to produce a photon. The rate at which this occurs takes the following form:

R_rad(z) = B_rad*n(z)*p(z)

where B_rad is the radiative recombination coefficient and is material dependent. The carrier concentrations for electrons and holes are given by n and p respectively. Here we use z to describe the depth within the device. The device is considered to have no lateral variations. The device can be considered one dimensional.

Once a photon is emitted at a given depth, z, it can travel through the device. It's path is governed by the bandgaps of the epi stack. It can escape the device or be reabsorbed within the device. Multiple reflections may play a significant role. The transfer matrix method is used to calculate how many photons are reabsorbed at a given point in the device for an emission event.
