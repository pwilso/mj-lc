An Introduction to Luminescent Coupling and Photon Recycling
===================================

In high-quality semiconducting materials, radiative recombination becomes the dominant recombination mechanism for electron-hole pairs within the device. These recombination events result in the production of a photon which can be absorbed again by the semiconducting material. This is particularly important in direct bandgap materials such as III-Vs, where the absorption coefficient is large. The emission and absorption events are position dependent within the device. A pair of emission and absorption events may occur at the same point within the device or they may occur at different locations within the device. The position dependence between the two events impacts the carrier statistics throughout the device and can lead to changes in measurable properties such as the short-circuit current. open-circuit voltage, and fill factor of the device.

In multi-junction devices, the emission-absorption process can be divided into two categories depending on the relative locations of the two events.

**Photon recycling**: The emission and absorption events occur within the same subcell of the multijunction device. This can result in an improvement to the V\ :sub:`oc`\  of the device.

**Luminescent coupling**: The emission and absorption events occur in different subcells of the multijunction device. This can result in a rebalancing of the photocurrent between the subcells if a current mismatch is present. If the mismatch is significant, luminescent coupling can improve the J\ :sub:`sc`\  of the device.

It is important to note that the physical process does not make the distinction between where one subcell ends and the next begins. These are simply definitions useful for discussing the physical process that takes place. The matrix produced by this code includes both types of processes implicitly as the coupling between every combination of emission and absorption locations are calculated and recorded.
