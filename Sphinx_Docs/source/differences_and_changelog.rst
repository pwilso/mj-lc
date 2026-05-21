Differences and Changelog
===================================

Since changes to the newer versions of Python resulted in changes to various packages being used between Python v3.6 and v3.13 we provide two versions of the code, an older version and a newer version. It is possible that one version may receive updates that the other does not or that they behave slightly differently from one another. We give here a list of differences that between the two versions of the project.

* Python v3.13 no longer uses the *np.trapz()* function. The program now calls *np.trapezoid()*.
* The *blist* package could not be installed through pip in v3.13. *blist.sortedset()* was replaced with the standard *set()* and sorted.
* *np.complex()* was changed to *np.complex128()*
* *scipy.cos()* (and sin, arcsin, arcos) were replaced with *np.cos()*
* *scipy.seterr(invalid = 'ignore')* was replaced with *scipy.special.seterr(all = 'ignore')*
