all: ../OpticalCoupling_Recombination.so.amd64

INCLUDEDIR=${STROOT}/tcad/current/lib/sdevice/include/

../OpticalCoupling_Recombination.so.amd64: OpticalCoupling_Recombination.cc
	${STROOT}/bin/cmi -v -O -I${INCLUDEDIR} -lm OpticalCoupling_Recombination.cc
	mv OpticalCoupling_Recombination.so.* ../


clean:
	rm ../OpticalCoupling_Recombination
