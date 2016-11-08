etraj:
	cd fix && ../mkhist.py --col=3 --t0=1000000 --dx=0.0001 -o null && cd ..
	cd reg && ../mkhist.py --col=3 --t0=1000000 --dx=0.0001 -o null && cd ..
	cd adp && ../mkhist.py --col=3 --t0=1000000 --dx=0.0001 -o null && cd ..


uhist:
	cd adp && ../mkhist.py --col=dif --t0=1000000 --dx=2 -o Kadp.his && cd ..
	cd fix && ../mkhist.py --col=dif --t0=1000000 --dx=2 -o Kfix.his && cd ..
	cd can && ../mkhist.py --col=dif --t0=1000000 --dx=2 -o Kcan.his && cd ..
	cd can_vrs && ../mkhist.py --col=dif --t0=1000000 --dx=2 -o Kcan_vrs.his && cd ..
	cd can2 && ../mkhist.py --col=dif --t0=1000000 --dx=2 -o Kcan.his && cd ..
