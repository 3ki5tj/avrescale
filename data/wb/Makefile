t0 = 1000000

etraj::
	cd fix     && ../mkhist.py --col=3   --t0=$(t0) --dx=0.0001 -o null && cd ..
	cd reg     && ../mkhist.py --col=3   --t0=$(t0) --dx=0.0001 -o null && cd ..
	cd reg100  && ../mkhist.py --col=3   --t0=$(t0) --dx=0.0001 -o null && cd ..
	cd adp     && ../mkhist.py --col=3   --t0=$(t0) --dx=0.0001 -o null && cd ..

uhist::
	cd adp     && ../mkhist.py --col=2   --t0=$(t0) --dx=2 -o Uadp.his     && cd ..
	cd fix     && ../mkhist.py --col=2   --t0=$(t0) --dx=2 -o Ufix.his     && cd ..
	cd can_vrs && ../mkhist.py --col=2   --t0=$(t0) --dx=2 -o Ucan_vrs.his && cd ..
	cd reg     && ../mkhist.py --col=2   --t0=$(t0) --dx=2 -o Ureg.his     && cd ..
	cd reg100  && ../mkhist.py --col=2   --t0=$(t0) --dx=2 -o Ureg.his     && cd ..

khist::
	cd adp     && ../mkhist.py --col=dif --t0=$(t0) --dx=2 -o Kadp.his && cd ..
	cd fix     && ../mkhist.py --col=dif --t0=$(t0) --dx=2 -o Kfix.his && cd ..
	cd can_vrs && ../mkhist.py --col=dif --t0=$(t0) --dx=2 -o Kcan_vrs.his && cd ..
	cd reg     && ../mkhist.py --col=dif --t0=$(t0) --dx=2 -o Kreg.his     && cd ..
	cd reg100  && ../mkhist.py --col=dif --t0=$(t0) --dx=2 -o Kreg.his     && cd ..

wberrz::
	gcc -O3 prederrz.c -o prederrz -lm
	./prederrz 100000  0.01 > wberrz_t100000.dat
	./prederrz 1000000 0.01 > wberrz_t1M.dat

