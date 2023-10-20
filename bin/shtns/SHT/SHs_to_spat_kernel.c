// This file was automatically generated by 'make' from file 'kernel_SH_to_spat.gen.c'.
// To modify it, please consider modifying kernel_SH_to_spat.gen.c
/*
 * Copyright (c) 2010-2020 Centre National de la Recherche Scientifique.
 * written by Nathanael Schaeffer (CNRS, ISTerre, Grenoble, France).
 * 
 * nathanael.schaeffer@univ-grenoble-alpes.fr
 * 
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * 
 */

//////////////////////////////////////////////////

	#ifdef HI_LLIM
	#ifndef SHT_GRAD
	#define BASE _sy2_hi
	#else
	#define BASE _sy1s_hi
	#endif
	#else
	#ifndef SHT_GRAD
	#define BASE _sy2
	#else
	#define BASE _sy1s
	#endif
	#endif

  #ifndef SHT_GRAD
	void GEN3(BASE,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1) {
  #else
	void GEN3(BASE,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1) {
  #endif

	#if !defined( _GCC_VEC_ ) && (NWAY & 1)
	#error "NWAY must be even when compiled without explicit vectorization."
	#endif
	#if VSIZE2*NWAY > 32
	#error "VSIZE2*NWAY must not exceed 32"
	#endif

  #ifndef SHT_AXISYM
   #ifndef SHTNS_ISHIOKA
   #else
   #endif
	#define vr(l) vall( ((double*) VWl)[4*(l)]   )
	#define vi(l) vall( ((double*) VWl)[4*(l)+1] )
	#define wr(l) vall( ((double*) VWl)[4*(l)+2] )
	#define wi(l) vall( ((double*) VWl)[4*(l)+3] )
  #endif
	long int nk,k,l,m;
	double *alm, *al;
	double *ct, *st;
	int robert_form;
	v2d VWl[llim*2+4];
  #ifdef SHTNS_ISHIOKA
  #endif

	ct = shtns->ct;		st = shtns->st;
	nk = it1;	//NLAT_2;
	#if _GCC_VEC_
		nk = ((unsigned)(nk+VSIZE2-1)) / VSIZE2;
		it0 = ((unsigned)(it0+VSIZE2-1)) / VSIZE2;
	#endif
	robert_form = shtns->robert_form;

	if (im == 0)
	{	//	im=0;
		double* const Sl0 = (double*) VWl;
		#ifdef SHT_GRAD
			// TODO: fix k,nk bounds
			if (BpF != NULL) memset(BpF, 0, sizeof(v2d) * NLAT_2);
		#endif
 		l=1;
		alm = shtns->alm;
		do {		// for m=0, compress the complex Q,S,T to double
			Sl0[l-1] = creal( Slm[l] );	//	Sl[l] = (double) Slm[l+1];
			++l;
		} while(l<=llim);
		k=it0;
		do {
			l=0;	al = alm;
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
			rnd sint[NWAY], dy0[NWAY], dy1[NWAY];
			rnd te[NWAY], to[NWAY];
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(ct, j+k);
				sint[j] = -vread(st, j+k);
				y0[j] = vall(al[0]);
				dy0[j] = vall(0.0);
				to[j] = dy0[j];
			}
			if (robert_form) {
				for (int j=0; j<NWAY; ++j) sint[j] *= -sint[j];
			}
			for (int j=0; j<NWAY; ++j) {
				y1[j]  = vall(al[0]*al[1]) * cost[j];
				dy1[j] = vall(al[0]*al[1]) * sint[j];
			}
			for (int j=0; j<NWAY; ++j) {
				te[j] = dy1[j] * vall(Sl0[0]);
			}
			al+=2;	l+=2;
			while(l<llim) {
				for (int j=0; j<NWAY; ++j) {
					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*sint[j]) + vall(al[0])*dy0[j];
					y0[j]  = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
					to[j] += dy0[j] * vall(Sl0[l-1]);
				}
				for (int j=0; j<NWAY; ++j) {
					dy1[j] = vall(al[3])*(cost[j]*dy0[j] + y0[j]*sint[j]) + vall(al[2])*dy1[j];
					y1[j]  = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				for (int j=0; j<NWAY; ++j) {
					te[j] += dy1[j] * vall(Sl0[l]);
				}
				al+=4;	l+=2;
			}
			if (l==llim) {
				for (int j=0; j<NWAY; ++j) {
					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*sint[j]) + vall(al[0])*dy0[j];
					y0[j]  = vall(al[1])*cost[j]*y1[j] + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
					to[j] += dy0[j] * vall(Sl0[l-1]);
				}
			}
			// combine even/odd into north/south
			for (int j=0; j<NWAY; ++j) {			
				rnd ts = te[j] - to[j];	te[j] = te[j] + to[j];
				to[j] = ts;
			}
		#ifndef SHTNS4MAGIC
			for (int j=0; j<NWAY; ++j) {
				S2D_STORE(BtF, j+k, te[j], to[j])
			}
		#else
			for (int j=0; j<NWAY; ++j) {
				if ((k+j)>=nk) break;
				S2D_STORE_4MAGIC((double*)BtF, j+k, te[j], to[j]);
			}
		#endif
			k+=NWAY;
		} while (k < nk);
	}
  #ifndef SHT_AXISYM
	else
	{		// im > 0
		BtF += im*(shtns->nlat_padded >>1);
		BpF += im*(shtns->nlat_padded >>1);
		m = im*MRES;
		l = (im*(2*(LMAX+1)-(m+MRES)))>>1;		//l = LiM(shtns, 0,im);
		#ifndef SHTNS_ISHIOKA
		alm = shtns->alm + 2*(l+m);		// shtns->alm + im*(2*(LMAX+1) -m+MRES);
		#else
		alm = shtns->clm + (l+m);		// shtns->clm + im*(2*(LMAX+1) -m+MRES)/2;
		#endif

  #ifndef SHT_GRAD
		SH_vect_to_2scal(shtns->mx_stdt + 2*l, llim, m, &Slm[l], &Tlm[l], (cplx*) VWl);
  #else
		SHsph_to_2scal(shtns->mx_stdt + 2*l, llim, m, &Slm[l], (cplx*) VWl);
  #endif

	#ifndef SHTNS_ISHIOKA
	#else
		// pre-processing for recurrence relation of Ishioka
		const double* restrict xlm = shtns->xlm + 3*im*(2*(LMAX+4) -m+MRES)/4;
		SH2_to_ishioka(xlm, VWl+2*m, llim-m+1);
	#endif

		// polar optimization
		k = shtns->tm[im];		// start index in theta (=0 if no polar optimization)
		#if _GCC_VEC_
		k = ((unsigned) k) / VSIZE2;	// in vector size units.
		#else
		k = (k>>1)*2;		// k must be even.
		#endif
		if (it0 < k) {
			const long ofsm = (NPHI-2*im)*(shtns->nlat_padded >>1);
		#ifndef SHTNS4MAGIC
			#if _GCC_VEC_
			const long ofs1 = NLAT_2 - k*(VSIZE2/2);
			#else
			const long ofs1 = NLAT_2 - k/2;
			#endif
			zero_poles4_vect(BtF+it0*(VSIZE2/2), ofsm, ofs1, k-it0);
			zero_poles4_vect(BpF+it0*(VSIZE2/2), ofsm, ofs1, k-it0);
		#else
			zero_poles2_vect(BtF+it0*VSIZE2, ofsm, 2*(k-it0));
			zero_poles2_vect(BpF+it0*VSIZE2, ofsm, 2*(k-it0));
		#endif
		} else k=it0;

		do {
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
			rnd ter[NWAY], tei[NWAY], tor[NWAY], toi[NWAY];
			rnd per[NWAY], pei[NWAY], por[NWAY], poi[NWAY];
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(st, k+j);
				y0[j] = vall(1.0);
			}
			l=m;
			if (robert_form == 0) l=m-1;
	#ifndef HI_LLIM
			const long int ny = 0;
			while(1) {		// sin(theta)^m
				if (l&1) for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
				l >>= 1;
				if (l==0) break;
				for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
			}
	#else
			long int ny = 0;
			for (long int nsint = 0;;) {		// sin(theta)^m		(use rescaling to avoid underflow)
				if (l&1) {
					for (int j=NWAY-1; j>=0; --j) y0[j] *= cost[j];
					ny += nsint;
					if (vlo(y0[NWAY-1]) < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {
						ny--;
						for (int j=NWAY-1; j>=0; --j) y0[j] *= vall(SHT_SCALE_FACTOR);
					}
				}
				l >>= 1;
				if (l==0) break;
				for (int j=NWAY-1; j>=0; --j) cost[j] *= cost[j];
				nsint += nsint;
				if (vlo(cost[NWAY-1]) < 1.0/SHT_SCALE_FACTOR) {
					nsint--;
					for (int j=NWAY-1; j>=0; --j) cost[j] *= vall(SHT_SCALE_FACTOR);
				}
			}
	#endif
			al = alm;
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(ct, j+k);
				#ifndef SHTNS_ISHIOKA
				y0[j] *= vall(al[0]);
				#else
				cost[j] *= cost[j];		// cos(theta)^2
				#endif
			}
			for (int j=0; j<NWAY; ++j) {
				#ifndef SHTNS_ISHIOKA
				y1[j]  = (vall(al[1])*y0[j]) *cost[j];		//	y1[j] = vall(al[1])*cost[j]*y0[j];
				#else
				y1[j] = (vall(al[1])*cost[j] + vall(al[0]))*y0[j];
				#endif
				por[j] = vall(0.0);		tei[j] = vall(0.0);
				tor[j] = vall(0.0);		pei[j] = vall(0.0);
				poi[j] = vall(0.0);		ter[j] = vall(0.0);
				toi[j] = vall(0.0);		per[j] = vall(0.0);
			}
			l=m;		al+=2;
	#ifdef HI_LLIM
		  if (ny < 0) {		// ylm treated as zero and ignored if ny < 0
			const rnd scale = vall(1.0/SHT_SCALE_FACTOR);
			while (l<llim) {
				#ifndef SHTNS_ISHIOKA
				for (int j=0; j<NWAY; ++j)	y0[j] = (vall(al[1])*cost[j])*y1[j] + vall(al[0])*y0[j];
				for (int j=0; j<NWAY; ++j)	y1[j] = (vall(al[3])*cost[j])*y0[j] + vall(al[2])*y1[j];
				al+=4;
				#else
				rnd a[NWAY];
				for (int j=0; j<NWAY; ++j)	a[j] = vall(al[1])*cost[j] + vall(al[0]);
				al+=2;
				for (int j=0; j<NWAY; ++j) {
					a[j] = a[j]*y1[j] + y0[j];
					y0[j] = y1[j];		y1[j] = a[j];
				}
				#endif
				l+=2;
				if (fabs(vlo(y0[NWAY-1])) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0) {		// rescale when value is significant
					for (int j=0; j<NWAY; ++j) {
						y0[j] *= scale;		y1[j] *= scale;
					}
					if (++ny == 0) break;
				}
			}
		  }
	#endif
		  if LIKELY(ny == 0) {
		#ifndef SHTNS_ISHIOKA
			while (l<llim) {	// compute even and odd parts
				for (int j=0; j<NWAY; ++j) {	ter[j] += y0[j]  * vr(l);		tei[j] += y0[j] * vi(l);	}
				for (int j=0; j<NWAY; ++j) {	per[j] += y0[j]  * wr(l);		pei[j] += y0[j] * wi(l);	}
				for (int j=0; j<NWAY; ++j) {
					y0[j] = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {	tor[j] += y1[j]  * vr(l+1);		toi[j] += y1[j] * vi(l+1);	}
				for (int j=0; j<NWAY; ++j) {	por[j] += y1[j]  * wr(l+1);		poi[j] += y1[j] * wi(l+1);	}
				for (int j=0; j<NWAY; ++j) {
					y1[j] = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				l+=2;	al+=4;
			}
				for (int j=0; j<NWAY; ++j) {	ter[j] += y0[j]  * vr(l);		tei[j] += y0[j] * vi(l);	}
				for (int j=0; j<NWAY; ++j) {	per[j] += y0[j]  * wr(l);		pei[j] += y0[j] * wi(l);	}
			if (l==llim) {
				for (int j=0; j<NWAY; ++j) {	tor[j] += y1[j]  * vr(l+1);		toi[j] += y1[j] * vi(l+1);	}
				for (int j=0; j<NWAY; ++j) {	por[j] += y1[j]  * wr(l+1);		poi[j] += y1[j] * wi(l+1);	}
			}
			// combine even/odd into north/south
			for (int j=0; j<NWAY; ++j) {
				rnd sr = ter[j] - tor[j];	ter[j] = ter[j] + tor[j];
			  	rnd si = tei[j] - toi[j];	tei[j] = tei[j] + toi[j];
				tor[j] = sr;		toi[j] = si;
				sr = per[j] - por[j];	per[j] = per[j] + por[j];
			  	si = pei[j] - poi[j];	pei[j] = pei[j] + poi[j];
				por[j] = sr;		poi[j] = si;
			}
		#else
			while (l<llim) {	// compute even and odd parts
				for (int j=0; j<NWAY; ++j)	ter[j] += y0[j] * vr(l);
				for (int j=0; j<NWAY; ++j)	tei[j] += y0[j] * vi(l);
				for (int j=0; j<NWAY; ++j)	per[j] += y0[j] * wr(l);
				for (int j=0; j<NWAY; ++j)	pei[j] += y0[j] * wi(l);
				for (int j=0; j<NWAY; ++j)	tor[j] += y0[j] * vr(l+1);
				for (int j=0; j<NWAY; ++j)	toi[j] += y0[j] * vi(l+1);
				for (int j=0; j<NWAY; ++j)	por[j] += y0[j] * wr(l+1);
				for (int j=0; j<NWAY; ++j)	poi[j] += y0[j] * wi(l+1);
				for (int j=0; j<NWAY; ++j) {
					rnd tmp = (vall(al[1])*cost[j] + vall(al[0]))*y1[j] + y0[j];
					y0[j] = y1[j];
					y1[j] = tmp;
				}
				l+=2;	al+=2;
			}
			for (int j=0; j<NWAY; ++j) cost[j] = vread(ct, k+j);		// read ahead to correct the odd part below
			for (int j=0; j<NWAY; ++j)	ter[j] += y0[j] * vr(l);
			for (int j=0; j<NWAY; ++j)	tei[j] += y0[j] * vi(l);
			for (int j=0; j<NWAY; ++j)	per[j] += y0[j] * wr(l);
			for (int j=0; j<NWAY; ++j)	pei[j] += y0[j] * wi(l);
			if LIKELY(l==llim) {
				for (int j=0; j<NWAY; ++j)	tor[j] += y0[j] * vr(l+1);
				for (int j=0; j<NWAY; ++j)	toi[j] += y0[j] * vi(l+1);
				for (int j=0; j<NWAY; ++j)	por[j] += y0[j] * wr(l+1);
				for (int j=0; j<NWAY; ++j)	poi[j] += y0[j] * wi(l+1);
			}
			// correct the odd part:
		//	for (int j=0; j<NWAY; ++j) {  tor[j] *= cost[j];	toi[j] *= cost[j]; }
		//	for (int j=0; j<NWAY; ++j) {  por[j] *= cost[j];	poi[j] *= cost[j]; }

			// combine even/odd into north/south, and correct the odd part for free with FMA
			for (int j=0; j<NWAY; ++j) {
				rnd sr = ter[j] - tor[j]*cost[j];	ter[j] = ter[j] + tor[j]*cost[j];
			  	rnd si = tei[j] - toi[j]*cost[j];	tei[j] = tei[j] + toi[j]*cost[j];
				tor[j] = sr;		toi[j] = si;
				sr = per[j] - por[j]*cost[j];	per[j] = per[j] + por[j]*cost[j];
			  	si = pei[j] - poi[j]*cost[j];	pei[j] = pei[j] + poi[j]*cost[j];
				por[j] = sr;		poi[j] = si;
			}
		#endif
		  }
		#ifndef SHTNS4MAGIC
		  #ifdef _GCC_VEC_
			const long ofs = (NPHI-2*im)*shtns->nlat_padded;
			for (int j=0; j<NWAY; ++j) cstore_north_south((double*) BtF, ((double*) (BtF)) +ofs, k+j, NLAT, ter[j], tor[j], tei[j], toi[j]);
			for (int j=0; j<NWAY; ++j) cstore_north_south((double*) BpF, ((double*) (BpF)) +ofs, k+j, NLAT, per[j], por[j], pei[j], poi[j]);
		  #else
		  	// NWAY is even when _GCC_VEC_ is not defined
			for (int j=0; j<NWAY/2; ++j) {	S2D_CSTOREX(BtF, k/2+j, 2*j, ter, tor, tei, toi)  }
			for (int j=0; j<NWAY/2; ++j) {	S2D_CSTOREX(BpF, k/2+j, 2*j, per, por, pei, poi)  }
		  #endif
		#else
			for (int j=0; j<NWAY; ++j) {
				if ((k+j)>=nk) break;
				S2D_CSTORE_4MAGIC((double*) BtF, (double*) (BtF + (NPHI-2*im)*(shtns->nlat_padded>>1)), k+j, ter[j], tor[j], tei[j], toi[j]);
				S2D_CSTORE_4MAGIC((double*) BpF, (double*) (BpF + (NPHI-2*im)*(shtns->nlat_padded>>1)), k+j, per[j], por[j], pei[j], poi[j]);
			}
		#endif
			k+=NWAY;
		} while (k < nk);
	}
  #endif
}

	#undef sr
	#undef si

	#undef BASE
