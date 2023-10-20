// This file was automatically generated by 'make' from file 'fly_SH_to_spat.gen.c'.
// To modify it, please consider modifying fly_SH_to_spat.gen.c
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
 

	#ifndef SHT_AXISYM

  #ifndef SHT_GRAD
	void GEN3(_sy2,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
  #else
	void GEN3(_sy1s,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
	void GEN3(_sy1t,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
  #endif

  #ifndef SHT_GRAD
	void GEN3(_sy2_hi,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
  #else
	void GEN3(_sy1s_hi,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
	void GEN3(_sy1t_hi,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, const unsigned im, int it0, int it1);
  #endif
  
  #endif

  #ifndef SHT_GRAD
	void GEN3(_sy2,NWAY,_m0l)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, int it0, int it1);
  #else
	void GEN3(_sy1s,NWAY,_m0l)(shtns_cfg shtns, cplx *Slm, v2d *BtF, v2d *BpF, const long int llim, int it0, int it1);
	void GEN3(_sy1t,NWAY,_m0l)(shtns_cfg shtns, cplx *Tlm, v2d *BtF, v2d *BpF, const long int llim, int it0, int it1);
  #endif


  #ifndef SHT_GRAD
	static void GEN3(SHsphtor_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, cplx *Tlm, double *Vt, double *Vp, const long int llim) {
  #else
	static void GEN3(SHsph_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Slm, double *Vt, double *Vp, const long int llim) {
	static void GEN3(SHtor_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, cplx *Tlm, double *Vt, double *Vp, const long int llim) {
  #endif

	unsigned imlim = 0;
	v2d* BtF = (v2d*) Vt;	v2d* BpF = (v2d*) Vp;

  #ifndef SHT_AXISYM
	imlim = MTR;
	#ifdef SHT_VAR_LTR
		if (imlim*MRES > (unsigned) llim) imlim = ((unsigned) llim)/MRES;		// 32bit mul and div should be faster
	#endif
	if (shtns->fftc_mode > 0) {		// alloc memory for the FFT
		unsigned long nv = shtns->nspat;
		BtF = (v2d*) VMALLOC( 2*nv * sizeof(double) );
		BpF = BtF + nv/2;
	}
  #endif

	const int it0 = 0;
	const int it1 = NLAT_2;
	#ifndef SHT_AXISYM
	if (llim < SHT_L_RESCALE_FLY) {
	#endif
		for (int im=0; im <= imlim; im++)
		{
	#ifndef SHT_GRAD
			GEN3(_sy2,NWAY,_l)(shtns, Slm, Tlm, BtF, BpF, llim, im, it0, it1);
	#else
			GEN3(_sy1s,NWAY,_l)(shtns, Slm, BtF, BpF, llim, im, it0, it1);
			GEN3(_sy1t,NWAY,_l)(shtns, Tlm, BtF, BpF, llim, im, it0, it1);
	#endif
		}
  #ifndef SHT_AXISYM
	} else {

		for (int im=0; im <= imlim; im++)
		{
	#ifndef SHT_GRAD
				GEN3(_sy2_hi,NWAY,_l)(shtns, Slm, Tlm, BtF, BpF, llim, im, it0, it1);
	#else
				GEN3(_sy1s_hi,NWAY,_l)(shtns, Slm, BtF, BpF, llim, im, it0, it1);
				GEN3(_sy1t_hi,NWAY,_l)(shtns, Tlm, BtF, BpF, llim, im, it0, it1);
	#endif
		}
	}

	// padding for high m's
	if (NPHI-1 > 2*imlim) {
		const int m_inc = shtns->nlat_padded >> 1;
		memset(BtF + m_inc*(imlim+1), 0, sizeof(cplx)* m_inc * (NPHI-1-2*imlim));
		memset(BpF + m_inc*(imlim+1), 0, sizeof(cplx)* m_inc * (NPHI-1-2*imlim));
	}

    // NPHI > 1 as SHT_AXISYM is not defined.
  	if (shtns->fftc_mode >= 0) {
		if (shtns->fftc_mode != 1) {
			fftw_execute_dft(shtns->ifftc, (cplx *) BtF, (cplx *) Vt);
			fftw_execute_dft(shtns->ifftc, (cplx *) BpF, (cplx *) Vp);
		} else {		// split dft
			fftw_execute_split_dft(shtns->ifftc,((double*)BtF)+1, ((double*)BtF), Vt+NPHI, Vt);
			fftw_execute_split_dft(shtns->ifftc,((double*)BpF)+1, ((double*)BpF), Vp+NPHI, Vp);
			VFREE(BtF);		// this frees also BpF.
		}
	}
  #endif

  }




  #ifndef SHT_AXISYM

  #ifndef SHT_GRAD
	static void GEN3(SHsphtor_m_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, int im, cplx *Slm, cplx *Tlm, cplx *Vt, cplx *Vp, const long int llim) {
  #else
	static void GEN3(SHsph_m_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, int im, cplx *Slm, cplx *Vt, cplx *Vp, const long int llim) {
	static void GEN3(SHtor_m_to_spat_fly,NWAY,SUFFIX)(shtns_cfg shtns, int im, cplx *Tlm, cplx *Vt, cplx *Vp, const long int llim) {
  #endif

	v2d *BtF, *BpF;
	#define vr(l) vall( ((double*) VWl)[4*(l)]   )
	#define vi(l) vall( ((double*) VWl)[4*(l)+1] )
	#define wr(l) vall( ((double*) VWl)[4*(l)+2] )
	#define wi(l) vall( ((double*) VWl)[4*(l)+3] )
	long int nk, k,l,m;
	double *alm, *al;
	s2d *ct, *st;
	int robert_form;

	BtF = (v2d*) Vt;	BpF = (v2d*) Vp;

	nk = NLAT_2;
	#if _GCC_VEC_
		nk = ((unsigned)(nk+VSIZE2-1)) / VSIZE2;
	#endif
	ct = (s2d*) shtns->ct;		st = (s2d*) shtns->st;
	robert_form = shtns->robert_form;

	if (im == 0) {
		double Sl0[llim];
		double Tl0[llim];

		#ifdef SHT_GRAD
			k=0; do { BpF[k]=vdup(0.0); } while(++k<NLAT);
			k=0; do { BtF[k]=vdup(0.0); } while(++k<NLAT);
		#endif

 		l=1;
		alm = shtns->alm;
		do {		// for m=0, compress the complex Q,S,T to double
			Sl0[l-1] = (double) Slm[l];	//	Sl[l] = (double) Slm[l+1];
			Tl0[l-1] = (double) Tlm[l];	//	Tl[l] = (double) Tlm[l+1];
			++l;
		} while(l<=llim);
		k=0;
		do {
			l=0;	al = alm;
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
			rnd sint[NWAY], dy0[NWAY], dy1[NWAY];
			rnd te[NWAY], to[NWAY];
			rnd pe[NWAY], po[NWAY];
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(ct, j+k);
				sint[j] = -vread(st, j+k);
				y0[j] = vall(al[0]);
				dy0[j] = vall(0.0);
				to[j] = dy0[j];
				po[j] = dy0[j];
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
				pe[j] = -dy1[j] * vall(Tl0[0]);
			}
			al+=2;	l+=2;
			while(l<llim) {
				for (int j=0; j<NWAY; ++j) {
					dy0[j] = vall(al[1])*(cost[j]*dy1[j] + y1[j]*sint[j]) + vall(al[0])*dy0[j];
					y0[j]  = vall(al[1])*(cost[j]*y1[j]) + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
					to[j] += dy0[j] * vall(Sl0[l-1]);
					po[j] -= dy0[j] * vall(Tl0[l-1]);
				}
				for (int j=0; j<NWAY; ++j) {
					dy1[j] = vall(al[3])*(cost[j]*dy0[j] + y0[j]*sint[j]) + vall(al[2])*dy1[j];
					y1[j]  = vall(al[3])*(cost[j]*y0[j]) + vall(al[2])*y1[j];
				}
				for (int j=0; j<NWAY; ++j) {
					te[j] += dy1[j] * vall(Sl0[l]);
					pe[j] -= dy1[j] * vall(Tl0[l]);
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
					po[j] -= dy0[j] * vall(Tl0[l-1]);
				}
			}
			// combine even/odd into north/south
			for (int j=0; j<NWAY; ++j) {			
				rnd ts = te[j] - to[j];	te[j] = te[j] + to[j];
				rnd ps = pe[j] - po[j];	pe[j] = pe[j] + po[j];
				to[j] = ts;
				po[j] = ps;
			}
		#ifndef SHTNS4MAGIC
			for (int j=0; j<NWAY; ++j) {
				S2D_CSTORE2((double*)BtF, k+j, NLAT, te[j], to[j], vall(0), vall(0));
				S2D_CSTORE2((double*)BpF, k+j, NLAT, pe[j], po[j], vall(0), vall(0));
			}
		#else
			for (int j=0; j<NWAY; ++j) {
				if ((k+j)>=nk) break;
				S2D_CSTORE2_4MAGIC((double*)BtF, k+j, te[j], to[j], vall(0), vall(0));
				S2D_CSTORE2_4MAGIC((double*)BpF, k+j, pe[j], po[j], vall(0), vall(0));
			}
		#endif
			k+=NWAY;
		} while (k < nk);

	} else {	// im > 0
		const int ms = im*MRES;		// signed m.
		im = abs(im);			// positive im
		v2d VWl[llim*2+4];
		m = im*MRES;
		l = im*(2*(LMAX+1) -m);		// to compute position in NLM array.
		alm = shtns->alm + l+m;
		Slm -= m;	// virtual pointer for l=0
		Tlm -= m;

		{	// convert from vector SH to scalar SH
			double* mx = shtns->mx_stdt + l-m;
			double em = ms;
			v2d sl = ((v2d*)Slm)[m];
			v2d tl = ((v2d*)Tlm)[m];
			v2d vs = vdup( 0.0 );
			v2d wt = vdup( 0.0 );
			for (int l=m; l<=llim; l++) {
				s2d mxu = vdup( mx[2*l] );
				s2d mxl = vdup( mx[2*l+1] );		// mxl for next iteration
				vs += IxKxZ(em, tl);	// vs = addi( vs ,  em*tl );
				wt += IxKxZ(em, sl);	// wt = addi( wt ,  em*sl );
				v2d vs1 = mxl*sl;			// vs for next iter
				v2d wt1 = -mxl*tl;			// wt for next iter
				if (l<llim) {
					sl = ((v2d*)Slm)[l+1];		// kept for next iteration
					tl = ((v2d*)Tlm)[l+1];
					vs += mxu*sl;
					wt -= mxu*tl;
				}
				VWl[2*l]   = vs;
				VWl[2*l+1] = wt;
				vs = vdup( 0.0 );		wt = vdup( 0.0 );
				vs = vs1;
				wt = wt1;
			}
			VWl[2*llim+2] = vs;
			VWl[2*llim+3] = wt;
		}

		l = shtns->tm[im];
		l = ((unsigned) l) / VSIZE2;
		#ifndef SHTNS4MAGIC
		 	zero_poles4_vect(BtF, NLAT-l, BpF-BtF, 2*l);
		#else
			zero_poles2_vect(BtF, BpF-BtF, 4*l);
		#endif

		k = l;
		do {
			al = alm;
			rnd cost[NWAY], y0[NWAY], y1[NWAY];
			rnd ter[NWAY], tei[NWAY], tor[NWAY], toi[NWAY];
			rnd per[NWAY], pei[NWAY], por[NWAY], poi[NWAY];
			for (int j=0; j<NWAY; ++j) {
				cost[j] = vread(st, k+j);
				y0[j] = vall(1.0);
			}
			l=m;
			if (robert_form == 0) l=m-1;
			long int ny = 0;
		  if ((int)llim <= SHT_L_RESCALE_FLY) {
			do {		// sin(theta)^m
				if (l&1) for (int j=0; j<NWAY; ++j) y0[j] *= cost[j];
				for (int j=0; j<NWAY; ++j) cost[j] *= cost[j];
			} while(l >>= 1);
		  } else {
			long int nsint = 0;
			do {		// sin(theta)^m		(use rescaling to avoid underflow)
				if (l&1) {
					for (int j=NWAY-1; j>=0; --j) y0[j] *= cost[j];
					ny += nsint;
					if (vlo(y0[NWAY-1]) < (SHT_ACCURACY+1.0/SHT_SCALE_FACTOR)) {
						ny--;
						for (int j=NWAY-1; j>=0; --j) y0[j] *= vall(SHT_SCALE_FACTOR);
					}
				}
				for (int j=NWAY-1; j>=0; --j) cost[j] *= cost[j];
				nsint += nsint;
				if (vlo(cost[NWAY-1]) < 1.0/SHT_SCALE_FACTOR) {
					nsint--;
					for (int j=NWAY-1; j>=0; --j) cost[j] *= vall(SHT_SCALE_FACTOR);
				}
			} while(l >>= 1);
		  }
			for (int j=0; j<NWAY; ++j) {
				y0[j] *= vall(al[0]);
				cost[j] = vread(ct, j+k);
			}
			for (int j=0; j<NWAY; ++j) {
				y1[j]  = (vall(al[1])*y0[j]) *cost[j];		//	y1[j] = vall(al[1])*cost[j]*y0[j];
				por[j] = vall(0.0);		tei[j] = vall(0.0);
				tor[j] = vall(0.0);		pei[j] = vall(0.0);
				poi[j] = vall(0.0);		ter[j] = vall(0.0);
				toi[j] = vall(0.0);		per[j] = vall(0.0);
			}
			l=m;		al+=2;
			while ((ny<0) && (l<llim)) {		// ylm treated as zero and ignored if ny < 0
				for (int j=0; j<NWAY; ++j) {
					y0[j] = (vall(al[1])*cost[j])*y1[j] + vall(al[0])*y0[j];
				}
				for (int j=0; j<NWAY; ++j) {
					y1[j] = (vall(al[3])*cost[j])*y0[j] + vall(al[2])*y1[j];
				}
				l+=2;	al+=4;
				if (fabs(vlo(y0[NWAY-1])) > SHT_ACCURACY*SHT_SCALE_FACTOR + 1.0) {		// rescale when value is significant
					++ny;
					for (int j=0; j<NWAY; ++j) {
						y0[j] *= vall(1.0/SHT_SCALE_FACTOR);		y1[j] *= vall(1.0/SHT_SCALE_FACTOR);
					}
				}
			}
		  if (ny == 0) {
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
		  }
		#ifndef SHTNS4MAGIC
			for (int j=0; j<NWAY; ++j) {
				S2D_CSTORE2((double*)BtF, k+j, NLAT, ter[j], tor[j], tei[j], toi[j]);
				S2D_CSTORE2((double*)BpF, k+j, NLAT, per[j], por[j], pei[j], poi[j]);
			}
		#else
			for (int j=0; j<NWAY; ++j) {
				if ((k+j)>=nk) break;
				S2D_CSTORE2_4MAGIC((double*)BtF, k+j, ter[j], tor[j], tei[j], toi[j]);
				S2D_CSTORE2_4MAGIC((double*)BpF, k+j, per[j], por[j], pei[j], poi[j]);
			}
		#endif
			k+=NWAY;
		} while (k < nk);
	}

	#undef sr
	#undef si
	#undef tr
	#undef ti
  }

  #endif
