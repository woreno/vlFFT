
// A diferença entre o FFT e o IFFT está no sinal
// utilizado no somatório no fator W=e(+/-(2TTjkn)/N)

#if defined(VL_FFT_UNROLL_COMPLEX)||defined(VL_IFFT_UNROLL_COMPLEX)

	FFTParam<DT> h(N_);
	vl::Complex<DT>* buffer;
	VL_FFT_INIT_BUFFER(x_, y_, buffer, N_, DT)

	W32 bpg = N_>>1;
	W32 dk = 1;
	const vl::Complex<DT>* const e1 = buffer + N_;
	while(bpg > 0)
	{
		//printf("\n-stage %i-",stage);
		Complex<DT> *A = buffer;
		Complex<DT> *B = buffer + bpg;
		while(A < e1)
		{
			const Complex<DT>* W = (Complex<DT>*)(h.arcTab->table);
			const Complex<DT>* const e2 = B;
			while(A < e2)
			{
				const Complex<DT> b = B[0];
				const Complex<DT> Z = A[0]-b;
				B[0] = Z*
					#ifdef VL_FFT_UNROLL_COMPLEX
						W[0];
					#else
						W[0].conj();
					#endif
				A[0] = (A[0] + b);
				A++; B++; W+=dk;
			}
			A = B;
			B += bpg;
		}
		bpg >>= 1;
		dk<<=1;
	}

	// normalizar se IFFT
	#ifdef VL_IFFT_UNROLL_COMPLEX
		Complex<DT> *A = buffer;
		const DT sc = (DT)1.0f/(DT)N_;
		while(A < e1)
		{
			A[0].re*=sc;
			A[0].im*=sc;
			A++;
		}		
	#endif

	#undef VL_FFT_UNROLL_COMPLEX
	#undef VL_IFFT_UNROLL_COMPLEX
#endif

#if defined(VL_FFT_UNROLL_COMPLEX_REC)||defined(VL_IFFT_UNROLL_COMPLEX_REC)
	if (N_==1)
	{
		y_[0] = x_[0];
		return;	
	}
	const VL_SZ s2 = s_<<1;
	const VL_SZ N2 = N_>>1;

	#ifdef VL_FFT_UNROLL_COMPLEX_REC
		complexRecFFT_(x_,    N2, y_   , h_, s2);
		complexRecFFT_(x_+s_, N2, y_+N2, h_, s2);
	#else
		complexRecIFFT_(x_,    N2, y_   , h_, s2);
		complexRecIFFT_(x_+s_, N2, y_+N2, h_, s2);
	#endif

	const Complex<DT>* W = (Complex<DT>*)(h_.arcTab->table);
	for (VL_SZ k=0; k< N2; k++)
	{
		//const Complex<DT> W = Complex::cis(-(2*VL_M_PI2)*((FP)k/N));
		const Complex<DT> z = y_[N2+k]*
			#ifdef VL_FFT_UNROLL_COMPLEX_REC
				W[0];
			#else
				W[0].conj();
			#endif
		y_[N2+k] = y_[k]-z;
		y_[k] = y_[k]+z;
		W += s_;
	}
	#undef VL_FFT_UNROLL_COMPLEX_REC
	#undef VL_IFFT_UNROLL_COMPLEX_REC
#endif

// Transformação necessária ao FFT real
#if defined(VL_FFT_UNROLL_REAL_TRANSFORM)||defined(VL_IFFT_UNROLL_REAL_TRANSFORM)

	const VL_SZ N = h.N;
	Complex<DT> *A = buffer + 1;
	Complex<DT> *B = buffer + No2 - 1;
	const Complex<DT>* W = ((Complex<DT>*)(h.arcTab->table)) + 1;
	while(A<=B)
	{
		const Complex<DT> M = A[0] + B[0].conj();
		const Complex<DT> m = A[0] - B[0].conj();
		const Complex<DT> V = m * Complex<DT>(W[0].im, 
			#ifdef VL_FFT_UNROLL_REAL_TRANSFORM
				-W[0].re);
			#else	
				 W[0].re);
			#endif
		A[0] = (M + V) * (DT)(0.5f) ;
		B[0] = (M - V).conj() * (DT)(0.5f);
		W++; A++; B--;
	}

	#ifdef VL_FFT_UNROLL_REAL_TRANSFORM
		// compactar o elemento 0 e N/2 que são reais, 
		// nas primeiras posições (reais) porque não 
		// existe espaço no vetor
		const DT temp0 = buffer[0].re+buffer[0].im;
		const DT temp1 = buffer[0].re-buffer[0].im;
		buffer[0].re = temp0;
		buffer[0].im = temp1;
	#else
		const DT temp1 = (DT)(0.5f)*(buffer[0].re-buffer[0].im);
		buffer[0].im = temp1;
		buffer[0].re -= temp1;
	#endif

	#undef VL_FFT_UNROLL_REAL_TRANSFORM
	#undef VL_IFFT_UNROLL_REAL_TRANSFORM
#endif

// Equalização
#if defined(VL_FFT_UNROLL_CREATE_ENVELOPE_COMPLEX)||defined(VL_FFT_UNROLL_CREATE_ENVELOPE_AMPS)


	if (set_.size()==0) return;

	//if (!pointsAreOrdered_)
	//	vl::polys::sortPointsByAxis(set_, 0);

	// criar vetor segundo os pontos
	const VL_SZ envSz = getFrequencyEnvelopeSize(winSz_); 
	
	const DT xm[] = { (DT)frame_.x,  (DT)frame_.width / (DT)(envSz - 1) },
			 ym[] = { (DT)frame_.y, (DT)frame_.height};
 
	#ifdef VL_FFT_UNROLL_CREATE_ENVELOPE_NORMY_NO
		#define VL_FFT_NORM_PT(P_)  ( (DT)(P_.x-xm[0])/xm[1]), \
									  (P_.y)
	#else
		#define VL_FFT_NORM_PT(P_)  ( (DT)(P_.x-xm[0])/xm[1]), \
									 ((DT)(P_.y-ym[0])/ym[1])

		#define VL_FFT_CHECK_PT(P_) \
				if (P_.x < 0) P_.x = 0; else if (P_.x >= envSz) P_.x = (envSz - 1);\
				if (P_.y <0) P_.y=0; else if (P_.y > 1) P_.y = 1;

	#endif

	// fazer o envelope de frequencias
	// segundo os pontos 
	if (pointsAreBars_)
	{
		W32 dstX0 = 0, dstX1 = 0;
		for (VL_SZ i=0; i< set_.size(); i++)
		{
			FFTPt2<DT> p0(  VL_FFT_NORM_PT(set_[i]) );
			VL_FFT_CHECK_PT(p0);

			DT t;
			if ((useOneBar == -1) || (useOneBar == i))
			{
				#ifdef VL_FFT_UNROLL_CREATE_ENVELOPE_NORMY_NO
					t = scaleY_ * p0.y;
				#else
					t = scaleY_ * (1.0f - p0.y);
				#endif
			}
			else
			{
				t = 0;
			}

			if (i == set_.size() - 1)
				dstX1 = envSz - 1;
			else
				dstX1 = p0.x;
			for (W32 j= dstX0; j<= dstX1; j++)
			{
				#ifdef VL_FFT_UNROLL_CREATE_ENVELOPE_COMPLEX
					env_[j].re = t;
					env_[j].im = 0;
				#else 
					env_[j] = t;
				#endif
			}
			dstX0 = dstX1;
		}
	}
	else
	{
		FFTPt2<DT> p0(  VL_FFT_NORM_PT(set_[0]) );
		VL_FFT_CHECK_PT(p0);
		W32 dstX0 = p0.x, dstX1 = 0;
		DT magn;

		#ifdef VL_FFT_UNROLL_CREATE_ENVELOPE_NORMY_NO
			magn = scaleY_*p0.y;
		#else
			magn = scaleY_*(1.0f - p0.y);
		#endif

		for (VL_SZ i = 0; i < dstX0; i++)
		{
			#ifdef VL_FFT_UNROLL_CREATE_ENVELOPE_COMPLEX
				env_[i].re = magn;
				env_[i].im = 0;
			#else
				env_[i] = magn;
			#endif
		}

		for (VL_SZ i=1; i< set_.size(); i++)
		{
			FFTPt2<DT> p1(  VL_FFT_NORM_PT(set_[i]) );
			VL_FFT_CHECK_PT(p1);

			if (set_[i-1].x!=set_[i].x) 
			{
				const FFTPt2<DT> v = p1-p0;
				const DT m = (v.y==0)? 0: v.y/v.x;

				dstX1 = (W32)p1.x;
				#ifdef VL_FFT_UNROLL_CREATE_ENVELOPE_NORMY_NO
					magn = scaleY_*p0.y;
				#else
					magn = scaleY_*(1.0f - p0.y);
				#endif
				for (W32 j=0; j<= (dstX1 - dstX0); j++, 
						#ifdef VL_FFT_UNROLL_CREATE_ENVELOPE_NORMY_NO
							magn+=m
						#else
							magn-=m
						#endif
					)
				{
					const VL_SZ idx = dstX0+j;
					#ifdef VL_FFT_UNROLL_CREATE_ENVELOPE_COMPLEX
						env_[idx].re = magn;
						env_[idx].im = 0;
					#else
						env_[idx] = magn;
					#endif
				}
			}
			dstX0 = dstX1;
			p0 = p1;
		}
		#ifdef VL_FFT_UNROLL_CREATE_ENVELOPE_NORMY_NO
			magn = scaleY_*p0.y;
		#else
			magn = scaleY_*(1.0f - p0.y);
		#endif
		for (VL_SZ i = dstX0; i < envSz; i++)
		{
			#ifdef VL_FFT_UNROLL_CREATE_ENVELOPE_COMPLEX
				env_[i].re = magn;
				env_[i].im = 0;
			#else
				env_[i] = magn;
			#endif
		}

	}

	#undef VL_FFT_NORM_PT
	#undef VL_FFT_UNROLL_CREATE_ENVELOPE_NORMY_NO
	#undef VL_FFT_UNROLL_CREATE_ENVELOPE_COMPLEX
	#undef VL_FFT_UNROLL_CREATE_ENVELOPE_AMPS
#endif 


#ifdef VL_FFT_UNROLL_EQUALIZATION_BEGIN 
	const VL_SZ shiftSz_ = winSz_ - overlapSz_;
	const VL_SZ envSz = getFrequencyEnvelopeSize(winSz_); 
	DT *buffer = new DT[winSz_];
	
	// varrer o sinal por deslocamentos de <shiftSz_>
	for (VL_SZ i=0; i< N_-shiftSz_; i+=shiftSz_)
	{
		// Copiar sinal da janela
		const VL_SZ toCopy = VL_MIN(winSz_, N_-i);
		memcpy(buffer, x_ + i, toCopy * sizeof(DT));

		if (toCopy < winSz_) // zerar se sobrar
			for (VL_SZ j=toCopy; j< winSz_; j++) buffer[j]=0;

#undef VL_FFT_UNROLL_EQUALIZATION_BEGIN
#endif
#ifdef VL_FFT_UNROLL_EQUALIZATION_END 
		// Sobrepor e copiar
		if (i==0)
		{
			// na 1ª iteração todo a janela é copiada
			memcpy(y_+i, buffer, sizeof(DT)*toCopy);
		}
		else
		{
			// nas seguintes iterações, a parte sobreposta
			// é adicionada ao vetor, e a restante é 
			// colocada à frente

			// mas se a janela está chegando ao fim
			// o que resta para copiar fica menor
			// que intervalo de sobreposição pelo que se vai
			// fechando

			if (overlapSz_>toCopy)
				overlapSz_ = toCopy;

			// adicionar ao anterior
			VL_SZ j;
			for (j=0; j<overlapSz_; j++)
				y_[i+j] += buffer[j];

			// colocar restante
			for (; j<toCopy; j++)
				y_[i+j]  = buffer[j];
		}
	}
	delete[] buffer;
	#undef VL_FFT_UNROLL_EQUALIZATION_END
#endif

