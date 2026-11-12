#ifndef VL_FFT_HPP
#define VL_FFT_HPP
#include <vector>
#include <map>
#include "vl_types.hpp"
#define VL_ARRAY_INCLUDES_CURVES
#define VL_ARRAY_INCLUDES_COMPLEX 
#include "vl_array.hpp"
//#include "vl_structs.hpp"

// estas estruturas devem ser substituidas por defines
#ifndef FFTPt2
template<typename DT>
struct FFTPt2
{
	DT x, y;
};
#endif
#ifndef FFTRect2
template<typename DT>
struct FFTRect2
{
	DT x, y, width, height;
};
#endif

///////////////////////////////////////////////////////////////////////////////////
// Velociraptor: FFT
//
//  <o tamanho dos vetores tem de ser potência de 2,
//   senão será necessário extender com zeros o vetor
//   para usar estas funções>
//
//  fft de vetor complexo:
//	   complexFFT/complexIFFT
//	   complexRecFFT/complexRecIFFT
//
//  fft de vetor real:
//	   realFFT/realIFFT
//     complexSymmetricExpand
//
//  FFParam:
//	   uma estrutura que guarda dois apontadores:
//			um vetor de N indices de bits "ao contrário"
//			um tabela sin/cos para N divisões do circulo
//
//	   os apontadores são conseguidos apartir dos mapas globais
//     globalRevBitsMap, globalArcTabs32/64Map, e são
//     criados senão encontrados; os indices são usados
//     para reordenar os vetores do fft iterativo,
//     as tabelas sin/cos são usadas pelo fft iterativo e 
//     recursivo
//
// -----------------------------------------------------
//  utilidades:
//
//     powerSpectrum
//        devolve um vetor que é a representação
//        da magnitude do sinal a várias frequências
//
//          para sinais complicados apenas fft sobre
//          todo o sinal não se conseguiria chegar
//          a esta representação, porque segundo as fases
//          dos sinais misturados, a magnitude tende-se a
//          espalhar nas frequencias ao redor
//
//          a solução é varrer todo sinal com uma janela,
//          multiplicando estes troços com uma função
//          que torna o sinal mais importante no centro 
//          e menos quando se afasta para os extremos 
//          (por padrao se usa a função de Hanning) - 
//          isto se faz para analizar o sinal localmente com 
//          um FFT pequeno, e agregar as magnitudes locais
//          num vetor final de magnitudes, magnitudes por frequencia;
//			a janela geralmente não se desloca amostra-a-amostra 
//          sobre o sinal mas, por um fator do seu tamanho
//          geralmente 1/2, ...1/16 do seu tamanho;
//          a isto se chama de sobreposição (overlapping)
//
//          a representação é normalizada e posta na escala 
//          10*log10 
//
//          uma das implicações desta solução, é que 
//          o tamanho da janela deve ser ajustado
//          segundo a gama de frequências que se pretende
//          representar
//			
//				por exemplo,
//              se 44100 amostras por segundo => 1Hz
//				
//              a frequência máxima é limitada pelas 44100 
//              amostras por segundo, e é de 44.1KHz
//				
//				se se pretende analizar o espetro áudio,
//              cujo intervalo audivel está nos [20Hz:20 kHz]
//              para conseguir atingir os 20Hz é 
//				é necessário 2*(44100/20) = 4410 amostras 
//				(100 milisegundos de audio)
//              o dobro pelo o teorema de amostragem de Niquist;
//				ajustado à potência de 2, seria de 8192 amostras(8K) 
//              
//				rebaixando o intervalo audivel para
//                      10Hz, 16384 amostras(16K) 
//				         5Hz, 32768 amostras(32K) 
//
//   equalizer	
//		dada uma curva de beizer ou de segmentos de reta, 
//      refaz o sinal ajustando as respetivas magnitudes
//      por das frequências, segundo a curva
//
//			se considera uma janela de tamanho ajustado
//			identico ao anterior
//
//			a anterior curva:
//           - é passada para o dominio do tempo por iFFT
//			 - se multiplica por função (padrao Blackman)
//			 - se passa para o dominio das frequencias por FFT,
//             daqui para a frente chama-se envelope
//
//			uma janela varre com sobreposição o sinal,
//			cada troço:
//			  - é passado para o dominio das frequencias
//			  -	se multiplica pelo o envelope 
//			  - é passado para o dominio do tempo
//			  - é adicionado ao vetor final do sinal, que se
//              adiciona à anterior parte do sinal que se sobrepoem
//              a atual janela


		
namespace vl{ namespace fft{

class RevBits
{
	public:
	VL_SZ N, *table;
	RevBits(VL_SZ n_): N(n_)
	{
		table = new VL_SZ[N];
		const W32 highestBit = N>>1;
		for (W32 i=0; i< N; i++)
		{
			W32 idx = 0;
			W32 bit = highestBit;
			while (bit)
			{
				idx>>=1;
				idx += (i&bit? highestBit:0);
				bit>>=1;
			}
			table[i] = idx;
		}
	}
	template<typename DT>
	// reordena para um vetor
	W0 reorderTo(Complex<DT>* from_, Complex<DT>* to_)
	{
		for (VL_SZ i=0; i < N; i++) to_[i] = from_[table[i]];
	}
};
template<typename DT>
class ArcTab
{
	public:
	VL_SZ N;
	DT    *table;
	ArcTab(VL_SZ N_): N(N_)
	{
		table = new DT[2*N];
		for(VL_SZ i = 0, j=0; j < N; i++, j+=2)
		{
			const DT ang = -(VL_M_PI2*(DT)i)/(DT)N;
			table[j]   = (DT)cos(ang);
			table[j+1] = (DT)sin(ang);
		}
	}
};

typedef std::map<VL_SZ, RevBits*>	   RevBitsMap;
typedef std::map<VL_SZ, ArcTab<FP32>*> ArcTab32Map;
typedef std::map<VL_SZ, ArcTab<FP64>*> ArcTab64Map;
//
RevBitsMap		globalRevBitsMap;
ArcTab32Map		globalArcTabs32Map;
ArcTab64Map		globalArcTabs64Map;

template<typename OBJ>
OBJ* getMapObj(VL_SZ N_, std::map<VL_SZ, OBJ*>& map_)
{
	std::map<VL_SZ, OBJ*>::iterator it = map_.find(N_);
	if (it!=map_.end()) return it->second;
	OBJ* obj = new OBJ(N_);
	map_[N_] = obj;
	return obj;
}

template<typename DT>
class FFTParam{};
template<>
class FFTParam<FP32>
{
	public:
	VL_SZ N;
	RevBits* revBitsTab;
	ArcTab<FP32>* arcTab;
	FFTParam(VL_SZ N_): N(N_)
	{
		revBitsTab = getMapObj<RevBits>(N_, globalRevBitsMap);
		arcTab = getMapObj<ArcTab<FP32>>(N_, globalArcTabs32Map);
	}
};
template<>
class FFTParam<FP64>
{
	public:
	VL_SZ N;
	RevBits* revBitsTab;
	ArcTab<FP64>* arcTab;
	FFTParam(VL_SZ N_): N(N_)
	{
		revBitsTab = getMapObj<RevBits>(N_, globalRevBitsMap);
		arcTab = getMapObj<ArcTab<FP64>>(N_, globalArcTabs64Map);
	}
}; 

#define VL_FFT_COPY_BUFFER(NAME_, X_, N__, DT_)		\
		memcpy(NAME_, X_, sizeof(Complex<DT_>)*N__);\

#define VL_FFT_INIT_BUFFER(X_, Y_, NAME_, N__, DT_)	\
	if (Y_)											\
	{												\
		memcpy(Y_, X_, sizeof(Complex<DT_>)*N__);	\
		NAME_ = Y_;									\
	}												\
	else											\
		NAME_ = X_;

// ** iterativo dif - decimation in frequency
//    entra: N elementos complexos
//      sai: N elementos complexos (ordem por bits ao avesso)
template<typename DT>
W0 complexFFT(Complex<DT> *x_, VL_SZ N_, Complex<DT>* y_=0)
{
	#define VL_FFT_UNROLL_COMPLEX
	#include "unroll/vl_fft_unroll.hpp"
}
template<typename DT>
W0 complexIFFT(Complex<DT> *x_, VL_SZ N_, Complex<DT>* y_=0)
{
	#define VL_IFFT_UNROLL_COMPLEX
	#include "unroll/vl_fft_unroll.hpp"
}

// ** complexo recursivo
//    entra: N elementos complexos
//      sai: N elementos complexos (por ordem)
template<typename DT>
W0 complexRecFFT_(Complex<DT>* x_,  VL_SZ N_, Complex<DT>* y_, FFTParam<DT>& h_, VL_SZ s_ = 1)
{
	#define VL_FFT_UNROLL_COMPLEX_REC
	#include "unroll/vl_fft_unroll.hpp"
}
template<typename DT>
W0 complexRecIFFT_(Complex<DT>* x_, VL_SZ N_, Complex<DT>* y_,  FFTParam<DT>& h_, VL_SZ s_ = 1)
{
	#define VL_IFFT_UNROLL_COMPLEX_REC
	#include "unroll/vl_fft_unroll.hpp"
}
template<typename DT>
W0 complexRecFFT(Complex<DT>* x_,  VL_SZ N_, Complex<DT>* y_=0)
{
	FFTParam<DT> h(N_);
	if (y_)
		complexRecFFT_<DT>(x_, N_, y_, h);
	else
		complexRecFFT_<DT>(x_, N_, x_, h);
}
template<typename DT>
W0 complexRecIFFT(Complex<DT>* x_, VL_SZ N_, Complex<DT>* y_=0)
{
	FFTParam<DT> h(N_);
	Complex<DT>* buffer;
	buffer = (y_)? y_: x_;
	complexRecIFFT_<DT>(x_, N_, buffer, h);
	const DT sc = (DT)1.0f/(DT)N_;
	for (VL_SZ i=0; i< N_; i++)
	{
		buffer[i].re*=sc;
		buffer[i].im*=sc;
	}
}

// ** iterativo real 
//
// - é possível fazer o DFT dum vetor real, tratando o vetor de N
//   elementos reais como um vetor de N/2 elementos complexos,
//   se após se fizer uma tranformação que é invertivel
//
// - dado que o resultado dum DFT dum vetor real - N complexos - 
//   tem simetria conjugada após o elemento N/2 [N/2+1:N-1] e
//   os elementos 0 e N/2 são reais se pode guardar ambos
//   no elemento 0, e assim compactar a representação
//   para um vetor de N/2 complexos
//
//	pelo que aqui se usa a transformação entre:
//		H (fft de N/2 complexos, vetor de N/2 complexos)
//		F (fft de N reais, vetor de N/2 complexos, compacto)
//
//  se depois for necessário o resultado verdadeiro do DFT real, usar 
//  complexSymmetricExpand tendo como entrada o vetor complexo compacto

//
// ** realFFT 
//      entra: N elementos reais
//        sai: N/2 elementos complexos (compacto, por ordem)
//         
template<typename DT>
W0 realFFT(DT *x__, VL_SZ N_, Complex<DT>* y_=0)
{
	const VL_SZ No2 = N_>>1; // numero de elementos como complexo

	//fazer H=fft(x_) como se vetor complexo (de N/2 complexos)
	//para temp
	Complex<DT> *x_ = (Complex<DT>*)(x__);
	Complex<DT> *temp = new Complex<DT>[No2];
	complexFFT<DT>(x_, No2, temp);

	//reordernar para x_ ou y_
	FFTParam<DT> h0(No2);
	Complex<DT> *buffer = y_? y_:x_;
	h0.revBitsTab->reorderTo(temp, buffer);

	// apagar temp
	delete[] temp;

	FFTParam<DT> h(N_);
	// transformação H->F
	//  H (fft de N/2 complexos, vetor N/2 complexos)
	//   para
	//  F (fft de N reais, vetor N/2 complexos compacto** ver acima)
	#define VL_FFT_UNROLL_REAL_TRANSFORM
	#include "unroll/vl_fft_unroll.hpp"
}
// ** realIFFT
//      entra: N/2 elementos complexos (compacto, por ordem)
//        sai: N elementos reais (por ordem)
//
template<typename DT>
W0 realIFFT(Complex<DT>* x_, VL_SZ N_, DT* y_=0)
{
	const VL_SZ No2 = N_>>1; // numero de elementos como complexo
	FFTParam<DT> h(N_);	

	Complex<DT> *buffer;
	VL_FFT_INIT_BUFFER(x_,((Complex<DT>*)y_),buffer, No2, DT)

	// transformação F->H
	//  F (fft de N reais, vetor N/2 complexos compacto** ver acima)
	//   para
	//  H (fft de N/2 complexos, vetor N/2 complexos)
	#define VL_IFFT_UNROLL_REAL_TRANSFORM
	#include "unroll/vl_fft_unroll.hpp"

	Complex<DT> *temp = new Complex<DT>[No2];

	//fazer ifft(H)
	complexIFFT<DT>(buffer, No2, temp);

	//reordenar
	FFTParam<DT> h0(No2);
	h0.revBitsTab->reorderTo(temp, buffer);

	// apagar temp
	delete[] temp;
}


// ** pega num vetor compacto resultado de realFFT
//    e completa a simetria, resultando no que seria
//    o verdadeiro DFT resultado dum DFT real
//
//      entra: N/2 elementos complexos (compacto, por ordem)
//        sai: N elementos complexos (por ordem)
//
//    - se shiftZeroToCenter_=0, é o vetor que seria o 
//                               verdadeiro resultado 
//                               de realFFT
//
//    - se shiftZeroToCenter_=1, os elementos são empurrados
//       para a frente, e a simetria é feita atrás deste, 
//       pelo que o elemento que representa a frequência 
//       mais baixa (DC) passa ao centro e os que representam
//       as altas nos ficarão nos extremos
//
template<typename DT>
W0 complexSymmetricExpand(Complex<DT>* x_, VL_SZ N_, Complex<DT>* y_, W1 shiftZeroToCenter_=true)
{
	if (shiftZeroToCenter_)
	{
		memcpy(y_+N_+1, x_+1, (N_-1)*sizeof(Complex<DT>));
		for (VL_SZ i=1; i < N_; i++)
		{
			y_[N_-i] = x_[i];
		}
		y_[N_].re = x_[0].re;
		y_[N_].im = 0;
		y_[0].re = x_[0].im;
		y_[0].im = 0;
	}
	else	
	{
		memcpy(y_+1, x_+1, (N_-1)*sizeof(Complex<DT>));
		for (VL_SZ i=1; i < N_; i++)
		{
			y_[N_+i] = x_[N_-i];
		}
		y_[0].re = x_[0].re;
		y_[0].im = 0;
		y_[N_].re = x_[0].im;
		y_[N_].im = 0;
	}
}

template<typename DT>
DT ceilPower2(DT n_)
{
	const DT n = n_;
	if (n_==0) return 0;
	SW32 i=-1;
	while(n_) { n_>>=1; i++; }
	if ((DT)(((W64)1)<<i) == n) return n;
	return (DT)(1<<(i+1));
}


VL_SZ getWindowSize(VL_SZ sampleRate_, VL_SZ minFreq_)
{	
	// a largura da janela necessária sobre um sinal real
	// para efetuar fft no intervalo de frequências
	// F = [minFreq_ : sampleRate_]
	return 2*(vl::fft::ceilPower2(sampleRate_/minFreq_));
}
VL_SZ getFrequencyEnvelopeSize(VL_SZ winSz_)
{	
	// a largura do envelope de frequencias 
	// (complexo compacto) é metade 
	return (winSz_>>1);
}
VL_SZ getSpectrumSize(VL_SZ winSz_)
{	
	// a largura do espetro de dada pela a janela
	// (real) é metade + 1 
	// para incluir a frequencia mais alta
	// F = [minFreq_ : sampleRate_]
	return (winSz_>>1)+1;
}

// ** powerSpectrum 
//      entra:     x_: sinal de N_ elementos reais
//             winSz_: largura da janela 
//                     (alargada para próxima potência de 2)
//         overlapSz_: largura da sobreposição
//
//        sai:     y_: com uma largura (winSz_/2)+1
//        
template<typename DT>
W0 powerSpectrum(DT* x_, 
				 VL_SZ N_, 
				 DT* y_, 
				 VL_SZ winSz_ = 1<<14,  //16K amostras
				 VL_SZ overlapSz_ = 0,  //metade
				 W1 toDecibels_ = true,
				 vl::curves::CurveId curveId_  = vl::curves::curveHamming
				)
{
	winSz_ = ceilPower2(winSz_);
	if (overlapSz_==0) overlapSz_ = winSz_>>1;

	const VL_SZ shiftSz = winSz_ - overlapSz_;
	const VL_SZ spectrumSz = getSpectrumSize(winSz_);

	DT *buffer = new DT[winSz_];
	memset(y_, 0, sizeof(DT)*spectrumSz);

	for (VL_SZ i=0; i< N_-shiftSz; i+=shiftSz)
	{
		// Copiar sinal da janela
		const VL_SZ toCopy = VL_MIN(winSz_, N_-i);
		memcpy(buffer, x_ + i, toCopy * sizeof(DT));

		if (toCopy < winSz_) // zerar se sobrar
			for (VL_SZ j=toCopy; j< winSz_; j++) buffer[j]=0;

		// Adaptar janela à função
		vl::curves::mulFunction<DT>(curveId_, buffer, winSz_);

		// FFT sobre
		realFFT(buffer, winSz_); 

		// Adicionar as magnitudes ao vetor
		Complex<DT>* b = (Complex<DT>*)buffer;
		for (VL_SZ j=1; j< spectrumSz-1; j++)
			y_[j] += b[j].re * b[j].re + b[j].im * b[j].im;
		// ...tratar a parte compactada
		y_[0] += b[0].re * b[0].re;
		y_[spectrumSz-1] += b[0].im * b[0].im;
	}

	// Normalizar e escala dB
	const DT windowsNumber = (DT)N_/shiftSz;
	const DT n = (DT)1.0f/(windowsNumber*(DT)winSz_);

	if (toDecibels_)
	{
		for (VL_SZ j=0; j< spectrumSz; j++)	
		{
			const DT res = n * y_[j];
			y_[j] = (res>0)? 10.0f * std::log10f(res) : 0;
		}
	}
	else
	{
		for (VL_SZ j=0; j< spectrumSz; j++)	
		{
			const DT res = n * y_[j];
			y_[j] = (res>0)? res: 0;
		}
	}

	delete[] buffer;
}

// ** equalização
//
//	genFrequencyEnvelope/genFrequencyDirectEnvelope
//
//		a largura da janela WinSz_ deve ser decidida
//		consoante a gama de frequencias a analizar
//
//		constroi o envelope que sai em y_ e tem o tamanho dado por 
//      getSpectrumSize(WinSz_)
//
//		os pontos <set_> são normalizados pelo o quadro <frame_>
//      em que se inserem 
//         - a coordenada X é normalizada à largura do vetor <env_>,
//         - a coordenada Y:
//              - se omitido <normalizeYrange_>, a altura de cada ponto
//                é a usada 
//              - senão é normalizada pelo interior do quadro, e 
//                se assume ser dum ponto gráfico em que altura representa
//		   	      amplitude, pelo que se usa ao contrário em altura
//
//		pointsAreBars_=true  cada barra é dada pelo canto
//                           topo direita (TR); o ultimo
//							 ponto pode estar em qualquer X,
//                           mas a sua altura continua constante
//                           como o que acontece a (3) na figura:
//										 _____x___
//                           .-----X     |        |  
//							 |  1  |-----X    3   |
//							 |     |  2  |        |
//		pointsAreBars_=false cada 2 pontos são interpretados
//                           como um segmento de reta
//	                         antes de 1º ponto      = constante	desde o inicio
//	                         depois de ultimo ponto = constante	até ao fim
//	                         entre pontos           = segmento de reta
//
//		pointsAreOrdered_	 evita ordenar os pontos se já
//		curveId_             a função a multiplicar
//

//#define VL_FFT_DEBUG_FREQ_ENV

//	genFrequencyEnvelope:
//		env_: real, só tem as amplitudes por frequência [0:N/2-1]
template<typename DT, typename SDT>
W0 genFrequencyEnvelope(std::vector<FFTPt2<SDT>>& set_,
						FFTRect2<SDT> frame_,
						VL_SZ winSz_,
						DT *env_, 
						FP32 scaleY_ = 1.0f,
						W1 normalizeYrange_ = false,
						W1 pointsAreBars_ = false,
						W1 pointsAreOrdered_ = false,
						SW32 useOneBar = -1)
{
	if (normalizeYrange_)
	{
		#define VL_FFT_UNROLL_CREATE_ENVELOPE_AMPS
		#include "unroll\vl_fft_unroll.hpp"
	}
	else
	{
		#define VL_FFT_UNROLL_CREATE_ENVELOPE_NORMY_NO
		#define VL_FFT_UNROLL_CREATE_ENVELOPE_AMPS
		#include "unroll\vl_fft_unroll.hpp"
	}

	#ifdef VL_FFT_DEBUG_FREQ_ENV
		cv::namedWindow("curva envelope",0);
		cv::Mat im(cv::Size(640,380), CV_8UC3);
		//VL_FP mm[] = {0,1.0f,1.0f};
		vl::draw::drawCurve(im, env_, envSz, cv::Scalar(0), 1);
		cv::imshow("curva envelope", im);
	#endif
}

//#define VL_FFT_DEBUG_FREQ_ENV
//	genFrequencyDirectEnvelope: 
//						<!><!><!><!><!> Não parece funcionar!!!
//						testado em 17/03/2022
//		env_: complexo compacto, foi passado para o tempo e
//            multiplicado no tempo por curveId_ e passado
//            novamente para o dominio das frequencias
//
template<typename DT, typename SDT>
W0 genFrequencyDirectEnvelope(std::vector<FFTPt2<SDT>>& set_,
							  FFTRect2<SDT> frame_,
							  VL_SZ winSz_,
							  Complex<DT> *env_, 
							  W1 normalizeYrange_ = false,
							  W1 pointsAreBars_ = false,
							  W1 pointsAreOrdered_ = false, 
							  SW32 useOneBar = -1,
							  vl::curves::CurveId curveId_ = 
							  	vl::curves::curveBlackman)
{
	
	const VL_SZ envSz = getFrequencyEnvelopeSize(winSz_); 
	if (normalizeYrange_)
	{
		#define VL_FFT_UNROLL_CREATE_ENVELOPE_COMPLEX
		#include "unroll\vl_fft_unroll.hpp"
	}
	else
	{
		#define VL_FFT_UNROLL_CREATE_ENVELOPE_NORMY_NO
		#define VL_FFT_UNROLL_CREATE_ENVELOPE_COMPLEX
		#include "unroll\vl_fft_unroll.hpp"
	}

	// o envelope foi feito para a largura envSz=winSz/2
	// logo a amplitude para a frequencia mais alta
	// <envSz> não foi atribuida, vou assumir que é igual  
	// à do ultimo elemento (envSz-1) e compactar
	env_[0].im = env_[envSz-1].re;
	
	// passar para dominio tempo, vetor real de tamanho duplo
	const VL_SZ envRealSz = envSz<<1;
	vl::fft::realIFFT(env_, envRealSz);
	DT* envReal = (DT*)env_;

	#ifdef VL_FFT_DEBUG_FREQ_ENV
		std::string wt1 = "tempo",
					wt2 = "tempo trocado",
					wt3 = "adaptado inv, 1/4 zero restante",
					wt4 = "funcao";
		cv::namedWindow(wt1,0);
		cv::namedWindow(wt2,0);
		cv::namedWindow(wt3,0);
		cv::namedWindow(wt4,0);

		cv::Mat im1(cv::Size(640,380), CV_8UC3);
		cv::Mat im2(cv::Size(640,380), CV_8UC3);
		cv::Mat im3(cv::Size(640,380), CV_8UC3);
		cv::Mat im4(cv::Size(640,380), CV_8UC3);

		//vl::draw::drawCurve(im1, envReal, envRealSz, cv::Scalar(0), 1);
		//DT mm[2];
		//vl::arr::minmax<DT>(envReal, envRealSz, mm, mm+1);
	#endif


	// agora é real, copiar extremos
	W32 envRealSz2 = envRealSz>>1;
	W32 envRealSz4 = envRealSz>>2;
	W32 j;
	for (j=0; j< envRealSz4; j++)
	{
		envReal[envRealSz4+j] = envReal[j];
		envReal[j]			  = envReal[envRealSz-envRealSz4+j];
	}
	// e zerar restante
	for (j=envRealSz2; j< envRealSz; j++) envReal[j] = 0;

	#ifdef VL_FFT_DEBUG_FREQ_ENV
		//vl::draw::drawCurve(im2, envReal, envRealSz, cv::Scalar(0), 1);
	#endif 

	// adaptar com a função ao contrario (blackman padrao)
	DT* envCurve = new DT[envRealSz2];
	vl::curves::makeFunction<DT>(curveId_, envCurve, envRealSz2, true);
	for (j=0; j< envRealSz2; j++)
		envReal[j] *= envCurve[j];
		
	#ifdef VL_FFT_DEBUG_FREQ_ENV
		//vl::draw::drawCurve(im3, envReal, envRealSz, cv::Scalar(0), 1);
		//vl::draw::drawCurve(im4, envCurve, envRealSz2, cv::Scalar(0), 1);
	#endif 
	
	#ifdef VL_FFT_DEBUG_FREQ_ENV
		cv::imshow(wt1,im1);
		cv::imshow(wt2,im2);
		cv::imshow(wt3,im3);
		cv::imshow(wt4,im4);
		//cv::waitKey(0);
	#endif

	delete[] envCurve;

	// passar para dominio frequências
	vl::fft::realFFT(envReal, envRealSz);

	// DEBUG: envelope = 1 equalização tem de ficar parecida
	//for (W32 i=0; i< envSz; i++)
	//{
	//	env_[i].re = 1.0f;
	//	env_[i].im = 1.0f;
	//}
	//env_[0].im = env_[envSz-1].re;
}

// Filtra um sinal x_ para y_ segundo um envelope 
//		
//	 O varrimento é feito de forma a sobrepor as janelas de
//	 largura winSz_ pela largura overlapSz_
//
//   a largura do envelope deve corresponder a winSz_
//
//	 equalizationDirect:
//		multiplica a janela no dominio das frequencias 
//      pelo o envelope complexo gerado antes por 
//      "genFrequencyDirectEnvelope"
//
//   equalization:
//		mais lento, multiplica a janela no dominio das frequencias 
//		pelo o envelope real gerado antes por 
//		"genFrequencyEnvelope", e depois no dominio do
//      tempo multiplica pelo envelope padrao Barllet (triangulo)
//		este ultimo envelope é do tamanho da janela <WinSz_> e pode
//		ser por entregue <signalEnv_> se já foi criado
//

template<typename DT>
W0 equalizationDirect
			  (DT* x_, 
				VL_SZ N_, 
				DT* y_, 
				Complex<DT>* freqEnvelope_,
				VL_SZ winSz_,
				VL_SZ overlapSz_)
{
	winSz_ = ceilPower2(winSz_);

	// inicio da iteração da janela
	#define VL_FFT_UNROLL_EQUALIZATION_BEGIN
	#include "unroll\vl_fft_unroll.hpp"

		// FFT local
		realFFT(buffer, winSz_); 
		Complex<DT> *b = (Complex<DT>*) buffer;
		
		// multiplicar por envelope
		// o envelope é real na frequencia mais baixa e na
		// mais alta; assumir que frequencia mais alta vai levar
		// a mm amplitude que ultimo elemento e compactar

		b[0].re *= freqEnvelope_[0].re;
		b[0].im *= freqEnvelope_[envSz-1].re;
		for (VL_SZ j=1; j< envSz; j++)
		{
			b[j].re *= freqEnvelope_[j].re;
			b[j].im *= freqEnvelope_[j].im;
		}
		
		// IFFT local
		realIFFT(b, winSz_); 
		
	// fim da iteração da janela
	#define VL_FFT_UNROLL_EQUALIZATION_END
	#include "unroll\vl_fft_unroll.hpp"
}

template<typename DT>
W0 equalization(DT* x_, 
				 VL_SZ N_, 
				 DT* y_, 
				 DT* ampEnvelope_,
				 VL_SZ winSz_,
				 VL_SZ overlapSz_, 
				 DT* signalEnv_ = 0,
					vl::curves::CurveId curveId_ = 
					vl::curves::curveBartlett)
{
	winSz_ = ceilPower2(winSz_);
	
	// blackman
	DT* signalEnv;

	if (!signalEnv_) 
	{
		// criar a função (blackman padrao)
		signalEnv = new DT[winSz_];
		vl::curves::makeFunction<DT>(curveId_, signalEnv, winSz_);
	}
	else
		signalEnv = signalEnv_;

	// inicio iteração da janela
	#define VL_FFT_UNROLL_EQUALIZATION_BEGIN
	#include "unroll\vl_fft_unroll.hpp"

		// FFT local
		realFFT(buffer, winSz_); 
		Complex<DT> *b = (Complex<DT>*) buffer;
		
		// envelope frequencias
		b[0].re *= ampEnvelope_[0];
		b[0].im *= ampEnvelope_[envSz-1];
		for (VL_SZ j=1; j< envSz; j++)
		{
			b[j].re *= ampEnvelope_[j];
			b[j].im *= ampEnvelope_[j];
		}
		
		// IFFT local
		realIFFT(b, winSz_); 

		// envelope sinal
		for (W32 j=0; j< winSz_; j++)
			buffer[j]*=signalEnv[j];
		
	// fim da iteração da janela
	#define VL_FFT_UNROLL_EQUALIZATION_END
	#include "unroll\vl_fft_unroll.hpp"

	if (!signalEnv_)
		delete[] signalEnv;
}

};


};

#undef VL_FFT_COPY_BUFFER
#undef VL_FFT_INIT_BUFFER
#endif