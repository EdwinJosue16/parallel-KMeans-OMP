#include "KMeansP.h"
#include <math.h> 
#include "Aleatorizador.h"
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <omp.h>
#define EPSILONS 1

using namespace std;

KMeansP::KMeansP() {
	//constructor vacio
}

KMeansP::KMeansP(int n, int m, int k, double l, vector < vector<double> > & datos, int thread_count, double epsilon) {
	this->n = n; // dimension de vectores R^n
	this->m = m; //cantidad de elementos
	this->k = k;
	this->l = l*k;
	this->epsilon_p=epsilon;
	this->thread_count = thread_count;
	this->U = datos;
	U.resize(m);
	llenarX();
}

int KMeansP::pos_del_minimo(vector<double>&vec) {
	int pos_min = 0;
	double valor = vec[0];
	for (unsigned int i = 1; i < vec.size(); ++i) {
		if (vec[i] < valor) {
			valor = vec[i];
			pos_min = i;
		}
	}
	return pos_min;
}

void KMeansP::llenarX() {
	X.resize(m);
	for (int i = 0; i < m; ++i) {
		X[i].setId(i);
		X[i].setPeso(0); // llena el espacio de elementos X
	}
}

double KMeansP::distanciaEuclidea(Elemento & e1, Elemento & e2, vector< vector <double> >& U1, vector< vector <double> >& U2) {//ok
	double distancia = 0.0;
	double diferencia = 0.0;
	vector <double> vec1= U1[e1.getId()];
	vector <double> vec2= U2[e2.getId()];
	for (int i = 0; i < n; ++i) {
		diferencia = vec1[i]-vec2[i]; // entrada Vi-Ui
		distancia += diferencia * diferencia;
		// U[e1.getId()] es el vector que representa a e1
		// U[e2.getId()] es el vector que representa a e2
	}
	distancia = sqrt(distancia);
	return distancia;
}

double KMeansP::distanciaElementoConjunto(Elemento & elem_x,vector<Elemento> & conjY, vector < vector <double> > & Ux, vector < vector <double> > & Uy) {//ok
	double distancia = distanciaEuclidea(elem_x, conjY[0], Ux, Uy); //asumimos que la distancia minima es de x vs el primer elemento de Y
	int cardY = conjY.size();
	double distanciaAct = 0.0;
	for (int i = 1; i < cardY; ++i) {
		distanciaAct = distanciaEuclidea(elem_x, conjY[i], Ux, Uy);
		if (distanciaAct<distancia) {
			distancia = distanciaAct;
		}
	}
	return distancia; 
}

double KMeansP::costoXRespectoY(vector <Elemento> & conjX , vector<Elemento> & conjY, vector< vector <double> >& Ux, vector< vector <double> >& Uy) { //ok
	double costo = 0.0;
	double distanciaElemX2Y = 0.0;
	int cardX = conjX.size();
	int i;
#	pragma omp parallel num_threads(thread_count) private(distanciaElemX2Y,i) reduction(+:costo)
	{
		#	pragma omp for schedule(static) 
		for (i = 0; i < cardX; ++i) {
			distanciaElemX2Y = distanciaElementoConjunto(conjX[i], conjY, Ux, Uy);
			costo += distanciaElemX2Y * distanciaElemX2Y;
		}
	}
	return costo;
}

int KMeansP::qEMCy_i(vector<Elemento>& X, Elemento& yi, vector < vector<double> >& Ux, vector < vector<double> >& Uyi) {
	double d_actual = 0.0;
	double d_min = distanciaEuclidea(X[0], yi, Ux, Uyi);
	int posDelMC = 0;
	int cardX = X.size();
	for (int j = 1; j < cardX; ++j) {
		d_actual= distanciaEuclidea(X[j], yi, Ux, Uyi);
		if (d_actual < d_min) {
			d_min = d_actual;
			posDelMC = j;
		}
	}
	return posDelMC;
}// metodo para saber cual es el elemento en X más cercano a yi (yi está en Y)

int KMeansP::centroideMC(vector<Elemento>& cent, Elemento& yi, vector < vector<double> >& Uc, vector < vector<double> >& Uyi) {
	double d_actual = 0.0;
	double d_min = distanciaEuclidea(cent[0], yi, Uc, Uyi);
	int posDelMC = 0;
	int cardCent = cent.size();
	for (int j = 1; j < cardCent; ++j) {
		d_actual = distanciaEuclidea(cent[j], yi, Uc, Uyi);
		if (d_actual < d_min) {
			d_min = d_actual;
			posDelMC = j;
		}
	}
	return posDelMC;
}// metodo para saber cual es el centroide más cercano a yi (yi está en Y)

void KMeansP::pesoElemX2Y(vector<Elemento>& X, vector<Elemento>& Y, vector < vector<double> >& Ux, vector < vector<double> >& Uy) {
	int cardY = Y.size();
	int masCercano = 0;
	#	pragma omp parallel num_threads(thread_count) private(masCercano)
	{
		#	pragma omp for schedule(static)
		for (int i = 0; i < cardY; ++i) {
			masCercano = qEMCy_i(X, Y[i], Ux, Uy);
			X[masCercano].setPeso(X[masCercano].getPeso() + 1);
		}
	}
}

vector<double> KMeansP::centroide(vector <Elemento>  & conjY) { // revisar
	vector <double> centroide(n, 0.0);
	int cardY = conjY.size();
	double entrada_i = 0.0;
	if (cardY > 1 ) {
		for (int i = 0; i < n; ++i) {
			entrada_i = 0.0;
			for (int j = 0; j < cardY; ++j) {
				entrada_i += U[conjY[j].getId()][i];// U[conjY[j].getId()] esto me da el vector correspondiente al elemento j-esimo de J 
				// por otra parte, si V = U[conjY[j].getId()] es un vector entonces U[conjY[j].getId()][i] es equivalente a V[i]
			}
			centroide[i] = entrada_i / cardY; // esto llena la entrada i de centroide con el promedio de las entradas i de cada vector de J
		}
	}
	else {
		if (cardY > 0) {
			centroide = U[conjY[0].getId()];
		}
		else{
			srand(time(NULL));
			int r = (int)rand() % C.size();
			centroide=U[C[r].getId()];
		}
	}

	return centroide;
}

void KMeansP::setPesosAElementosDeC() {
	//cout << "asignando pesos a los elementos de C..." << endl;
	pesoElemX2Y(C, X, U, U); //pesos de los elementos de C con respecto al Conjunto X
	//inicialmente los element de C y X todos están en el universo U
	//cout << "pesos asignados" << endl;
}

void KMeansP::unir(vector<Elemento>& A, vector<Elemento>& B) { //hace A=AUB
	for (unsigned int i = 0; i < B.size(); ++i) {
		A.push_back(B[i]);
	}
}

void KMeansP::initC() { 
	srand(time(NULL));
	int r = (int) rand() % m;
	C.push_back(X[r]);
	X[r].setAdd2C(true);
	double P = 1.0;
	int iteraciones = 0;
	vector<Elemento> CC;
	Aleatorizador::inicializar_generador_random();
	while (P!=0 && iteraciones < 5) {
		//cout << "Calculando costo de it #" << iteraciones << " ..."<< endl;
		P = costoXRespectoY(X, C, U, U);
		//cout << "Costo calculado, ahora se hacen las probabilidades... " << endl;
		int x;
		#	pragma omp parallel num_threads(thread_count) shared(x,CC)
		{
			#	pragma omp for schedule(static)
			for (x = 0; x < m; ++x) {
				double aleatorio = Aleatorizador::random_uniform_real(Aleatorizador::generador);
				double dx2C = distanciaElementoConjunto(X[x], C, U, U);
				double probax = (this->l*dx2C*dx2C) / P;
				if (!X[x].getAdd2C() && probax >= aleatorio) {
					#	pragma omp critical
					{
						CC.push_back(X[x]);
					}
					X[x].setAdd2C(true);
				}
			}
		}
		//cout << "Haciendo una union+clear..." << endl;
		unir(C, CC); //C=C U CC
		CC.clear();
		//cout << "Finaliza init C iter #" << iteraciones << endl;
		++iteraciones;
	}
	//cout << "Candidatos listos, la cardinalidad de C es " << C.size() << endl;
}

void KMeansP::generarCentroidesDeC() {
	//cout << "generando los centroides para hacer agrupamientos en C..." << endl;
	srand(time(NULL));
	int r = (int)rand() % C.size();
	centroidesInicDeC.push_back(C[r]);
	int cardCentroides = centroidesInicDeC.size();
	int i = 0;
	double probac = 0;
	double distancia = 0;
	double aleatorio = 0;
	double costo = costoXRespectoY(C, centroidesInicDeC, U, U);
	while (cardCentroides < k) {
		distancia = distanciaElementoConjunto(C[i], centroidesInicDeC, U, U);
		costo = costoXRespectoY(C, centroidesInicDeC, U, U);
		probac = (distancia * distancia) / costo;
		aleatorio = Aleatorizador::random_uniform_real(Aleatorizador::generador);
		if (!C[i].getAdd2C() && probac <= aleatorio) {
			C[i].setAdd2C(true);
			centroidesInicDeC.push_back(C[i%C.size()]);
		}
		cardCentroides = centroidesInicDeC.size();
		i = (i + 1) % C.size();
	}
	centroidesUniverso.resize(k);
	centroides.resize(k);
	//cout << "centroides para hacer agrupamientos en C listos" << endl;
	//cout << "cantidad de centroides es " << cardCentroides << endl;
	for (int c = 0; c < k; ++c) {
		centroidesUniverso[c] = U[centroidesInicDeC[c].getId()]; // llenar el universo de centroides
		centroides[c].setGrupo(c);
		centroides[c].setId(c);

	}
}

void  KMeansP::KMeansSerial() { 
	//cout << "Ejecutando Kmeans ++ Serial sobre el conjunto <C> de candidatos a centroides..." << endl;
	int cardC = C.size();
	int grupo = 0;
	gruposKMeansS.resize(k); // k grupos de elementos 
	double costoAct = costoXRespectoY(C, centroides, U, centroidesUniverso);
	double costoAnt = 0;
	int it = 1;
	while (abs(costoAct - costoAnt) > EPSILONS) {
		costoAnt = costoAct;
		for (int e = 0; e < cardC; ++e) {
			grupo = centroideMC(centroides, C[e], centroidesUniverso, U);
			C[e].setGrupo(grupo);
			gruposKMeansS[grupo].push_back(C[e]);
		}
		for (int c = 0; c < k; ++c) {
			centroidesUniverso[c] = centroide(gruposKMeansS[c]); // generar un nuevo centroide para cada grupo
			gruposKMeansS[c].clear();
		}
		costoAct = costoXRespectoY(C, centroides, U, centroidesUniverso);
		++it;
	}
	for (int e = 0; e < cardC; ++e) {
		gruposKMeansS[C[e].getGrupo()].push_back(C[e]); // llenar los grupos con el resultado de la ultima (optima)
	}
	//cout << "Fin de KMeans ++, total de iteraciones: " << it << endl;
}

void KMeansP::ordenarGruposXpeso() {
	//cout << "Ordenando (x peso) cada grupo elaborado por KM ++ Serial... " << endl;
	int total_elem = 0;
	for (int grupo = 0; grupo < k; ++grupo) {
		total_elem += gruposKMeansS[grupo].size();
		sort(gruposKMeansS[grupo].begin(), gruposKMeansS[grupo].end());
	}
	//cout << "fin de ordenamiento..."  << endl;
}

void KMeansP::elegirCentroidesKMP() {
	//cout << "Eligiendo centroides de la 1 it de KMeans ||... " << endl;
	for (int c = 0; c < k; ++c) {
		if (gruposKMeansS[c].size() > 0) {
			centroidesUniverso[c] = U[gruposKMeansS[c][gruposKMeansS[c].size() - 1].getId()];
		}
		else {
			srand(time(NULL));
			int r = (int)rand() % X.size();
			centroidesUniverso[c]=U[X[r].getId()];
		}
	}
	//cout << "Centroides de la 1 it de KMeans || seleccionados" << endl;
}


void KMeansP::agrupar() {
	//cout << "ejecutando KMeans ||... " << endl;
	int cardX = X.size();
	int grupo = 0;
	gruposKMeansS.clear();
	gruposKMeansS.resize(k); // k grupos de elementos 
	double costoAct = costoXRespectoY(X, centroides, U, centroidesUniverso);
	double costoAnt = 0;
	int it = 1;
	int e;
	while (abs(costoAct - costoAnt) > epsilon_p) {
		costoAnt = costoAct;
		#	pragma omp parallel num_threads(thread_count)  private (grupo,e)
		{
			#	pragma omp for schedule(static)
			for (e = 0; e < cardX; ++e) {
				grupo = centroideMC(centroides, X[e], centroidesUniverso, U);
				X[e].setGrupo(grupo);
				#	pragma omp critical
				{
					gruposKMeansS[grupo].push_back(X[e]);
				}
			}
		}
		for (int c = 0; c < k; ++c) {
			centroidesUniverso[c] = centroide(gruposKMeansS[c]); // generar un nuevo centroide para cada grupo
			gruposKMeansS[c].clear();
		}
		costoAct = costoXRespectoY(C, centroides, U, centroidesUniverso);
		++it;
	}
	costoFinal=costoAct;
	for (int e = 0; e < cardX; ++e) {
		gruposKMeansS[X[e].getGrupo()].push_back(X[e]); // llenar los grupos con el resultado de la ultima (optima)
	}
	iteraciones_km=it;
	//cout << "Fin de KMeans ||, total de iteraciones: " << it << endl;
}