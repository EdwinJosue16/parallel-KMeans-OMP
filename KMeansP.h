#pragma once 
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <utility>
#include "Elemento.h"
#include <fstream>
using namespace std;
class KMeansP {
	private:
		int n; //dimension de los vectores que se van a agrupara cada v pert a R^n
		int m; // cantidad de vectores
		int k; // cantidad de clusters o grupos que se desean formar
		int thread_count; // cantidad de hilos (debe ser pasado en el main)
		double l; // constante usada en la init de C (l=0.5)
		double epsilon_p;
		double costoFinal;
		double tiempoPared;
		int iteraciones_km;
		vector < vector<double> > U; // este es el espacio de cardinalidad m, que contiene todos los elementos que se desean agrupar, cada elemento es un vector de dim n
		vector< vector <double> > centroidesUniverso; // es un vector que contiene los centroides que se usan dentro de KMS y KMP, seria una Universo de centroides y es de cardinalidad k
		vector < Elemento > X; // este es el conjunto que contiene los indices cada elementos (indice y peso de cada vector en U)
		vector < Elemento > C; // es el conjunto de centroides candidatos (almacena el indice de cada vector cantidato y su respectivo peso)
		vector< Elemento > centroidesInicDeC; // centroides manejados como elementos
		vector < vector <Elemento> > gruposKMeansS;
		vector< Elemento > centroides; // centroides manejados como elementos
		void llenarX();
		double distanciaElementoConjunto(Elemento &, vector<Elemento> &, vector < vector <double> >&, vector < vector <double> >&); // metodo para calcular la distancia de cada x en X a C
		double costoXRespectoY(vector <Elemento> &, vector<Elemento> &, vector< vector <double> >&, vector< vector <double> >&);
		double distanciaEuclidea(Elemento &, Elemento &, vector< vector <double> >&, vector< vector <double> >&); // los vectores son los universos a los cuales pertenece cada elemento
		vector<double> centroide(vector <Elemento> &);
		int qEMCy_i(vector<Elemento>&, Elemento&, vector < vector<double> >&, vector < vector<double> >&); // metodo para saber cual es el elemento en X más cercano a yi (yi está en Y)
		int centroideMC(vector<Elemento>&, Elemento&, vector < vector<double> >&, vector < vector<double> >&);
		void pesoElemX2Y(vector<Elemento>&, vector<Elemento>&, vector < vector<double> >&, vector < vector<double> >&);
		int pos_del_minimo(vector<double>&);
		void unir(vector<Elemento>&, vector<Elemento>&);
public:
		KMeansP();
		KMeansP(int,int,int,double, vector < vector<double> >&,int,double);
		void initC(); // es el algoritmo que se ejecuta logLETRAG veces
		void setPesosAElementosDeC(); // For each c in C, set Wc to be the number of points in X closer to c than any other point in C
		void generarCentroidesDeC(); // llenar a centroidesInicDeC
		void  KMeansSerial(); // una vez C init listo y los centroidesInicDeC listos, se hace KmeansSerial dentro de C 
		// KMeansSerial llena al vector de grupos hechos en C (gruposKMeansS), donde cada grupo está formado por los indices de los elementos agrupados
		void ordenarGruposXpeso(); //ordenar cada grupo de gruposKMeansS según su peso
		void elegirCentroidesKMP(); // seleccionar a los elementos de mayor peso en cada grupo de gruposKMeansS 
		//este algoritmo, entonces, va llenar a centroidesInicDeC
		void agrupar(); //KMeansP
		
		int getIterKM(){
			return iteraciones_km;
		}
		
		void verInfo(){
			cout << "el valor de k*l es " << this->l << endl;
			cout << "el valor de epsilon es " << this->epsilon_p << endl;
			cout << "la cantidad de clusters k es " << this->k << endl;
			cout << "la dimension de los vectores es " << this->n<< endl;
			cout << "la cantidad de elementos es " << this->m << endl;
			cout << "la cantidad de hilos es " << this->thread_count << endl;
		}
		
		void setDuracion(double tp){
			this->tiempoPared=tp;
		}
		
		string vector2string(vector<double> & vec){
			string hilera="";
			for(unsigned int i=0;i<vec.size();++i){
				hilera+=to_string(vec[i]);
				hilera.append(",");
			}
			hilera.erase(hilera.size()-1);
			return hilera;
		}
		
		void verResultados() {

			cout << endl << endl << endl;
			for (int grupo = 0; grupo < k; ++grupo) {
				cout << "el grupo #" << grupo << "  tiene "<< gruposKMeansS[grupo].size() << " elementos" << endl;
			}

		}
		
		void generarArchivo(){
			ofstream salida("salida.csv");
			for (int grupo = 0; grupo < k; ++grupo) {
				salida << vector2string(centroidesUniverso[grupo]) << endl;
				for(unsigned int e=0;e<gruposKMeansS[grupo].size();++e){
					salida << vector2string(U[gruposKMeansS[grupo][e].getId()]) << endl;
				}
				salida << to_string(gruposKMeansS[grupo].size()) << endl;
			}
			salida << "Tiempo pared: " << tiempoPared << " s -> " << "Costo: " << costoFinal;
		}

};
