#pragma warning( disable : 4290 )
#pragma warning( disable : 4290 )
#pragma warning( disable : 5040 )
#include "KMeansP.h"
#include "Elemento.h"
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <exception>
#include <iostream>
#include <omp.h>
using namespace std;
// /Zc:twoPhase- 
//COMPILAR DESDE LINEA DE COMANDOS: g++ -g -Wall -fopenmp -o km *.cpp
double stod_wrapper(string v) throw (invalid_argument, out_of_range) { return std::stod(v); }
template < typename T, class F >
vector< vector< T > > carga_valida_datos(ifstream& archivo, F t) throw (invalid_argument, out_of_range);
vector< vector< double > > recuperarArch();
void recuperarNumeros(double & l, double & epsilon); 


int main(int argc, const char * argv[]) {

	int cp = omp_get_num_procs();
	int hilos = cp*3;
	double epsilon = 0.0;
	double l = 0.0;
	vector < vector < double  > > datos=recuperarArch(); 
	recuperarNumeros(l,epsilon);
	int m = datos.size();
	int n= datos[0].size();
	//parametros orden: n,m,k,l,datos,hilos
	int k = 100;
	KMeansP algoritmo(n,m,k,l,datos,hilos,epsilon);
	cout << "***INFO DEL PROGRAMA***" << endl;
	algoritmo.verInfo();
	
	double inicioG= omp_get_wtime();

	algoritmo.initC(); //inicializa conjunto de candidatos a centroides
	algoritmo.setPesosAElementosDeC(); // asigna el peso a tales candidatos
	algoritmo.generarCentroidesDeC(); // genera centroides iniciales para hacer kmeans ++ sobre el conjunto de candidatos
	algoritmo.KMeansSerial(); //ejecuta kmeans++ sobre el conjunto de candidatos
	algoritmo.ordenarGruposXpeso(); // ordena por peso los elementos de cada grupo creado en el paso anterior
	algoritmo.elegirCentroidesKMP(); // elige los centroides de KMEANS || (elementos + pesados) de cada grupo anterior
	algoritmo.agrupar(); // EJECUTA KMEANS || sobre TODO el conjunto de datos ingresados 
	
	double finG = omp_get_wtime();
	cout << "Se ha ejecutado KMeans || las iteraciones de KMeans || para convergencia son: " << algoritmo.getIterKM() << endl;
	cout << "Tiempo pared: " << finG - inicioG << " segundos " << endl;
	cout << "Tiempo pared: " << (finG - inicioG)/60 << " minutos " << endl;
	algoritmo.setDuracion(finG - inicioG);
	cout << "Generando archivo..." << endl;
	algoritmo.generarArchivo();
	cout << "Archivo generado" << endl;
	
	int salir;
	cout << endl << endl << "digite cualquier tecla + enter para salir... ";
	cin>>salir;
	return 0;
}


template < typename T, class F >
vector< vector< T > > carga_valida_datos(ifstream& archivo, F t) throw (invalid_argument, out_of_range)
{
	vector< vector< T > > valores;
	vector< T > linea_valores;
	string linea;
	while (getline(archivo, linea)) {
		linea_valores.clear();
		stringstream ss(linea);
		string numero_S;
		T numero_T;
		while (getline(ss, numero_S, ',')) {
			try {
				numero_T = t(numero_S);
			}
			catch (exception e) {
				throw e;
			}
			linea_valores.push_back(numero_T);
		}
		valores.push_back(linea_valores);
	}
	return valores;
}

void recuperarNumeros(double & l, double & epsilon) {
	while (l <= 0 || epsilon <= 0) {
		cout << "digite el valor de epsilon (se recomienda 100) ";
		cin >> epsilon;
		cout << endl << "digite el valor de l (se recomienda 0.3) ";
		cin >> l;
		cout << endl;
	}
}

vector< vector< double > > recuperarArch(){
	cout << "El programa por default (como dicen las instrucciones) usa k=100 si desea cambiar <k> modifique en linea 34 de Main.cpp" << endl;
	bool capturado=true;
	string file_name="";
	cout << "digite el nombre del archivo (ejemplo: <vectores_desordenados.csv>) :";
	cin >> file_name;
	cout << endl;
	ifstream d(file_name, ios::in);
	if (!d){
		cout << "no encuentra el archivo de datos...Vuelva a intentar" << endl;
		capturado=false;
	}
	while(!capturado){
		cout << "digite el nombre del archivo (ejemplo: <vectores_desordenados.csv>) :";
		cin >> file_name;
		cout << endl;
		ifstream d(file_name);
		if (!d){
			cout << "no encuentra el archivo de datos...Vuelva a intentar" << endl;
			capturado=false;
		}
	}
	
	vector< vector< double > > vd;
	try {
		vd = carga_valida_datos< double >(d, stod_wrapper); 
	}
	catch (exception e) {
		cout << "valor invalido o fuera de limite" << endl;
	}
	return vd;
}