#pragma once
#include <iostream>
using namespace std;
class Elemento {
	private:
		int peso;
		int id;
		int mi_grupo;
		bool add2C; // es un booleano para no agregar repetidos en la init de C
	public:
		Elemento();
		Elemento(int, int, int);
		Elemento(const Elemento&);
		void setPeso(int);
		void setId(int);
		void setGrupo(int);
		void setAdd2C(bool);
		int getPeso();
		int getId();
		int getGrupo();
		bool getAdd2C();
		bool operator == (const Elemento &);
		bool operator != (const Elemento &);
		bool operator < (const Elemento &);
		bool operator > (const Elemento &);
		bool operator <= (const Elemento &);
		bool operator >= (const Elemento &);
		Elemento& operator = (const Elemento &);
};