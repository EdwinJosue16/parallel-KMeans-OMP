#include "Elemento.h"

Elemento::Elemento() {
	id = 0;//
	peso = 0;
	mi_grupo = 0;//
	add2C = false;
}

Elemento::Elemento(int id, int peso, int grupo) {
	this->id = id;
	this->peso = peso;
	this->mi_grupo = grupo;
	add2C = false;
}

Elemento::Elemento(const Elemento& otro) {
	setId(otro.id);
	setPeso(otro.peso);
	setGrupo(otro.mi_grupo);
	setAdd2C(otro.add2C);
}

void Elemento::setId(int id) {
	this->id = id;
}

void Elemento::setPeso(int peso) {
	this->peso = peso;
}

void Elemento::setGrupo(int grupo) {
	this->mi_grupo = grupo;
}

void Elemento::setAdd2C(bool add) {
	this->add2C = add;
}

int Elemento::getPeso() {
	return peso;
}

int Elemento::getId() {
	return   id;
}
int Elemento::getGrupo() {
	return mi_grupo;
}

bool Elemento::getAdd2C() {
	return add2C;
}

bool Elemento::operator == ( const Elemento & otro) {
	return this->peso == otro.peso;
}

bool Elemento::operator != (const Elemento & otro) {
	return this->peso != otro.peso;
}

bool Elemento::operator < (const Elemento & otro) {
	return this->peso < otro.peso;
}

bool Elemento::operator > (const Elemento & otro) {
	return this->peso > otro.peso;
}

bool Elemento::operator <= (const Elemento & otro) {
	return this->peso < otro.peso || this->peso == otro.peso;
}

bool Elemento::operator >= (const Elemento & otro) {
	return this->peso > otro.peso || this->peso == otro.peso;
}
Elemento& Elemento::operator = (const Elemento & otro) {
	setId(otro.id);
	setPeso(otro.peso);
	setGrupo(otro.mi_grupo);
	setAdd2C(otro.add2C);
	return *this;
}