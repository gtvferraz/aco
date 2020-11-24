#include "../include/Aresta.h"

using namespace std;

///INICIA A ARESTA COM O VÉRTICE PARA O QUAL ELA APONTA, O PESO DA ARESTA E SUA ID
Aresta::Aresta(Vertice* sucessor, double peso)
{
    this->sucessor = sucessor;
    this->peso = peso;
    prox = NULL;
}

Aresta::~Aresta()
{
    if(prox != NULL)
        delete prox; /// DELETA TODA A LISTA ENCADEADA DE ARESTAS
}

double Aresta::getPeso()
{
    return peso;
}

Vertice* Aresta::getSucessor()
{
    return sucessor;
}

Aresta* Aresta::getProx()
{
   return prox;
}

Aresta* Aresta::addProx(Aresta* prox)
{
    if(this->prox == NULL)
    {
        this->prox = prox;
        return prox;
    }
    return this->prox->addProx(prox);
}

void Aresta::setProx(Aresta *prox)
{
    this->prox = prox;
}

void Aresta::setPeso(double peso)
{
    this->peso = peso;
}

void Aresta::imprimeSucessor()
{
    cout << " - ";
    cout << sucessor->getId();
    if(prox != NULL)
            prox->imprimeSucessor();
}
