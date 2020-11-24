#include "../include/Vertice.h"

using namespace std;

Vertice::Vertice(int id)
{
    this->id = id;
    adjacente = NULL;
    prox = NULL;
}

Vertice::~Vertice()
{
    if(adjacente != NULL)
        delete adjacente; ///DELETA A LISTA ENCADEADA DE ARESTAS DE SAÍDA DO VÉRTICE
    if(prox != NULL)
        delete prox; ///DELETA A LISTA ENCADEADA DE VÉRTICES
}

int Vertice::getId()
{
    return id;
}

Vertice* Vertice::getProx()
{
    return prox;
}

Aresta* Vertice::getAdjacente()
{
    return adjacente;
}

Aresta* Vertice::buscaMenorPeso()
{
    if(adjacente != NULL)
    {
        Aresta* auxAresta = adjacente;
        Aresta* retorno = adjacente;
        int menorPeso = adjacente->getPeso();
        while(auxAresta->getProx() != NULL)
        {
            auxAresta = auxAresta->getProx();
            if(auxAresta->getPeso() < menorPeso)
            {
                menorPeso = auxAresta->getPeso();
                retorno = auxAresta;
            }
        }
        return retorno;
    }
    return NULL;
}

Aresta* Vertice::addAresta(Vertice* sucessor, float peso)
{
    if(adjacente == NULL)
    {
        adjacente = new Aresta(sucessor, peso);
        return adjacente;
    }
    return adjacente->addProx(new Aresta(sucessor, peso));
}

void Vertice::setProx(Vertice* prox)
{
    this->prox = prox;
}

void Vertice::addAresta(Aresta* aresta)
{
    if(adjacente == NULL)
        adjacente = aresta;
    else
        adjacente->addProx(aresta);
}

void Vertice::setAdjacente(Aresta* aresta)
{
    adjacente = aresta;
}

void Vertice::imprimeSucessor()
{
    cout << id;
    if(adjacente != NULL)
        adjacente->imprimeSucessor();
}
