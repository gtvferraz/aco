#include "../include/Grafo.h"
#include <vector>
#include <time.h>
#include <math.h>
#include <random>
#include <ctime>
#include <algorithm>

using namespace std;

Grafo::Grafo(int numVertices)
{
    this->numVertices = 0;
    numArestas = 0;
    tamVetor = numVertices;
    vertices = new Vertice*[numVertices];
    for(int i=0; i<tamVetor; i++)
        vertices[i] = NULL;
}

Grafo::~Grafo()
{
    for(int i=0; i<tamVetor; i++)
        if(vertices[i] != NULL)
            delete vertices[i];
    delete [] vertices;
}

void Grafo::addAresta(int idPredecessor, int idSucessor, float peso)
{
    Vertice* predecessor = buscaVertice(idPredecessor);
    Vertice* sucessor = buscaVertice(idSucessor);

    if(predecessor == NULL)
        predecessor = addVertice(idPredecessor);
    if(sucessor == NULL)
        sucessor = addVertice(idSucessor);

    predecessor->addAresta(sucessor, peso);
    sucessor->addAresta(predecessor, peso);

    numArestas++;
}

bool Grafo::removeAresta(int idPredecessor, int idSucessor)
{
    Vertice* predecessor = buscaVertice(idPredecessor);
    Vertice* sucessor = buscaVertice(idSucessor);
    if(predecessor == NULL || sucessor == NULL)
        return false;

    Aresta* auxAresta = predecessor->getAdjacente();
    Aresta* arestaAnterior = NULL;
    while(auxAresta != NULL)
    {
        if(auxAresta->getSucessor()->getId() == idSucessor)
        {
            if(arestaAnterior != NULL)
                arestaAnterior->setProx(auxAresta->getProx());
            else
            {
                predecessor->setAdjacente(auxAresta->getProx());
            }
            auxAresta->setProx(NULL);
            delete auxAresta;
            numArestas--;
            return true;
        }
        arestaAnterior = auxAresta;
        auxAresta = auxAresta->getProx();
    }
    return false;
}

Vertice* Grafo::addVertice(int id)
{
    Vertice* auxVertice = vertices[id%tamVetor];
    numVertices++;
    if(auxVertice == NULL)
    {
        vertices[id%tamVetor] = new Vertice(id);
        return vertices[id%tamVetor];
    }
    while(auxVertice->getProx() != NULL)
        auxVertice = auxVertice->getProx();
    auxVertice->setProx(new Vertice(id));

    return auxVertice->getProx();
}

bool Grafo::removeVertice(int id)
{
    Vertice* auxVertice = vertices[id%tamVetor];
    Vertice* auxVertice2;
    Vertice* anterior = NULL;
    int aux = 0;
    while(auxVertice != NULL)
    {
        if(auxVertice->getId() == id)
        {
            for(int i=0; i<tamVetor; i++)
            {
                auxVertice2 = vertices[i];
                while(auxVertice2 != NULL)
                {
                    if(aux == numVertices)
                        break;
                    aux++;
                    removeAresta(auxVertice2->getId(), auxVertice->getId());
                    auxVertice2 = auxVertice2->getProx();
                }
            }
            if(anterior == NULL)
            {
                vertices[id%tamVetor] = auxVertice->getProx();
                auxVertice->setProx(NULL);
                delete auxVertice;
                numVertices--;
                return true;
            }
            anterior->setProx(auxVertice->getProx());
            numVertices--;
            auxVertice->setProx(NULL);
            delete auxVertice;
            return true;
        }
        anterior = auxVertice;
        auxVertice = auxVertice->getProx();
    }
    return false;
}

Vertice* Grafo::buscaVertice(int id)
{
    Vertice* auxVertice = vertices[id%tamVetor];
    while(auxVertice != NULL && auxVertice->getId() != id)
        auxVertice = auxVertice->getProx();
    return auxVertice;
}

int Grafo::getNumVertices()
{
    return numVertices;
}

Aresta* Grafo::buscaAresta(int idPredecessor, int idSucessor)
{
    Vertice* predecessor = buscaVertice(idPredecessor);
    Vertice* sucessor = buscaVertice(idSucessor);
    if(predecessor == NULL || sucessor == NULL)
        return NULL;

    Aresta* auxAresta = predecessor->getAdjacente();
    while(auxAresta != NULL)
    {
        if(auxAresta->getSucessor()->getId() == idSucessor)
            return auxAresta;
        auxAresta = auxAresta->getProx();
    }
    return NULL;
}

bool Grafo::taNoVetor(vector<int> vetor, int valor)
{
    for(int i=0; i<vetor.size(); i++)
        if(vetor.at(i) == valor)
            return true;
    return false;
}

void Grafo::apagaDoVetor(vector<int>vetor, int valor)
{
    int i;
    for(i=0; i<vetor.size(); i++)
        if(vetor.at(i) == valor)
            break;
    if(i < vetor.size())
        vetor.erase(vetor.begin()+i);
}

Vertice** Grafo::getVertices()
{
    return vertices;
}

void Grafo::imprimeArestas()
{
    int aux = 0;
    Vertice* auxVertice;
    for(int i=0; i<tamVetor; i++)
    {
        if(aux == numVertices)
            break;
        if(vertices[i] != NULL)
        {
            aux++;
            vertices[i]->imprimeSucessor();
            auxVertice = vertices[i]->getProx();
            cout << endl;
            while(auxVertice != NULL)
            {
                aux++;
                auxVertice->imprimeSucessor();
                auxVertice = auxVertice->getProx();
                cout << endl;
            }
        }
    }
}