#ifndef VERTICE_H
#define VERTICE_H

#include <iostream>
#include <cstddef>

class Vertice;
class Aresta;

#include "Aresta.h"

class Vertice
{
    public:
        Vertice(int id);
        ~Vertice();
        int getId();
        Vertice* getProx();
        Aresta* getAdjacente();
        Aresta* buscaMenorPeso(); ///RETORNA A ARESTA DE SAÍDA DO VÉRTICE COM O MENOR PESO
        Aresta* addAresta(Vertice* sucessor, float peso); ///ADICIONA UMA ARESTA DE SAÍDA NO VÉRTICE
        void setProx(Vertice* prox);
        void addAresta(Aresta* aresta); ///ADICIONA UMA ARESTA DE SAÍDA CRIADA ANTERIORMENTE NO VÉRTICE
        void setAdjacente(Aresta* aresta);
        void imprimeSucessor();

    private:
        int id;
        Aresta* adjacente; ///LISTA ENCADEADA DAS ARESTAS DE SAÍDA DO VÉRTICE
        Vertice* prox; ///LISTA ENCADEADA DE VÉRTICES CUJO RESULTADO DA FUNÇÃO DE HASHING DEU O MESMO (COLISÃO)
};

#endif // VERTICE_H
