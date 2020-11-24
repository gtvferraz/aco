#ifndef ARESTA_H
#define ARESTA_H

#include <cstddef>
#include <iostream>

#include "Vertice.h"

class Aresta
{
    public:
        Aresta(Vertice* sucessor, double peso);
        ~Aresta();
        int getId();
        double getPeso();
        Vertice* getSucessor();
        Aresta* getProx(); ///RETORNA A PRÓXIMA ARESTA DA LISTA ENCADEADA
        Aresta* addProx(Aresta *prox); ///ADICIONA UMA ARESTA NA LISTA ENCADEADA
        void setProx(Aresta *prox);
        void setPeso(double peso);
        void imprimeSucessor();

    private:
        double peso;
        Vertice* sucessor;
        Aresta* prox; ///LISTA ENCADEADA DE ARESTAS
};

#endif // ARESTA_H
