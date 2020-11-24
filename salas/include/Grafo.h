#ifndef GRAFO_H
#define GRAFO_H

#include <cstddef>
#include <iostream>
#include <vector>

#include "Aresta.h"

using namespace std;

class Grafo
{
    public:
        Grafo(int numVertices);
        ~Grafo();
        int getNumVertices();
        Vertice** getVertices(); ///RETORNA O VETOR DE ARMAZENAMENTO DE VÉRTICES
        Vertice* addVertice(int id); ///ADICIONA UM VÉRTICE
        Vertice* buscaVertice(int id); ///RETORNA O VÉRTICE COM O ID PASSADO POR PARÂMETRO
        void addAresta(int idPredecessor, int idSucessor, float peso); ///ADICIONA UMA ARESTA
        void imprimeArestas(); ///IMPRIME O GRAFO
        bool removeVertice(int id); ///REMOVE UM VÉRTICE
        bool removeAresta(int idPredecessor, int idSucessor); ///REMOVE UMA ARESTA
        Aresta* buscaAresta(int idPredecessor, int idSucessor); ///RETORNA A ARESTA IDENTIFICADA PELOS VÉRTICES PASSADOS POR PARÂMETROS

    private:
        int numVertices; ///NÚMERO DE VÉRTICES DO GRAFO
        int numArestas; ///NÚMERO DE ARESTAS DO GRAFO
        int tamVetor; ///TAMANHO FIXO PARA O VETOR DE ARMAZENAMENTO DE VÉRTICES
        Vertice** vertices; ///VETOR DE ARMAZENAMENTO DE VÉRTICES. CADA ELEMENTO DO VETOR É UMA LISTA ENCADEADA DE VÉRTICES
        bool taNoVetor(vector<int> vetor, int valor); ///VERIFICA SE O INTEIRO VALOR ESTÁ PRESENTE NO VECTOR VETOR
        void apagaDoVetor(vector<int> vetor, int valor); ///REMOVE O ELEMENTO VALOR DO VECTOR VETOR
};

#endif // GRAFO_H
