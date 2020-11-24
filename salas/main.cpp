#include <iostream>
#include <algorithm> 
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <time.h>
#include "include/Grafo.h"

//Comando para compilar o programa
//g++ include/*.h src/*.cpp main.cpp -O3
//Comando para executar
//./a.out

#define NUMALPHAS 10

using namespace std;

typedef struct {
    int id;
    int numSalas;
    int largura;
    int *fluxo = nullptr;
} Sala;

struct AuxOrdena{
    int idCandidato;
    double fluxoCandidato;

    AuxOrdena(int idCandidato, double fluxoCandidato) {
        this->idCandidato = idCandidato;
        this->fluxoCandidato = fluxoCandidato;
    }
};

struct AuxOrdenaSolucao{
    double x;
    Sala* sala;

    AuxOrdenaSolucao(double x, Sala* sala) {
        this->x = x;
        this->sala = sala;
    }
};

vector<Sala>* carregaInstancia(string nomeArquivo) {
    vector<Sala> *salas;
    int numSalas;
    int *fluxo;
    ifstream arquivoEntrada;
    string str;
    stringstream ss;
    arquivoEntrada.open(nomeArquivo);


    if(arquivoEntrada.is_open()) {
        salas = new vector<Sala>;

        getline(arquivoEntrada, str);
        ss << str;
        ss >> numSalas;
        //cout << numSalas;
        getline(arquivoEntrada, str);
        ss.clear();
        ss << str;
        for(int i=0; i<numSalas; i++) {
            getline(ss, str, ' '); //Espaço para instâncias do artigo
            Sala aux;
            aux.id = i;
            aux.numSalas = numSalas;
            aux.largura = stoi(str);
            salas->push_back(aux);
        }

        for(int i=0; i<numSalas; i++) {
            getline(arquivoEntrada, str);
            ss.clear();
            ss << str;
            fluxo = new int[numSalas];
            for(int j=0; j<numSalas; j++) {
                getline(ss, str, ' ');
                fluxo[j] = stoi(str);
            }
            salas->at(i).fluxo = fluxo;
        }
    } else {
        cout << "Arquivo nao Encontrado" << endl;
        return nullptr;
    }
    return salas;
}

bool compara_sort(AuxOrdena* a, AuxOrdena* b) {
    return (a->fluxoCandidato > b->fluxoCandidato);
}

bool compara_sort_solucao(AuxOrdenaSolucao* a, AuxOrdenaSolucao* b) {
    return (a->x < b->x);
}

bool verificaSolucao(vector<Sala> *salas, double *solucao) {
    int numSalas = salas->size();
    int indiceSup;
    bool valida;
    double posicaoAtual;
    vector<AuxOrdenaSolucao*> auxSolucao;
    
    for(int i=0; i<numSalas; i++) 
        auxSolucao.push_back(new AuxOrdenaSolucao(solucao[i], &salas->at(i)));
    sort(auxSolucao.begin(), auxSolucao.end(), compara_sort_solucao);

    for(int i=0; i<numSalas; i++)
        if(auxSolucao[i]->x > 0) {
            indiceSup = i;
            break;
        }

    valida = true;

    //CORREDOR INFERIOR
    posicaoAtual = 0;
    for(int i=indiceSup-1; i>=0; i--) {
        if(posicaoAtual - (auxSolucao[i]->sala->largura / 2.0) != auxSolucao[i]->x) {
            valida = false;
            break;
        }
        posicaoAtual -= auxSolucao[i]->sala->largura;
    }

    //CORREDOR SUPERIOR
    if(valida) {
        posicaoAtual = 0;
        for(int i=indiceSup; i<numSalas; i++) {
            if(posicaoAtual + (auxSolucao[i]->sala->largura / 2.0) != auxSolucao[i]->x) {
                valida = false;
                break;
            }
            posicaoAtual += auxSolucao[i]->sala->largura;
        }
    }

    for(int i=0; i<numSalas; i++) 
        delete auxSolucao[i];

    return valida;
}

double calculaCusto(vector<Sala> *salas, double *solucao) {
    double custo = 0;

    if(solucao == nullptr)
        return -1;

    for(int i=0; i<salas->size(); i++) {
        for(int j=i+1; j<salas->size(); j++)
                custo += fabs(fabs(solucao[i]) - fabs(solucao[j]))*salas->at(i).fluxo[j];
    }
    return custo;
}

double estimaCusto(vector<Sala> *salas, double custo, double *solucao, double posCandidato, int idCandidato) {
    for(int i=0; i<salas->at(0).numSalas; i++) {
        if(solucao[i] != 0)
            custo += fabs(fabs(solucao[i]) - fabs(posCandidato))*salas->at(i).fluxo[idCandidato];
    }
    return custo;
}

void apresentaSolucao(vector<Sala> *salas, double* solucao) {
    cout << "Custo: " << calculaCusto(salas, solucao);
    if(verificaSolucao(salas, solucao))
        cout << ", Valida" << endl;
    else
        cout << ", Invalida" << endl;
    cout << endl;
}

int selecionaOrigem(vector<Sala> *salas) {
    int numSalas = salas->size();
    int idMenorFluxo;
    double somatorioFluxos;
    double menorFluxo = INFINITY;

    for(int i=0; i<numSalas; i++) {
        somatorioFluxos = 0;
        for(int j=0; j<numSalas; j++)
            somatorioFluxos += salas->at(i).fluxo[j];
        if(somatorioFluxos < menorFluxo) {
            menorFluxo = somatorioFluxos;
            idMenorFluxo = i;
        }
    }
    return idMenorFluxo;
}

void buscaLocalSwap(vector<Sala> *salas, double* solucao) {
    int numSalas = salas->size();
    int indiceSup;
    bool atualizou;
    double diferenca;
    double custo;
    double novoCusto;
    vector<AuxOrdenaSolucao*> novaSolucao;
    vector<vector<double>*> fluxos;
    Sala* auxSwap;

    custo = calculaCusto(salas, solucao);

    for(int i=0; i<numSalas; i++) {
        novaSolucao.push_back(new AuxOrdenaSolucao(solucao[i], &salas->at(i)));
        fluxos.push_back(new vector<double>());
        for(int j=0; j<numSalas; j++) {
            fluxos[i]->push_back(salas->at(i).fluxo[j]);
        }
    }

    atualizou = true;
    while(atualizou) {
        atualizou = false;
        sort(novaSolucao.begin(), novaSolucao.end(), compara_sort_solucao);

        for(int i=0; i<numSalas; i++)
            if(novaSolucao[i]->x > 0) {
                indiceSup = i;
                break;
            }
        
        //SWAP INFERIOR COM INFERIOR
        for(int i=0; i<indiceSup; i++) {
            for(int j=i+1; j<indiceSup; j++) {
                diferenca = novaSolucao[i]->sala->largura - novaSolucao[j]->sala->largura;
                auxSwap = novaSolucao[i]->sala;
                novaSolucao[i]->sala = novaSolucao[j]->sala;
                novaSolucao[j]->sala = auxSwap;
            
                novaSolucao[j]->x -= diferenca/2;
                for(int k=j-1; k>i; k--) 
                    novaSolucao[k]->x -= diferenca;
                novaSolucao[i]->x -= diferenca/2;

                novoCusto = 0;
                for(int k=0; k<novaSolucao.size(); k++) {
                    for(int l=k+1; l<novaSolucao.size(); l++)
                            novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                }
                
                if(novoCusto < custo) {
                    custo = novoCusto;
                    atualizou = true;
                    break;
                } else {
                    auxSwap = novaSolucao[i]->sala;
                    novaSolucao[i]->sala = novaSolucao[j]->sala;
                    novaSolucao[j]->sala = auxSwap;

                    diferenca *= -1;
                    novaSolucao[j]->x -= diferenca/2;
                    for(int k=j-1; k>i; k--) 
                        novaSolucao[k]->x -= diferenca;
                    novaSolucao[i]->x -= diferenca/2;

                    novoCusto = 0;
                    for(int k=0; k<novaSolucao.size(); k++) {
                        for(int l=k+1; l<novaSolucao.size(); l++)
                                novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                    }
                }
            }
            if(atualizou)
                break;
        }
        
        //SWAP SUPERIOR COM SUPERIOR
        for(int i=indiceSup; i<numSalas; i++) { 
            for(int j=i+1; j<numSalas; j++) {
                diferenca = novaSolucao[j]->sala->largura - novaSolucao[i]->sala->largura;
                auxSwap = novaSolucao[i]->sala;
                novaSolucao[i]->sala = novaSolucao[j]->sala;
                novaSolucao[j]->sala = auxSwap;
            
                novaSolucao[i]->x += diferenca/2;
                for(int k=i+1; k<j; k++)
                    novaSolucao[k]->x += diferenca;
                novaSolucao[j]->x += diferenca/2;

                novoCusto = 0;
                for(int k=0; k<novaSolucao.size(); k++) {
                    for(int l=k+1; l<novaSolucao.size(); l++)
                        novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                }

                if(novoCusto < custo) {
                    custo = novoCusto;
                    atualizou = true;
                    break;
                } else {
                    auxSwap = novaSolucao[i]->sala;
                    novaSolucao[i]->sala = novaSolucao[j]->sala;
                    novaSolucao[j]->sala = auxSwap;
                
                    diferenca *= -1;
                    novaSolucao[i]->x += diferenca/2;
                    for(int k=i+1; k<j; k++)
                        novaSolucao[k]->x += diferenca;
                    novaSolucao[j]->x += diferenca/2;

                    novoCusto = 0;
                    for(int k=0; k<novaSolucao.size(); k++) {
                        for(int l=k+1; l<novaSolucao.size(); l++)
                            novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                    }
                }
            }
            if(atualizou)
                break;
        }
        
        //SWAP INFERIOR COM SUPERIOR
        for(int i=0; i<indiceSup; i++) {
            for(int j=indiceSup; j<numSalas; j++) {
                diferenca = novaSolucao[i]->sala->largura - novaSolucao[j]->sala->largura;
                auxSwap = novaSolucao[i]->sala;
                novaSolucao[i]->sala = novaSolucao[j]->sala;
                novaSolucao[j]->sala = auxSwap;
            
                novaSolucao[i]->x += diferenca/2;
                novaSolucao[j]->x += diferenca/2;
                for(int k=i-1; k>=0; k--) 
                    novaSolucao[k]->x += diferenca;
                for(int k=j+1; k<numSalas; k++) 
                    novaSolucao[k]->x += diferenca;


                novoCusto = 0;
                for(int k=0; k<novaSolucao.size(); k++) {
                    for(int l=k+1; l<novaSolucao.size(); l++)
                            novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                }
                
                if(novoCusto < custo) {
                    custo = novoCusto;
                    atualizou = true;
                    break;
                } else {
                    auxSwap = novaSolucao[i]->sala;
                    novaSolucao[i]->sala = novaSolucao[j]->sala;
                    novaSolucao[j]->sala = auxSwap;
                
                    diferenca *= -1;
                    novaSolucao[i]->x += diferenca/2;
                    novaSolucao[j]->x += diferenca/2;
                    for(int k=i-1; k>=0; k--) 
                        novaSolucao[k]->x += diferenca;
                    for(int k=j+1; k<numSalas; k++) 
                        novaSolucao[k]->x += diferenca;

                    novoCusto = 0;
                    for(int k=0; k<novaSolucao.size(); k++) {
                        for(int l=k+1; l<novaSolucao.size(); l++)
                                novoCusto += fabs(fabs(novaSolucao[k]->x) - fabs(novaSolucao[l]->x))*fluxos[novaSolucao[k]->sala->id]->at(novaSolucao[l]->sala->id);
                    }
                }
            }
            if(atualizou)
                break;
        }
    }

    for(int i=0; i<numSalas; i++)
        solucao[novaSolucao[i]->sala->id] = novaSolucao[i]->x;

    for(int i=0; i<numSalas; i++) {
        delete fluxos[i];
        delete novaSolucao[i];
    }
}

double* guloso(vector<Sala> *salas) {
    int numSalas = salas->size();
    int numSalasInseridas = 0;
    int idAtual;
    int idMaiorFluxo;
    int idUltimaAdd;
    double maiorFluxo;
    double *solucao = new double[numSalas];
    double espOcupadoSup = 0;
    double espOcupadoInf = 0;
    double espOcupadoMenor = 0;
    double larguraMaiorFluxo;
    double custo = 0;
    double custoAux;
    bool *salasInseridas = new bool[numSalas];
    int *ordemInseridas = new int[numSalas];

    for(int i=0; i<numSalas; i++) {
        solucao[i] = 0;
        salasInseridas[i] = false;
    }

    //idAtual = selecionaOrigem(salas);
    idAtual = rand()%numSalas;
    espOcupadoSup = salas->at(idAtual).largura;
    solucao[idAtual] = espOcupadoSup/2;
    salasInseridas[idAtual] = true;
    numSalasInseridas++;
    ordemInseridas[0] = idAtual;

    idUltimaAdd = idAtual;
    while(numSalasInseridas < numSalas) {
        maiorFluxo = -1;
        for(int i=0; i<numSalas; i++) {
            if(!salasInseridas[i]) {
                //custoAux = fabs(fabs(solucao[idUltimaAdd]) - (espOcupadoMenor + salas->at(i).largura/2.0)) * salas->at(idUltimaAdd).fluxo[i];
                
                custoAux = 0;
                for(int j = 0; j < numSalasInseridas; j++){
                    custoAux += (1.0/pow(2,j))*fabs(fabs(solucao[ordemInseridas[j]]) - (espOcupadoMenor + salas->at(i).largura/2.0)) * salas->at(ordemInseridas[j]).fluxo[i];
                }
                if(custoAux > maiorFluxo) {
                    idMaiorFluxo = i;
                    maiorFluxo = custoAux;
                    larguraMaiorFluxo = salas->at(i).largura;
                }
            }
        }
        
        if(espOcupadoSup < espOcupadoInf) {
            solucao[idMaiorFluxo] = espOcupadoSup + larguraMaiorFluxo/2.0;
            espOcupadoSup += larguraMaiorFluxo;
        } else {
            solucao[idMaiorFluxo] = (espOcupadoInf + larguraMaiorFluxo/2.0)*-1;
            espOcupadoInf += larguraMaiorFluxo;
        }
        if(espOcupadoSup < espOcupadoInf)
            espOcupadoMenor = espOcupadoSup;
        else
            espOcupadoMenor = espOcupadoInf;
        
        salasInseridas[idMaiorFluxo] = true;
        numSalasInseridas++;

        idUltimaAdd = idMaiorFluxo;
        ordemInseridas[numSalasInseridas-1] = idMaiorFluxo;
    }

    delete[] salasInseridas;
    return solucao;
}

double* auxGulosoRandomizado(vector<Sala> *salas, float alpha, int seed) {
    int numSalas = salas->size();
    int aleatorio;
    int idSala;
    int idUltimaAdd;
    double *solucao = new double[numSalas];
    double espOcupadoSup = 0;
    double espOcupadoInf = 0;
    double espOcupadoMenor = 0;
    vector<AuxOrdena*> candidatos;
    int *ordemInseridas = new int[numSalas];

    for(int i=0; i<numSalas; i++) {
        solucao[i] = 0;
        candidatos.push_back(new AuxOrdena(i, 0));
    }

    idSala = rand()%numSalas;
    espOcupadoSup = salas->at(idSala).largura;
    solucao[idSala] = espOcupadoSup/2.0;
    delete candidatos[idSala];
    candidatos.erase(candidatos.begin()+idSala);
    ordemInseridas[0] = idSala;
    idUltimaAdd = idSala;
    while(candidatos.size() > 0) {
        for(int i=0; i<candidatos.size(); i++) {
            candidatos[i]->fluxoCandidato = (fabs(fabs(solucao[idUltimaAdd]) - (espOcupadoMenor + salas->at(candidatos[i]->idCandidato).largura/2.0)) * salas->at(idUltimaAdd).fluxo[candidatos[i]->idCandidato]) / salas->at(candidatos[i]->idCandidato).largura;
        }
        sort(candidatos.begin(), candidatos.end(), compara_sort);

        aleatorio = rand()%(int)ceil(alpha*candidatos.size());
        idSala = candidatos[aleatorio]->idCandidato;
        
        if(espOcupadoSup < espOcupadoInf) {
            solucao[idSala] = espOcupadoSup + salas->at(idSala).largura/2.0;
            espOcupadoSup += salas->at(idSala).largura;
        } else {
            solucao[idSala] = (espOcupadoInf + salas->at(idSala).largura/2.0)*-1.0;
            espOcupadoInf += salas->at(idSala).largura;
        }
        if(espOcupadoSup < espOcupadoInf)
            espOcupadoMenor = espOcupadoSup;
        else
            espOcupadoMenor = espOcupadoInf;

        ordemInseridas[numSalas-candidatos.size()] = idSala;

        delete candidatos[aleatorio];
        candidatos.erase(candidatos.begin()+aleatorio);
    }
    
    return solucao;
}

double* gulosoRandomizado(vector<Sala> *salas, float alpha, int seed) {
    double menorCusto = INFINITY;
    double custo;
    double *melhorSolucao = nullptr;
    double *solucao;
    for(int i=0; i<100; i++) {
        solucao = auxGulosoRandomizado(salas, alpha, seed+i);
        custo = calculaCusto(salas, solucao);
        if(custo < menorCusto) {
            if(melhorSolucao != nullptr)
                delete[] melhorSolucao;
            melhorSolucao = solucao;
            menorCusto = custo;
        } else
            delete[] solucao;
    }
    return melhorSolucao;
}

double* gulosoReativo(vector<Sala> *salas, int numInteracoes, int tamanhoBloco, int seed) {
    vector<double> alphas;
    vector<double> probAlphas;
    vector<double> somatorioCustos;
    vector<int> quantSelecoes;
    double menorCusto;
    double custo;
    double *melhorSolucao = nullptr;
    double *solucao;
    double auxTotalMedia;
    double acumulada;
    float melhorAlpha;
    int countBloco;
    int aleatorio;
    int indiceAlpha;

    melhorSolucao = guloso(salas);
    menorCusto = calculaCusto(salas, melhorSolucao);
    buscaLocalSwap(salas, melhorSolucao);
    melhorAlpha = 0;
    
    for(int i=1; i<=NUMALPHAS; i++)
        alphas.push_back(0.1*i);
    for(int i=0; i<NUMALPHAS; i++) {
        solucao = gulosoRandomizado(salas, alphas[i], seed+i);
        buscaLocalSwap(salas, solucao);
        custo = calculaCusto(salas, solucao);
        somatorioCustos.push_back(custo);
        quantSelecoes.push_back(1);
    }
    
    for(int i=0; i<NUMALPHAS; i++)
        probAlphas.push_back(100.0/NUMALPHAS);

    countBloco = 0;
    for(int i=NUMALPHAS; i<numInteracoes; i++) {
        if(countBloco == tamanhoBloco) {
            countBloco = 0;
            auxTotalMedia = 0;
            for(int j=0; j<alphas.size(); j++)
                auxTotalMedia += 1.0/(somatorioCustos[j]/quantSelecoes[j]);
            for(int j=0; j<alphas.size(); j++)
                probAlphas[j] = ((1.0/(somatorioCustos[j]/quantSelecoes[j])) / auxTotalMedia) * 100;
        }

        aleatorio = rand() % 100;
        acumulada = 0;
        for(int j=0; j<alphas.size(); j++) {
            acumulada += probAlphas[j];
            if(aleatorio < acumulada) {
                indiceAlpha = j;
                break;
            }
        }

        solucao = gulosoRandomizado(salas, alphas[indiceAlpha], seed+i);
        buscaLocalSwap(salas, solucao);
        custo = calculaCusto(salas, solucao);
        if(custo < menorCusto) {
            melhorSolucao = solucao;
            menorCusto = custo;
            melhorAlpha = alphas[indiceAlpha];
        } else
            delete[] solucao;

        indiceAlpha++;
    }
    //cout << "MELHOR: " << melhorAlpha << endl;
    return melhorSolucao;
}

double* graspGuloso(vector<Sala> *salas, int numInteracoes) {
    double *solucao;
    double *melhorSolucao;
    double menorCusto;
    double custo;

    melhorSolucao = nullptr;
    menorCusto = INFINITY;
    for(int i=0; i<numInteracoes; i++) {
        solucao = guloso(salas);
        buscaLocalSwap(salas, solucao);
        custo = calculaCusto(salas, solucao);
        if(custo < menorCusto) {
            menorCusto = custo;
            if(melhorSolucao != nullptr)
                delete[] melhorSolucao;
            melhorSolucao = solucao;
        } else
            delete[] solucao;
    }

    return melhorSolucao;
}

double* graspRandomizado(vector<Sala> *salas, int numInteracoes, float alpha) {
    double *solucao;
    double *melhorSolucao;
    double menorCusto;
    double custo;

    melhorSolucao = nullptr;
    menorCusto = INFINITY;
    for(int i=0; i<numInteracoes; i++) {
        solucao = gulosoRandomizado(salas, alpha, -1);//i);
        buscaLocalSwap(salas, solucao);
        custo = calculaCusto(salas, solucao);
        if(custo < menorCusto) {
            menorCusto = custo;
            if(melhorSolucao != nullptr)
                delete[] melhorSolucao;
            melhorSolucao = solucao;
        } else
            delete[] solucao;
    }

    return melhorSolucao;
}

double* graspReativo(vector<Sala> *salas, int numIteracoesReativo, int tamanhoBloco, int seed) {
    double *solucao;

    srand(seed);
    solucao = gulosoReativo(salas, numIteracoesReativo, tamanhoBloco, seed);
    return solucao;
}

bool taNoVetor(vector<int>* vetor, int valor) {
    for(int i=0; i<vetor->size(); i++)
        if(vetor->at(i) == valor)
            return true;
    return false;
}

double* ACO(Grafo* grafo, vector<Sala> *salas, int numIteracoes, int numAgentes, float taxaReducao, float alpha, float beta, int seed) {
    int numSalas;
    int posicoesIniciais[numAgentes];
    int melhorIteracao;
    double aleatorio;
    double somaProb;
    double somaFeromonio;
    double distancia;
    double tamanhoSup;
    double tamanhoInf;
    double* melhorCorredor;
    vector<double*> corredor;
    vector<vector<int>*> solucoes;
    vector<int>* melhorSolucao;
    vector<int> verificaOcorrencia;
    vector<float> prob;
    //Grafo* grafo;
    
    melhorSolucao = nullptr;
    melhorCorredor = nullptr;
    numSalas = salas->size();
    //grafo = new Grafo(numSalas);

    srand(seed);

    /*for(int i=0; i<numSalas; i++)
        for(int j=i+1; j<numSalas; j++) {
            grafo->addAresta(i, j, 0.001);
        }*/

    for(int i=0; i<numSalas; i++)
        for(int j=i+1; j<numSalas; j++) {
            grafo->buscaAresta(i, j)->setPeso(0.00001);
        }

    
    
    for(int i=0; i<numIteracoes; i++) {
        solucoes.clear();
        corredor.clear();
        //verificaOcorrencia.clear();

        for(int i=0; i<numAgentes; i++) { 
            posicoesIniciais[i] = rand() % numSalas;
            //while(taNoVetor(&verificaOcorrencia, posicoesIniciais[i])) {
                //posicoesIniciais[i] = (posicoesIniciais[i]+1) % numSalas;
            //}
            //verificaOcorrencia.push_back(posicoesIniciais[i]);
        }
        
        for(int j=0; j<numAgentes; j++) {
            tamanhoInf = 0;
            corredor.push_back(new double[numSalas]);
            solucoes.push_back(new vector<int>);
            
            solucoes[j]->push_back(posicoesIniciais[j]);
            corredor[j][posicoesIniciais[j]] = salas->at(posicoesIniciais[j]).largura/2.0;
            tamanhoSup = salas->at(posicoesIniciais[j]).largura;

            for(int s=1; s<numSalas; s++) {
                somaProb = 0;
                prob.clear();
                for(int k=0; k<numSalas; k++) {
                    if(taNoVetor(solucoes[j], k))
                        prob.push_back(0);
                    else {
                        double aux;
                        double antes;
                        if(salas->at(solucoes[j]->back()).fluxo[k] == 0) {
                            aux = pow((double)grafo->buscaAresta(solucoes[j]->back(), k)->getPeso(), alpha) * pow(0.00001, beta);
                            antes = pow(0.00001, beta);
                        }
                        else {
                            aux = pow((double)grafo->buscaAresta(solucoes[j]->back(), k)->getPeso(), alpha) * pow(salas->at(solucoes[j]->back()).fluxo[k], beta);
                            antes = pow(salas->at(solucoes[j]->back()).fluxo[k], beta);
                        }
                        //cout << i << ", " << alpha << ", " << beta << ": " << pow((double)grafo->buscaAresta(solucoes[j]->back(), k)->getPeso(), alpha) << " - " << antes << endl;
                        prob.push_back(aux);
                        somaProb += prob.back();
                    }
                }

                for(int k=0; k<numSalas; k++) {
                    prob[k] /= somaProb;
                }

                aleatorio = ((double)(rand() % 100))/100;
                somaProb = 0;
                for(int k=0; k<numSalas; k++) {
                    somaProb += prob[k];

                    if(aleatorio < somaProb) {
                        solucoes[j]->push_back(k);
                        prob.clear();

                        if(tamanhoSup <= tamanhoInf) {
                            corredor[j][k] = tamanhoSup + salas->at(k).largura/2.0;
                            tamanhoSup += salas->at(k).largura;
                        } else {
                            corredor[j][k] = (tamanhoInf + salas->at(k).largura/2.0) * -1;
                            tamanhoInf += salas->at(k).largura;
                        }
                        break;
                    }
                }
            }

            //if(!verificaSolucao(salas, corredor[j]))
            //    cout << "INVALIDA" << endl << endl << endl;
        }
        
        double menorCusto = INFINITY;
        double custo;
        int melhorFormiga;
        for(int a=0; a<numAgentes; a++) {
            custo = calculaCusto(salas, corredor[a]);
            if(custo < menorCusto) {
                melhorFormiga = a;
                menorCusto = custo;
            }
        }

        for(int j=0; j<numSalas; j++)
            for(int k=j+1; k<numSalas; k++)
                grafo->buscaAresta(j, k)->setPeso(grafo->buscaAresta(j, k)->getPeso()*(1.0-taxaReducao));

        for(int a=0; a<numAgentes; a++) {
            custo = menorCusto / calculaCusto(salas, corredor[a]);
            for(int j=1; j<numSalas; j++)
                grafo->buscaAresta(solucoes[a]->at(j-1), solucoes[a]->at(j))->setPeso(grafo->buscaAresta(solucoes[a]->at(j-1), solucoes[a]->at(j))->getPeso() + custo);
        }

        buscaLocalSwap(salas, corredor[melhorFormiga]);
        
        for(int j=0; j<numAgentes; j++) {
            if(melhorSolucao == nullptr) {
                melhorIteracao = i;
                melhorSolucao = solucoes[j];
                melhorCorredor = corredor[j];
            } else if(calculaCusto(salas, corredor[j]) < calculaCusto(salas, melhorCorredor)) {
                delete melhorSolucao;
                delete[] melhorCorredor;
                melhorSolucao = solucoes[j];
                melhorCorredor = corredor[j];
                melhorIteracao = i;
            } else {
                delete solucoes[j];
                delete[] corredor[j];
                solucoes[j] = nullptr;
                corredor[j] = nullptr;
            }
        }
        /*
        double desvioPadrao = 0;
        double mediaFeromonio = 0;
        int count = 0;
        for(int j=0; j<numSalas; j++) {
            for(int k=j+1; k<numSalas; k++) {
                mediaFeromonio += grafo->buscaAresta(j, k)->getPeso();
                count++;
            }
        }
        mediaFeromonio /= count;

        for(int j=0; j<numSalas; j++) {
            for(int k=j+1; k<numSalas; k++) {
                desvioPadrao += pow(mediaFeromonio - grafo->buscaAresta(j, k)->getPeso(), 2);
            }
        }
        desvioPadrao = sqrt(desvioPadrao/count);

        alpha = 1.0/exp(desvioPadrao);
        beta = 1.0/alpha;

        cout << i << ", "  << desvioPadrao << " " << alpha << " " << beta << endl;*/

        alpha += 1.0/numIteracoes;
        beta -= 1.0/numIteracoes;
    }
    //cout << endl;
    //delete grafo;

    cout << "Iteracao: " << melhorIteracao << endl;
    delete melhorSolucao;

    return melhorCorredor;
}

void cenarioUm(string arquivo, double custoObjetivo) {
    int melhorMetodo;
    double *solucao;
    double menorCusto;
    double custo;
    double duracao;
    double duracaoMenorCusto;
    double mediaCusto;
    double mediaTempo;
    clock_t tempo[2];
    vector<Sala> *salas;

    salas = carregaInstancia("../Instancias/" + arquivo);
    if(salas == nullptr)
        return;
    
    mediaCusto = 0;
    mediaTempo = 0;
    menorCusto = INFINITY;
    for(int i=0; i<30; i++) {
        tempo[0] = clock();
        solucao = graspReativo(salas, 50, 5, 55+i);
        tempo[1] = clock();
        duracao = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
        custo = calculaCusto(salas, solucao);

        mediaTempo += duracao;
        mediaCusto += custo;
        if(custo < menorCusto) {
            menorCusto = custo;
            duracaoMenorCusto = duracao;
        }
        delete[] solucao;

        //cout << "GRASP-Reativo(" << i << ")" << ": " << "Custo: " << custo << ", Erro: " << (custo - custoObjetivo)*100/custoObjetivo << ", Tempo: " << duracao << "ms" << endl;
    }

    mediaTempo /= 1000;
    mediaTempo /= 30;
    mediaCusto /= 30;

    cout << "GRASP-Reativo" << ": " << ", Menor Custo: " << menorCusto << ", Erro: " << (menorCusto - custoObjetivo)*100/custoObjetivo << ", Media Custo: " << mediaCusto << ", Media Tempo: " << mediaTempo << "s"  << endl;

    ofstream arqSaida;
    arqSaida.open("resultados.txt",ios::app);

    arqSaida << arquivo << ";";
    arqSaida << menorCusto << ";";
    arqSaida << (menorCusto - custoObjetivo)*100/custoObjetivo << ";";
    arqSaida << mediaCusto << ";";
    arqSaida << mediaTempo << ";";

    arqSaida << endl;
    arqSaida.close();

    for(int i=0; i<salas->size(); i++)
        delete[] salas->at(i).fluxo;

    delete salas;
}

void cenarioACO(string arquivo, double custoObjetivo) {
    int numSeeds;
    int numSalas;
    double desvioPadrao;
    double mediaCusto;
    double mediaTempo;
    double menorCusto;
    double *solucao;
    double custo;
    double duracao;
    clock_t tempo;
    vector<Sala> *salas;
    vector<double> qualidades;;
    Grafo* grafo;

    salas = carregaInstancia("../Instancias/" + arquivo);
    if(salas == nullptr)
        return;

    numSalas = salas->size();
    grafo = new Grafo(numSalas);
    for(int i=0; i<numSalas; i++)
        for(int j=i+1; j<numSalas; j++) {
            grafo->addAresta(i, j, 0.001);
        }
    
    mediaCusto = 0;
    mediaTempo = 0;
    numSeeds = 10;
    menorCusto = INFINITY;
    for(int i=0; i<numSeeds; i++) {
        tempo = clock();
        solucao = ACO(grafo, salas, 500, numSalas/2, 0.05, 0, 1, i);
        duracao = ((double)(clock() - tempo)) / CLOCKS_PER_SEC;
        custo = calculaCusto(salas, solucao);

        mediaTempo += duracao;
        mediaCusto += custo;
        if(custo < menorCusto) {
            menorCusto = custo;
        }

        ofstream arqSaida;
        arqSaida.open("saidasParciais/resultadosACO.txt",ios::app);
        arqSaida << arquivo << "-" << i << ";";
        arqSaida << custo << ";";
        arqSaida << duracao << ";";

        arqSaida << endl;
        arqSaida.close();

        qualidades.push_back(custo);
        delete[] solucao;
    }

    mediaTempo /= numSeeds;
    mediaCusto /= numSeeds;

    desvioPadrao = 0;
    for(int i=0; i<numSeeds; i++)
        desvioPadrao += pow((qualidades[i] - mediaCusto),2);

    desvioPadrao = sqrt(desvioPadrao/numSeeds);


    cout << "ACO" << ": " << "Menor Custo: " << menorCusto << ", Gap: " << (menorCusto - custoObjetivo)*100/custoObjetivo << "%, Media Custo: " << mediaCusto << ", Media Tempo: " << mediaTempo << "s" << ", Desvio Padrão: " << desvioPadrao << endl;

    ofstream arqSaida;
    arqSaida.open("resultadosACO.txt",ios::app);

    arqSaida << arquivo << ";";
    arqSaida << menorCusto << ";";
    arqSaida << (menorCusto - custoObjetivo)*100/custoObjetivo << ";";
    arqSaida << mediaCusto << ";";
    arqSaida << mediaTempo << ";";
    arqSaida << desvioPadrao << ";";

    arqSaida << endl;
    arqSaida.close();

    for(int i=0; i<salas->size(); i++)
        delete[] salas->at(i).fluxo;

    delete grafo;
    delete salas;
}

double teste_apr(string arquivo, double custoObjetivo, int rep){

    int numMetodos = 3;
    int melhorMetodo = 0;
    double *solucao;
    double menorCusto = INFINITY;
    double custo[numMetodos];
    double tempoMetodos[numMetodos];
    double melhorAlpha, scoreAlpha;
    //double melhorAlphaG, scoreAlphaG;
    clock_t tempo[2];
    string metodos[numMetodos] = {"Guloso", "Guloso Randomizado", "GRASP-Guloso"};//, "GRASP-Randomizado"};
    vector<Sala> *salas;

    double* melhorMetodos = new double[numMetodos];
    double* mediaMetodos = new double[numMetodos];
    double* mediaTempoMetodos = new double[numMetodos];
    for(int i = 0; i < 3; i++){
        melhorMetodos[i] = INFINITY;
        mediaMetodos[i] = 0;
        mediaTempoMetodos[i] = 0;
    }

    salas = carregaInstancia("../Instancias/" + arquivo);
        if(salas == nullptr) 
            return -1;

    for(int exec = 0; exec < rep; exec++){

        //Guloso - Construtivo

        tempo[0] = clock();
        solucao = guloso(salas);
        tempo[1] = clock();
        tempoMetodos[0] = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
        custo[0] = calculaCusto(salas, solucao);
        if(custo[0] < menorCusto){   
            menorCusto = custo[0];
            melhorMetodo = 0;
        }
        delete[] solucao;

        //Guloso Randomizado [vários parâmetros]

        melhorAlpha = 0;
        scoreAlpha = INFINITY;
        for(double alpha = 0.1; alpha < 0.95; alpha += 0.1){
            tempo[0] = clock();
            solucao = gulosoRandomizado(salas, alpha, -1);
            tempo[1] = clock();
            tempoMetodos[1] = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
            custo[1] = calculaCusto(salas, solucao);
            if(custo[1] < scoreAlpha){
                scoreAlpha = custo[1];
                melhorAlpha = alpha;
            }
            if(custo[1] < menorCusto) {
                menorCusto = custo[1];
                melhorMetodo = 1;
            }
            delete[] solucao;
        }

        //Grasp Guloso [vários parâmetros ?]

        tempo[0] = clock();
        solucao = graspGuloso(salas, 1);
        tempo[1] = clock();
        tempoMetodos[2] = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
        custo[2] = calculaCusto(salas, solucao);
        if(custo[2] < menorCusto) {
            menorCusto = custo[2];
            melhorMetodo = 2;
        }
        delete[] solucao;

        //Grasp Randomizado [vários parâmetros]
        /*
        melhorAlphaG = 0;
        scoreAlphaG = INFINITY;
        for(double alpha = 0.1; alpha < 1; alpha += 0.1){
            tempo[0] = clock();
            solucao = graspRandomizado(salas, 500, alpha);
            tempo[1] = clock();
            tempoMetodos[3] = (tempo[1] - tempo[0]) * 1000 / CLOCKS_PER_SEC;
            custo[3] = calculaCusto(salas, solucao);
            if(custo[1] < scoreAlphaG){
                scoreAlphaG = custo[1];
                melhorAlphaG = alpha;
            }
            if(custo[3] < menorCusto) {
                menorCusto = custo[3];
                melhorMetodo = 3;
            }
            delete[] solucao;
            cout << "GrR:";
        }
        cout << endl;*/
        // Relatório Exec
        for(int i=0; i<numMetodos; i++){
            mediaMetodos[i] += custo[i];
            mediaTempoMetodos[i] += tempoMetodos[i];
            if(custo[i] < melhorMetodos[i]){
                melhorMetodos[i] = custo[i];
            }
            cout << metodos[i] << ": " << "Custo: " << custo[i] << ", Erro: " << (custo[i] - custoObjetivo)*100/custoObjetivo << ", Tempo: " << tempoMetodos[i] << "ms" << endl;
        }
        //cout << "Melhor: " << menorCusto << ":" << metodos[melhorMetodo] << " : " << melhorAlpha << " : " << melhorAlphaG << endl;

    }

    for(int i=0; i<salas->size(); i++)
        delete[] salas->at(i).fluxo;
    delete salas;

    cout << "Fim da instância:" << endl;
    for(int i=0; i<numMetodos; i++){
        cout << metodos[i] << ": " << " Custo (médio): " << mediaMetodos[i]/rep << "Custo (melhor): " << melhorMetodos[i] << ", Erro: " << (melhorMetodos[i] - custoObjetivo)*100/custoObjetivo << ", Tempo médio: " << mediaTempoMetodos[i] << "ms" << endl;
    }
    cout << "Melhor => " << metodos[melhorMetodo] << " : " << "Custo (melhor): " << melhorMetodos[melhorMetodo] << ", Erro: " << (melhorMetodos[melhorMetodo] - custoObjetivo)*100/custoObjetivo << endl;
    cout << "##############################################" << endl;

    //escrevendo no arquivo de saída

    ofstream arqSaida;
    arqSaida.open("resultados.txt",ios::app);
    arqSaida << arquivo;
    for(int i=0; i<numMetodos; i++){
        arqSaida << "," << melhorMetodos[i];
    }
    for(int i=0; i<numMetodos; i++){
        arqSaida << "," << mediaMetodos[i]/((double)rep);
    }
    for(int i=0; i<numMetodos; i++){
        arqSaida << "," << mediaTempoMetodos[i]/((double)rep);
    }
    arqSaida << endl;
    arqSaida.close();

    delete[] mediaMetodos;
    delete[] mediaTempoMetodos;
    delete[] melhorMetodos;

    return menorCusto;
}

int main() {
    string arquivo;
    
    double melhorCusto[3] = {1374.5,3439.5,3195};
    string tiposMenores[3] = {"S10.txt","S11.txt","A15"};
    for(int i=0; i<3; i++){
        cout << "INSTANCIA: " << tiposMenores[i] << endl;
        cenarioACO(tiposMenores[i], melhorCusto[i]);
    }
    
    double melhor[15] = {12731.0,108016.5,86646.5,68708.0,124017.5,20470.0,208079.0,162196.0,118260.5,332834.0,31972.0,248225.0,85190.0,156666.0,296176.5};
    string tipos[3] = {"QAP_sko42_0*_n","QAP_sko49_0*_n","QAP_sko56_0*_n"};
    for(int i=0; i<3; i++){
        for(int inst_n = 1; inst_n <= 5; inst_n++){
            tipos[i][11] = to_string(inst_n)[0];
            cout << "INSTANCIA: " << tipos[i] << endl;
            cenarioACO(tipos[i], melhor[5*i+inst_n-1]);
        }
    }

    double melhor2[36] = {204103.0, 193218.5, 161486.5, 135133.5, 159114.0, 158976.5, 110507.5, 115321.0, 114165.0, 105948.0,
    109942.5, 91638.0, 445381.5, 407892.5, 416935.5, 313320.0, 320758.5, 363020.5, 273423.0, 269683.5, 295424.0, 227898.0, 246523.0, 
    206494.5, 628859.0, 561170.5, 587824.5, 474046.0, 479963.0, 512453.0, 479692.0, 445032.0, 495059.5, 385430.5, 331192.0, 411206.0};

    string tipos2[12] = {"CAP_n_60_d_30_L_30_*.txt","CAP_n_60_d_30_L_40_*.txt","CAP_n_60_d_30_L_50_*.txt", "CAP_n_60_d_30_L_60_*.txt", "CAP_n_60_d_60_L_30_*.txt", "CAP_n_60_d_60_L_40_*.txt", 
    "CAP_n_60_d_60_L_50_*.txt", "CAP_n_60_d_60_L_60_*.txt", "CAP_n_60_d_90_L_30_*.txt", "CAP_n_60_d_90_L_40_*.txt", "CAP_n_60_d_90_L_50_*.txt", "CAP_n_60_d_90_L_60_*.txt"};
    for(int i=0; i<12; i++){
        for(int inst_n = 1; inst_n <= 3; inst_n++){
            tipos2[i][19] = to_string(inst_n)[0];
            cout << "INSTANCIA: " << tipos2[i] << endl;
            cenarioACO(tipos2[i], melhor2[3*i+inst_n-1]);
        }
    }

    return 0;
}