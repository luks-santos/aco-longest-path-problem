#--------------------------------------------------%
#        Instituto Federal de Minas Gerais         %
#     Departamento de Engenharia e Computacao      %
#       Otimizacao por Colonia de Formigas         %
#            Professor Ciniro Nametala             %
#--------------------------------------------------%

#--------------------------------------------------%
# ROTINAS INICIAIS
#--------------------------------------------------%

#preparacao do ambiente
rm(list = ls())
cat("\014")  # clear console
dir_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir_path)

#bibliotecas usadas
library(igraph)

#funcoes
plotaGrafo <- function(x, vinicial, vfinal, melhor_caminho) {
  graph <- graph_from_data_frame(x, directed = FALSE)
  
  V(graph)$color <- ifelse(V(graph)$name %in% c(vinicial, vfinal), "orange", "white")
  
  edge_ids <- which((paste(get.edgelist(graph)[,1], get.edgelist(graph)[,2]) %in%
                       paste(melhor_caminho[-length(melhor_caminho)], melhor_caminho[-1])) |
                      (paste(get.edgelist(graph)[,2], get.edgelist(graph)[,1]) %in%
                         paste(melhor_caminho[-length(melhor_caminho)], melhor_caminho[-1])))
  
  E(graph)$color <- "black"
  E(graph)$width <- 1
  E(graph)[edge_ids]$color <- "red"
  E(graph)[edge_ids]$width <- 3
  
  plot(graph)
}

#--------------------------------------------------%
# CONFIGURACOES DO ALGORITMO
#--------------------------------------------------%

# Configuracoes do ACO
ferinicial <- 0.01     # feromonio inicial
rho <- 0.07             # taxa de evaporacao
maxit <- 500           # maximo de iteracoes
qtformigas <- 40         # quantidade de formigas

# Configuracoes do experimento
nexec <- 1                # quantidade de experimentos
nomedoexperimento <- 'experimento_A'
exportar_resultados <- TRUE

# Escolha do grafo a ser usado
# A = Grafo com 7 vertices e 11 arestas (exemplo do slide)
#melhor resultado: 21.4
# B = Grafo com 12 vertices e 25 arestas (grafo 1 - trabalho)
#melhor resultado: 33.92
# C = Grafo com 20 vertices e 190 arestas (grafo 2 - trabalho)
#melhor resultado: 168
# D = Grafo com 100 vertices e 8020 arestas (grafo 3 - trabalho)
#melhor resultado: 990

base <- "C"

#--------------------------------------------------%
# PREPARACAO DE DADOS
#--------------------------------------------------%

switch(base,
       "A" = {
         x <- read.table('exemplo_slides.csv', header = T)
         vinicial <- 1
         vfinal <- 4
         nomebase <- "Grafo A - 7 vertices e 11 arestas"
       },
       "B" = {
         x <- read.table('grafo1.csv', header = T)
         vinicial <- 1
         vfinal <- 12
         nomebase <- "Grafo B - 12 vertices e 25 arestas"
       },
       "C" = {
         x <- read.table('grafo2.csv', header = T)
         vinicial <- 1
         vfinal <- 20
         nomebase <- "Grafo C - 20 vertices e 190 arestas"
       },
       "D" = {
         x <- read.table('grafo3.csv', header = T)
         vinicial <- 1
         vfinal <- 100
         nomebase <- "Grafo D - 100 vertices e 8020 arestas"
       }
       )

#separa as colunas em vetores
origens <- x[, 1]     #ponta de origem de um vertice
destinos <- x[, 2]    #ponta de destino de um vertice
pesos <- x[, 3]       #peso associado a aresta

#--------------------------------------------------%
# ACO APLICADO AO CAMINHO MAIS LONGO
#--------------------------------------------------%

# armazena todos os custos obtidos em todos os experimentos
todoscustos <- matrix(0, nrow = nexec, ncol = qtformigas)

# matrix que armazena os custos apenas do primeiro experimento para geracao
# do grafico de convergencia
custos_expcorrente <- matrix(0, nrow = maxit, ncol = qtformigas)

#guarda o melhor caminho
melhor_caminho = ''
custo_melhor_caminho = 0

#PARA CADA EXPERIMENTO FACA
for (nexp in 1:nexec) {
  #inicializa o vetor de feromonios
  feros <- matrix(rep(ferinicial, each = length(origens)), nrow = length(origens), ncol = length(destinos), byrow = TRUE)
  #ENQUANTO NAO ATINGIR O MAXIMO DE ITERACOES FACA
  #inicia as iteracoes
  it <- 1
  while (it <= maxit) {
    #apresenta a % de execucao do algoritmo
    cat("\014")
    cat("Experimento ", nexp, ": Concluido ", (it*100)/maxit, "%\n")
   
    #variavel que armazenara os caminhos e os respectivos custos
    #criados pelas qtformigas formigas
    caminhos <- vector("list", qtformigas)
    custos <- numeric(qtformigas)

    #ENQUANTO NAO FECHAR A ULTIMA FORMIGA FACA
    #caminha com cada uma das formigas na geracao it
    n <- 1
    while (n <= qtformigas) {
      
      #determina o vertice que a formiga comecara o caminho
      vatual <- vinicial
      
      #enquanto a formiga nao gerar um caminho valido repete-se
      fim <- FALSE
      caminho <- c()
      
      #flag que indica se o caminho gerado e valido
      cvalido <- TRUE
      
      #ENQUANTO A FORMIGA NAO PERCORRER UM CAMINHO VALIDO FACA
      while (fim == FALSE) {
        #caminho que a formiga i esta percorrendo
        caminho <- c(caminho, vatual)
        
        #seleciona os indices onde estao os vertices que a formiga esta
        #naquele dado momento
        ind_vatual <- which(origens == vatual)
        
        #retira os vertices ja visitados
        ind_remove <- c()
        for (k in 1:length(ind_vatual)) {
          remove <- (caminho == destinos[ind_vatual[k]])
          if (any(remove) == 1) {
            #indices que precisam ser removidos
            ind_remove <- c(ind_remove, k)
          }
        }
        
        #remove os indices referentes aos vertices ja visitados listados
        if (!is.null(ind_remove) && length(ind_remove) > 0) {
          ind_vatual <- ind_vatual[-ind_remove]
        }
        
        #caso nao sobre nenhum indice no vetor e pq o caminho em
        #questao e invalido, pois todos os vertices de destino ja foram
        #visitado (ciclo), entao volta a formiga para o primeiro vertice e
        #tenta novamente ate gerar um caminho valido.
        if (length(ind_vatual) == 0) {
          n <- n - 1
          fim <- TRUE
          cvalido <- FALSE
        } else {
          #calculo de probabilidade de cada aresta
          probs <- numeric(length(ind_vatual))
          for (k in 1:length(ind_vatual)) {
            probs[k] <- pesos[ind_vatual[k]] * feros[ind_vatual[k]]
          }
          soma <- sum(probs)
          for (k in 1:length(ind_vatual)) {
            if (k == 1) {
              probs[k] <- probs[k] / soma
            } else {
              probs[k] <- probs[k-1] + (probs[k] / soma)
            }
          }
          
          #executa o sorteio de qual caminho sera seguido
          valor <- runif(1)
          ind_caminho <- 1
          while (valor > probs[ind_caminho]) {
            ind_caminho <- ind_caminho + 1
          }
          
          #define novo vertice escolhido
          vatual <- destinos[ind_vatual[ind_caminho]]
          
          #verifica se o vertice escolhido e o ultimo de interesse
          if (vatual == vfinal) {
            fim <- TRUE
            caminho <- c(caminho, vatual)
          }
        }
      }
      
      if (cvalido == TRUE) {
        custo <- 0
        for (k in 1:(length(caminho)-1)) {
          ind_origens <- which(origens == caminho[k])
          
          for (j in 1:length(ind_origens)) {
            if (destinos[ind_origens[j]] == caminho[k+1]) {
              ind_peso <- ind_origens[j]
              break
            }
          }
          custo <- custo + pesos[ind_peso]
        }
        # armazena o caminho e seu respectivo custo
        caminhos[[n]] <- caminho
        custos[n] <- custo
      }
      
      # passa para a proxima formiga (ou volta para anterior se o caminho gerado tiver sido um caminho invalido)
      n <- n + 1
      
    }
  
    #aplica taxa de evaporacao a todos os feromonios
    feros <- (1 - rho) * feros

    #aplica atualizacao de feromonio nos caminhos utilizados
    for (m in 1:length(caminhos)) {
      caminho <- caminhos[[m]]
      #caso o caminho tenha sido usado adiciona ao deltatau o custo
      for (l in 1:(length(caminho)-1)) {
        ind_origens <- which(origens==caminho[l])
        
        for (j in 1:length(ind_origens)) {
          if (destinos[ind_origens[j]]==caminho[l+1]) {
            ind_fero <- ind_origens[j]
            break
          }
        }
        
        feros[ind_fero] <- custos[m] + feros[ind_fero]

      }
    }
    
    #guarda o custo obtido naquela iteracao - util para geracao do
    #grafico de convergencia. Executa apenas na primeira rodada para
    #exportar apenas um grafico de convergencia
    if (nexp == 1) {
      for (z in 1:qtformigas) {
        custos_expcorrente[it,z] <- custos[z]
      }
    }
    
    #verifica e atualiza o melhor caminho encontrado e seu custo
    for (z in 1:qtformigas) {
      if (custo_melhor_caminho < custos[z]) {
        melhor_caminho <- caminhos[[z]]
        custo_melhor_caminho <- custos[z]
      }
    }
    
    #passa para a proxima iteracao
    it <- it + 1
  }
  
  #armazena todos os resultados finais de cada experimento realizado na bateria de testes
  for (z in 1:qtformigas) {
    todoscustos[nexp,z] <- custos[z]
  }
 
  #mostra a melhor solucao do experimento
  cat("\014")
  cat("Quantidade de experimentos: ", nexec, "\n",
      "Melhor caminho: ", melhor_caminho, "\n",
      "Custo do melhor caminho: ", custo_melhor_caminho, "\n",
      "Custo do pior caminho: ", min(todoscustos), "\n",
      "Custo médio: ", mean(todoscustos), "\n")
}

# PLOTA O GRAFO COM O MELHOR CAMINHO DE TODOS ENCONTRADO
plotaGrafo(x, vinicial, vfinal, melhor_caminho)

#grafico de convergencia da primeira execucao
medias <- apply(custos_expcorrente, 1, mean)
plot(medias, type = "l", xlab = "Iteracao", ylab = "Custo do caminho")

#GERANDO RESULTADOS--------------
if (exportar_resultados == TRUE) {
  dir_name <- paste0(nomedoexperimento, '/')
  
  if (dir.exists(dir_name)) {
    unlink(dir_name, recursive = TRUE)
    dir.create(dir_name)
  } else {
    dir.create(dir_name)
  }

  #exportando resultados para arquivo texto
  cat(paste0("CONFIGURACAO DO EXPERIMENTO----------------------------\n",
             "Base estudada: ", nomebase, "\n",
             "Feromonio inicial: ", ferinicial, "\n",
             "Qtde de formigas: ", qtformigas, "\n",
             "Taxa de evaporacao (rho): ", rho, "\n",
             "Qtde de iteracoes: ", maxit, "\n\n",
             "RESULTADOS PARA ", nexec, " EXECUCOES--------------------------\n",
             "Melhor resultado obtido: ", custo_melhor_caminho, "\n",
             "Pior resultado obtido: ", min(todoscustos), "\n",
             "Resultado medio obtido: ", mean(todoscustos), "\n",
             "Desvio padrao total: ", sd(todoscustos)),
      file = paste0(dir_name, '/RESULTADOS.txt'))
}
