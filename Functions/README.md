# Funções auxiliares à sinteze de HRTF individualizadas
Funções utilizadas pelas rotinas principais nas fases de pré e pós-processamento de HRTFs 


 ## Funções
### AA_ARI2CIPIC.m 
Realiza a transformação de hrtfs do banco de dados ARI em seu formato de padronização original para o formato de padronização adotado no banco CIPIC (ambos atualmente obsoletos pela adoção do formato SOFA). 


### hor2geo.m
Conversão de coordenadas horizontais polares para geodésicas.


### hro2sph.m
Conversão de coordenadas horizontais polares para esférica.


### hrtf2DTF.m
Transformação de hrtf em dtf pelo método de Middlebrooks (média logaritmica), Esta função é semelhante à SOFAhrtf2dtf() encontrada no SOFA-Toolbox, recomenda-se a utilização da função do Toolbox por possuir menos bugs. 


### itd_synthesis.m
Sintese itd a partir da largura e profundidade da cabeça. Opções para modelo Woodsworth otimizado por Algazi e modelo Khun otimizado por Ramona. 


### LSD.m 
Calculo da distorção espectral logaritmica entre duas hrtfs (recomenda-se o uso da função spec_dist.m, por ser um pouco mais robusta)


### minphase.m 
Calcula a fase mínima a partir da magnitude apenas e aplica a fase de excesso (ITD) sobre a RI [ITA-Toolbox] (o uso do Toolbox torna o processo bem mais lento, recomenda-se o uso da função phase_job.m para realizar a mesma operação).	


### nav2sph.m
Conversão de coordenadas navegacionais para coordenadas esfericas.


### nn_preprocess.m 
Normalização das linhas de uma matriz para terem média zero e variancia um. Permite a aplicação dos mesmos parametros utilizados em uma matriz X1 sobre uma matriz X2, e também permite a reverter o processo de normalização.


### PCA_fun.m
Realiza a analise de componentes principais pelo metodo de decomposição 
dos autovalores a partir da matriz de covariancia. Corresponde ao uso da função built-in do matlab pca() no modo 'eig'. 


### phase_job.m 
Calculo da fase mínima a partir de magnitude apenas e implemtação de atraso temporal correspondente ao ITD.


### pT_nOctaveBands.m 
Subdivisão de espectro em n-bandas de oitava. 


### sofaFit2Grid.m
Conversão de posições fonte para determinadas posições objetivo. Dois método estão disponíveis, 'move': escolha das posições mais proximas entre o grid de entrada e o grid objetivo. 'interp': faz a interpolação de hrirs em posições conhecidas no grid para a projeção em posições objetivo (uso da função interpolateHRTF do matlab audioToolbox).


### Erro_sofaFit2Grid.m
Apesar de não ser uma função, na verdade é uma análise da performance dos diferentes modos disponíveis da função sofaFit2Grid().



### sofaResample
Raz o resample de hrirs SOFA para dada taxa de amostragem objetivo. 


###  spec_dist.m 
Calculo da distorção espectral logaritmica entre duas hrirs.


### sph2hor.m 
Conversão entre coordenadas esféricas para horizontais polares


### sph2nav.m
Conversão entre coordenadas esféricas para navegacionais.

### miinterpolateHRTF
Função de interpolação de HRTFs pelos metodos bilinear e VBAP, adaptação da função interna do matlab r2020a interpolateHRTF()


### natsortfiles / natsort
Ordenar arquivos carregados alfabeticamente ou numericamente
